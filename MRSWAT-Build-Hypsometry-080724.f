C====================================================================C
C  This program reads in River Mile Cross Section - Version 0522-2024
C====================================================================C
      Program CrossSection 
C====================================================================C
      Parameter (NPTS=800,IMILE=800)
C====================================================================C
CCC   $$$$$$$$$$$$$$  GLOBAL VARIABLE  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C====================================================================C
CCC   NPTS - Number of points for each Cross Section
CCC   IMILE - Number of River Mile 
CCC   X(NPTS,IMILE) - X-Coordinate of each point on Cross Section
CCC   ELV(NPTS,IMILE) - Elevation of each point on Cross Section
CCC   RVM(IMILE) - Store the name of each River Mile
CCC   ISEG(IMILE) - Store number of points for each River Mile
C====================================================================C
      Dimension X(NPTS,IMILE),ELV(NPTS,IMILE)
      Dimension RVM(IMILE),ISEG(IMILE)
      COMMON/RVMile/X,ELV,RVM,RORM
      COMMON/Numseg/ISEG
C====================================================================C
      Character*8 name01(NPTS),name02(NPTS)
      Character*35 string
      Character*80 fname, fnameo
C====================================================================C
      write(*,*) 
      write(*,*) 'MRSWAT Build Hypsometry'
      write(*,*) 
      write(*,*) 'This code reads in a HEC-RAS geometry (.g file)'
      write(*,*) 'and creates hypsometric representations of each'
      write(*,*) 'cross-section.'
      write(*,*) 'That is, cross-sectional area as a function of'
      write(*,*) 'depth at the thalweg'
      write(*,*) 
      write(*,*) 'Enter the HEC-RAS geometry (.g) filename'
      read(*,'(A)') fname 
      open(10,file=TRIM(fname),form='formatted',status='old')
      write(*,*) 'Enter the offset to convert HEC-RAS River Miles to ',
     +'standard River Miles.'
      write(*,*) 'NOTE: value of this as of 062524 is -4.3'
      read(*,*) RORM 
      write(*,*) 'Enter the desired filename for the output hypsometry'
      read(*,'(A)') fnameo 
      open(21,file=TRIM(fnameo),form='formatted',status='unknown')

      write(*,*) 'Running...' 
C====================================================================C
      ict = 0    !!! Initialize number of Cross Sections to 0
C====================================================================C
 100  READ(10,101,END=999) string
 101  Format(32a)
      IF(string(1:4).eq.'Type') then
         Call cn(string(28:33),riverm,6) 
         ict = ict + 1
         RVM(ict) = riverm
      ELSEIF(string(1:4).eq.'#Sta') then
         Call cn(string(12:14),segnum,3) 
         ISEG(ict) = INT(segnum)
         Read(10,201)(name01(I),name02(I),I=1,ISEG(ict)) 
C====================================================================C
CCC  After reading characters --- convert them to real values
C====================================================================C
         DO I=1,ISEG(ict)
         Call cn(name01(I),X1,8) 
         Call cn(name02(I),Z1,8) 
             X(I,ict) = X1
           ELV(I,ict) = Z1
         ENDDO
      ENDIF
C====================================================================C
      GOTO 100
 999  CONTINUE
C====================================================================C
      Print *,''
      Print *,'Total # of Cross Sections = ', ict
      Print *,''
C====================================================================C
CCC  Call comarae to compute area under each segment
C====================================================================C
      Call comarea(ict)
C====================================================================C
 201  Format(10(a8))
C====================================================================C
      stop
      end
C====================================================================C
      subroutine cn(ch,num,icl)
C====================================================================C
C!character to real converter
C====================================================================C
      implicit none
      integer icl
      character(len=icl),intent(in)::ch
      real,intent(out)::num
      integer::istat
C====================================================================C
C     read (ch,"(f8.2)",iostat=istat) num
C====================================================================C
      read (ch,*) num
      if (istat/=0) then
        num=0
      end if
      return
      end
C====================================================================C
      subroutine comarea(ict)
C====================================================================C
      Parameter (NPTS=800,IMILE=800)
C====================================================================C
CCC   $$$$$$$$$$$$$$  GLOBAL VARIABLE  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C====================================================================C
      Dimension X(NPTS,IMILE),ELV(NPTS,IMILE)
      Dimension RVM(IMILE),ISEG(IMILE)
      Dimension XT(NPTS,IMILE),ELVT(NPTS,IMILE)
      Dimension RVMT(IMILE),ISEGT(IMILE)
      COMMON/RVMile/X,ELV,RVM,RORM
      COMMON/Numseg/ISEG
C====================================================================C
CCC   dx and dz are defined for each interval of a River Mile Cross
CCC   Section
C====================================================================C
      Dimension   dx(NPTS,IMILE),   dz(NPTS,IMILE)
      Dimension delx(NPTS,IMILE), delz(NPTS,IMILE)
C====================================================================C
      Dimension ELMAX(IMILE), ELMIN(IMILE), ev(10000,IMILE)
      Dimension XMAX(IMILE), XMIN(IMILE)
      Dimension zz(10000), delXX(10000), delZZ(10000)
      Dimension ipoint(IMILE), area(10000), izct(IMILE)
      Dimension Asum(500,IMILE), ddZ(IMILE), Width(500,IMILE)
      Dimension exarea(IMILE), XW(100)
      Dimension AREAI(100)
      Integer Zint
C====================================================================C
      open(20,file='g21-CRSelev.txt',form='formatted')
C====================================================================C
CCC   Print out X and ELV for checking after reading in
C====================================================================C
      ELMAX(:) = -999.0
      ELMIN(:) = +999.0
      XMAX(:) = -1E+8
      XMIN(:) = 1E+8
      XW(:) = 0.0
      XXMAX = -999.0
      EXMAX = -999.0
      EXMIN = +999.0
C====================================================================C
C RESET RM AND LIMIT RANGE

      DO k = 1, ict
!        i = ict - k + 1
        i = k
        RVMT(i) = RVM(k) + RORM
        ISEGT(i) = ISEG(k)
        do j = 1, ISEG(k)
          XT(j,i) = X(j,k)
          ELVT(j,i) = ELV(j,k)
        end do
      END DO

      i = 0
      DO k = 1, ict
        if (RVMT(k) .GT. -25. .AND. RVMT(k) .LT. 230.0) THEN
          i = i + 1
          RVM(i) = RVMT(k)
          ISEG(i) = ISEGT(k)
          do j = 1, ISEG(i)
            X(j,i) = XT(j,k)
            ELV(j,i) = ELVT(j,k)
          end do
        END IF
      END DO

      ict = i
        
         
      DO k = 1,ict
         Write(20,*)'RVMile=',RVM(k),'NumSeg=', ISEG(k)
         DO i = 1,ISEG(k)
            Write(20,*) X(i,k),ELV(i,k)
            If(X(i,k).gt.XXMAX) XXMAX = X(i,k)
            If(ELV(i,k).gt.EXMAX) EXMAX = ELV(i,k)
            If(ELV(i,k).lt.EXMIN) EXMIN = ELV(i,k)
         ENDDO
         DO i = 1,ISEG(k)
            If(ELV(i,k).gt.ELMAX(k)) ELMAX(k) = ELV(i,k)
            If(ELV(i,k).lt.ELMIN(k)) ELMIN(k) = ELV(i,k)
            If(X(i,k).gt.XMAX(k)) XMAX(k) = X(i,k)
            If(X(i,k).lt.XMIN(k)) XMIN(k) = X(i,k)
         ENDDO         
C====================================================================C
CCC Compute dx and dz for each interval between 2 points
C====================================================================C
         DO i = 1,ISEG(k)-1
            if(X(i+1,k).gt.X(i,k)) then
               dx(i,k) =   X(i+1,k) -   X(i,k)
            endif
            dz(i,k) = ELV(i+1,k) - ELV(i,k)
C====================================================================C
CCC Then subdivide this interval into 5 equal spacing sub-intervals
C====================================================================C
            delx(i,k) = dx(i,k)/5.0
            delz(i,k) = dz(i,k)/5.0
         ENDDO         
C====================================================================C
CCC Initialize number of discrete points for k segment
C====================================================================C
         ivt = 0
         DO i = 1,ISEG(k)-1
            zz(1) = ELV(i,k)  !!! Specify the 1st point of interval
C====================================================================C
            DO j = 2,6
              zz(j) = zz(j-1) + delz(i,k)
            ENDDO
C====================================================================C
            IF(delx(i,k).GT.0) THEN
            DO j = 1,5
              ivt = ivt + 1
              delZZ(ivt) = zz(j+1)-zz(j)
              delXX(ivt) = delx(i,k)
              slope = delZZ(ivt)/delXX(ivt)
              ev(ivt,k) = (zz(j) + (zz(j)+slope*delXX(ivt)))/2.0
            ENDDO
            ENDIF
         ENDDO         
C====================================================================C
         ipoint(k) = ivt
C====================================================================C
      Print 100,'k,ELMAX,ELMIN,ipoint = ',
     &           k,ELMAX(k),ELMIN(k),ipoint(k)
C====================================================================C
      ENDDO     !!! End of DO LOOP ict (number of River Mile Segment)
C====================================================================C
CCC   NOW Computing AREA under each Cross Section Segments
C====================================================================C
      Print *,'    '
      Print *,'XXMAX , EXMAX , EXMIN = ',XXMAX,EXMAX,EXMIN
      Print *,'    '
C====================================================================C
      area(:) = 0.0
      Asum(:,:) = 0.0
C====================================================================C
CCC   Specified 20 increments from Z-min to Z-max 
C====================================================================C
      Zint = 21
      RZMAX = 20.0
C====================================================================C
      DO k = 1,ict
         zmax = ELMAX(k)
         zmin = ELMIN(k)
         zmws = 300.0
         ddZ(k) = (zmax-zmin)/FLOAT(Zint - 2)
          zinc = zmin 
CCC     Print *,'### k, ddZ = ',k,ddZ(k)
C====================================================================C
         DO n = 1, Zint
            rarcum = 0.0
            DO i = 2,ISEG(k)
              IF (ELV(i,k) .lt. zinc) THEN
                IF (ELV(i-1,k) .ge. zinc) THEN 
                   xint = X(i,k) + (X(i-1,k) - X(i,k)) *
     *            (zinc - ELV(i,k)) / (ELV(i-1,k) - ELV(i,k))
                  rarcum = rarcum + (zinc -ELV(i,k)) * 
     +            (X(i,k) - xint) / 2.0 
                ELSE
                  rarcum = rarcum + ((zinc -ELV(i,k)) + 
     +            (zinc -ELV(i-1,k))) * (X(i,k) - X(i-1,k)) / 2.0
                END IF
              ELSE IF (ELV(i,k) .ge. zinc) THEN 
                IF (ELV(i-1,k) .lt. zinc) THEN
                   xint = X(i,k) + (X(i-1,k) - X(i,k)) *
     *            (zinc - ELV(i,k)) / (ELV(i-1,k) - ELV(i,k))
                  rarcum = rarcum + (zinc -ELV(i-1,k)) *
     +            (xint - X(i-1,k)) / 2.0
                END IF
              END IF
            END DO
            Asum(n,k) = rarcum
            zinc = zinc + ddZ(k)
            if (n .eq. Zint - 1) zinc = ZMWS
         END DO

C User numerical derivative to calculate width

         zinc = zmin + DDZ(k) - 1.E-4 
         DO  n = 2, Zint-1
            rarcum = 0.0
            DO i = 2,ISEG(k)
              IF (ELV(i,k) .lt. zinc) THEN
                IF (ELV(i-1,k) .ge. zinc) THEN
                   xint = X(i,k) + (X(i-1,k) - X(i,k)) *
     *            (zinc - ELV(i,k)) / (ELV(i-1,k) - ELV(i,k))
                  rarcum = rarcum + (zinc -ELV(i,k)) *
     +            (X(i,k) - xint) / 2.0
                ELSE
                  rarcum = rarcum + ((zinc -ELV(i,k)) +
     +            (zinc -ELV(i-1,k))) * (X(i,k) - X(i-1,k)) / 2.0
                END IF
              ELSE IF (ELV(i,k) .ge. zinc) THEN
                IF (ELV(i-1,k) .lt. zinc) THEN
                   xint = X(i,k) + (X(i-1,k) - X(i,k)) *
     *            (zinc - ELV(i,k)) / (ELV(i-1,k) - ELV(i,k))
                  rarcum = rarcum + (zinc -ELV(i-1,k)) *
     +            (xint - X(i-1,k)) / 2.0
                END IF
              END IF
            END DO
            Width(n,k) = (Asum(n,k) - rarcum) / (1.E-4)
            zinc = zinc + ddZ(k)
         END DO
         Width(Zint,k) = XMAX(k) - XMIN(k) 
         Width(1,k) = 0.0
      ENDDO
C Convert to metric

      DO k = 1,ict
         ELMIN(k) = ELMIN(k) * 0.3048
         ELMAX(k) = ELMAX(k) * 0.3048
         ddZ(k) = ddZ(k) * 0.3048
         DO n = 1,Zint
           Asum(n,k) = Asum(n,k) * 0.3048  * 0.3048
           Width(n,k) = Width(n,k) * 0.3048
         ENDDO
      ENDDO
      ZMWS = ZMWS * 0.3048


C====================================================================C
CCC      Print *,'### k, Asum = ', k,Asum(Zint,k)
C====================================================================C
CCC      Compute Width for each Cross Section River Mile
C====================================================================C
C====================================================================C
C====================================================================C
C====================================================================C
C      open(21,file='g21-CRSarea.txt',form='formatted')
C====================================================================C
      Write(21,*) 'Number-of-cross-sections Min-RM MAX-RM'
      Write(21,*) ict, RVM(ict), RVM(1)
      DO k = 1,ict
C====================================================================C
       Write(21,*) 'Cross-section RM, Thalweg-elev(m) No-of-values'
       write(21,*) k, RVM(k), ELMIN(k), '21' 
C====================================================================C
C        zinc = ELMIN(k) 
C====================================================================C
         zinc = 0.0 
         Write(21,*) 'Depth(m) CS-Area(m2) Width(m)'
         DO n = 1,Zint
           Write(21,555) zinc,Asum(n,k),Width(n,k)
           zinc = zinc + ddZ(k)
           if (n .eq. Zint - 1) zinc = ZMWS
         ENDDO
      ENDDO     !!! End of DO LOOP ict (number of River Mile Segment)
C====================================================================C
 100  Format(a,I5,2x,f9.4,2x,f9.4,2x,I8)
 444  Format(a,3I6,2x,f9.4,1x,f6.2,1x,f6.2,2x,f9.3)
 555  Format(f10.4,2x,f14.4,2x,f14.4)
C====================================================================C
CCC   This section calculate the exact area between 2 points
C====================================================================C
      open(22,file='g21-EXAarea.txt',form='formatted')
C====================================================================C
      Write(22,*)'  RVmile ------ Exact A ------ Appro A'
      DO k = 1,ict
         exarea(k) = 0.0    !!! Initialized area for this k R-Mile
         DO i = 1,ISEG(k)-1
CCC Calculate Area for rectangle
           rect = (ELMAX(k)-ELV(i,k))*dx(i,k)
           rect = (ZMWS-ELV(i,k))*dx(i,k)
CCC Calculate Area for triangle
           tria = (dz(i,k)*dx(i,k))/2.0
           exarea(k) = exarea(k) + rect + tria
         ENDDO         
C====================================================================C
C      Print *,'### k, exrare = ',k,exarea(k)
C====================================================================C
       Write(22,555) RVM(k),exarea(k),Asum(Zint,k)
C====================================================================C
      ENDDO     !!! End of DO LOOP ict (number of River Mile Segment)
C====================================================================C
      return
      end

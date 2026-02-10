!*********************************************************************
!*********************************************************************
!*********************************************************************
! MRSWAT - Mississippi River Salt Wedge Analysis Tool
! This is a 1d (cros-ssection averaged) 2 layer model, where the interface
! between layers is explicitly resolved.
! It uses a Macormick scheme to solve the 2-layer shallow water equaitons
! This model includes flow diversions, that are given as a percent
! of the local river flow
!
! Also has sediment transport capability.  If requested, MRSWAT will simulate
! transport of silts and clays, including deposition and resuspension
! of deposited sediment.  The code will accumulate sediment deposition
! and will estimate required dredging for a given elevation threshold
!
! revision 121924
!
! Gary L Brown and Phu V Luong
! ERDC Coastal and Hydraulics Laboratory
!*********************************************************************
!*********************************************************************
PROGRAM MRSWAT
IMPLICIT NONE
!*********************************************************************

REAL :: DT, DX, depthT, grav, mu, time, ttime
REAL :: delh, flux, evt, evb
REAL :: tflux, tflxa, tflxb, fluxa, fluxb, fluxi
REAL :: delflux, delv, telva, elva, elvi, telvb, elvb
REAL :: sdcon, tsdcon, tsdca, tsdcb, sdcona, sdconb, sdconi, delsdc
REAL :: RGHT, WDCSA, SWSAL, RMAF, ENTFAC
REAL :: USTAR,  VTAZ 
REAL :: VDIFF,  trst, toeloc, dwvloc, freq, m2ft
REAL :: delcsa, RDUM, RDUM2
REAL :: LP0, LP1, RQT, RQB, UBAF, RGAP
REAL :: wsd, cse, erc, csd, spgrav, por, drelv, drvol, swsed
INTEGER :: NM, N, I, J, II, nstep,irst, ifreq, NLVLSM, NLVLS ,NCC, NCCM, IDUM, NOBSP
INTEGER :: NUMDIV, ISDFLG, NUMDVM


!*********************************************************************
PARAMETER(NM=10000)
PARAMETER(NLVLSM=30)
PARAMETER(NCCM=1000)
PARAMETER(NUMDVM=1000)
PARAMETER(grav=9.806)  ! acceleration of gravity
PARAMETER(DELH=1.E-8) ! specify wet/dry limit in meter
PARAMETER(m2ft=1. / 0.3048) ! metric conversion

!*********************************************************************
REAL ::   X(0:NM) ! coordinates
REAL ::   Z(0:NM) ! elevation of domain bottom 
REAL ::  UT(0:NM),UTS(0:NM),UTSS(0:NM) ! top layer fluid velocity
REAL ::  LT(0:NM),LTS(0:NM),LTSS(0:NM) ! top layer thickness 
REAL ::  AT(0:NM),ATS(0:NM),ATSS(0:NM) ! top layer cross-sectional area
REAL ::  WIF(0:NM) ! river width at layer interface
REAL ::  WWS(0:NM) ! river width at water surface
REAL ::  UB(0:NM),UBS(0:NM),UBSS(0:NM) ! bottom layer fluid velocity
REAL ::  LB(0:NM),LBS(0:NM),LBSS(0:NM) ! bottom layer thickness 
REAL ::  AB(0:NM),ABBS(0:NM),ABSS(0:NM) ! bottom layer cross-sectional are
REAL :: UTL(0:NM),UBL(0:NM)        ! layer fluid velocity
REAL ::  ST(0:NM),STS(0:NM),STSS(0:NM) ! top layer salinity
REAL ::  SB(0:NM),SBS(0:NM),SBSS(0:NM) ! top layer salinity
REAL ::  STL(0:NM),SBL(0:NM)  
REAL ::  SDT(0:NM),SDTS(0:NM),SDTSS(0:NM) ! top layer sediment
REAL ::  SDB(0:NM),SDBS(0:NM),SDBSS(0:NM) ! top layer sediment
REAL ::  SDTL(0:NM),SDBL(0:NM)
REAL ::  SSDIT(0:NM),SSDIB(0:NM),ANDRG(0:NM)
REAL ::  SDEPFT(0:NM),SDEPFB(0:NM),SEROFT(0:NM),SEROFB(0:NM)
REAL ::  TAUGB(0:NM), TAUGT(0:NM), DREDGE(0:NM)
REAL ::  DVT(0:NM),DVB(0:NM),RICH(0:NM)  
REAL ::  CDBED(0:NM),CDSRF(0:NM)
REAL ::  EDIF(0:NM),DDIF(0:NM),ENTR(0:NM),PHI(0:NM)
REAL ::  LTOLD(0:NM),LBOLD(0:NM)
!*********************************************************************
! This section for flipping data & convert meter to mile and feet
!*********************************************************************
REAL :: XMILE(NM), DIVF(0:NM), STU(0:NM), SBU(0:NM), DIVIA(0:NM)
REAL :: SDTU(0:NM), SDBU(0:NM)
REAL ::  XCC(1:NCCM),ZCC(1:NCCM),DH(1:NCCM,1:NLVLSM),ACSEC(1:NCCM,1:NLVLSM),WCSEC(1:NCCM,1:NLVLSM)
REAL :: DIVRM(NUMDVM),DIVFT(NUMDVM),DIVW(NUMDVM),DIVIE(NUMDVM),DIVSAL(NUMDVM),ROBRM(NUMDVM)
REAL :: zws, zt, zb, ute, ube
REAL :: SBANK, HMAX, TOPFRAC,XCMAX, ABAIV
INTEGER :: ILC(0:NM),IOBRM(NUMDVM),IDIVB(NUMDVM,2),IDIV(NUMDVM)
CHARACTER*80 CHARDUM
CHARACTER*80 filein,filerd,filews,fileflx,filegeo,filere,filop
CHARACTER*80 filtsol,filtloc,filic,filsn,fildiv
CHARACTER*80 filsdi,filsds,filsdn,filsdb,filsre,filsdr
!*********************************************************************
! Read in initial conditions and time series boundary conditions.
!*********************************************************************

sdcon = 0.0
do I = 1, NM
  SSDIT(I) = 0.0
  SSDIB(I) = 0.0
  SDEPFT(I) = 0.0
  SDEPFB(I) = 0.0
  SEROFT(I) = 0.0
  SEROFB(I) = 0.0
end do

!*********************************************************************
write(*,*)
write(*,*) '***************************************************'
write(*,*) ' MRSWAT - Mississippi River Salt Wedge Analysis Tool'
write(*,*) ' version 3  - 12-24 '
write(*,*) ' written by Gary L. Brown and Phu V. Luong'
write(*,*) ' Coastal and Hydraulics Laboratory'
write(*,*) ' Engineer Research and Development Center'
write(*,*) '***************************************************'
write(*,*)
write(*,*) 'enter the input file containing the requested run parameters'
read(*,'(A)') filein
open(4,file=TRIM(filein),form='formatted',status='unknown')
!*********************************************************************
Read(4,*)
Read(4,*) DT
Read(4,*)
Read(4,*) DX
DX = DX/m2ft
Read(4,*)
Read(4,*) N
Read(4,*)
Read(4,*) NLVLS
Read(4,*)
Read(4,*) time
time = time*86400.
Read(4,*)
Read(4,*) freq
freq = freq*86400.
ifreq = INT(freq/DT)
Read(4,*)
Read(4,*) RGHT
RGHT = RGHT/100.0
Read(4,*)
Read(4,*) evt
evt = evt/(m2ft*m2ft)
evb = evt
Read(4,*)
Read(4,*) SWSAL
Read(4,*)
Read(4,*) RMAF
Read(4,*) 
Read(4,*) irst
trst = 0.0
if (irst .eq. 1) then
  Read(4,*)
  Read(4,'(A)') filere
  Read(4,*)
  Read(4,*) trst
  trst = trst*86400.
end if
Read(4,*)
Read(4,'(A)') filews
Read(4,*)
read(4,'(A)') fileflx
Read(4,*)
read(4,'(A)') fildiv
Read(4,*)
read(4,'(A)') filegeo
Read(4,*)
read(4,'(A)') filop
Read(4,*)
read(4,'(A)') filic
Read(4,*)
read(4,'(A)') filtsol
Read(4,*)
read(4,'(A)') filsn
Read(4,*)
read(4,'(A)') filtloc
Read(4,*)
read(4,*) ISDFLG
IF (ISDFLG .EQ. 1) THEN
  Read(4,*)
  read(4,'(A)') filsdi
  Read(4,*)
  read(4,'(A)') filsdb
  Read(4,*)
  read(4,'(A)') filsds
  Read(4,*)
  read(4,'(A)') filsdn
  Read(4,*)
  read(4,'(A)') filsdr
  if (irst .eq. 1) then
    Read(4,*)
    Read(4,'(A)') filsre
  end if
END IF

close (4)

spgrav = 2.65
wsd = 0.0
cse = 0.0
erc = 0.0
csd = 0.0
por = 0.3
drelv = 0.0
swsed = 0.0

IF (ISDFLG .EQ. 1) THEN
  open(4,file=TRIM(filsdi),form='formatted',status='unknown')
  Read(4,*)
  Read(4,*) spgrav
  Read(4,*)
  Read(4,*) wsd
  wsd = wsd / 1000.0
  Read(4,*)
  Read(4,*) cse
  Read(4,*)
  Read(4,*) erc
  erc = erc / 1000.0
  Read(4,*)
  Read(4,*) csd
  Read(4,*)
  Read(4,*) por
  Read(4,*)
  Read(4,*) drelv
  drelv = drelv/m2ft
  Read(4,*)
  Read(4,*) swsed
  close (4)
END IF

print *,'     '
print 1100,'You requested the following: '
print *,'     '
print 1110,'Total Number of Days of Simulation ',time/86400.
print 1110,'Frequency of output in days ',freq/86400.
print 1110,'Roughness Height in Centimeters  ',RGHT*100.
print *,'     '
print 1120,'Filename for downstream wsel info -- ', TRIM(filews)
print 1120,'Filename for upstream discharge info -- ', TRIM(fileflx)
print 1120,'Filename for diversion info -- ', TRIM(fildiv)
print 1120,'Filename for geometry/hypsometry data -- ', TRIM(filegeo)
print 1120,'Filename for observation locations -- ', TRIM(filop)
print *,'     '
if (irst .eq. 1) then
  print 1100,'Restart Conditions: '
  print 1120,'Filename for restart conditions -- ', TRIM(filere)
  print 1110,'Restart time in days  ',trst/86400.0
  print *,'     '
end if
print 1120,'Filename for initial conditions output -- ', TRIM(filic)
print 1120,'Filename for full output -- ', TRIM(filtsol)
print 1120,'Filename for output at selected locations -- ', TRIM(filsn)
print 1120,'Filename for toe location output -- ', TRIM(filtloc)
print *,'     '
IF (ISDFLG .EQ. 1) THEN
  print 1100,'SEDIMENT TRANSPORT IS ACTIVE: '
  print *,'     '
  print 1120,'Filename for sediment parameters -- ', TRIM(filsdi)
  print 1120,'Filename for upstream sediment concentrations  -- ', TRIM(filsdb)
  print 1120,'Filename for sediment output at selected locations -- ', &
        TRIM(filsdn)
  print 1120,'Filename for full sediment output -- ', TRIM(filsds)
  print 1120,'Filename for dredge volume output -- ', TRIM(filsdr)
  if (irst .eq. 1) then
    print 1120,'Filename for sediment restart conditions -- ', TRIM(filsre)
  end if
  print *,'     '
END IF
print 1100,'Beginning simulation now '
print *,'     '
1100 Format(a)
1110 Format(a,f14.4)
1120 Format(a,a)
!*********************************************************************
nstep = (time - trst)/DT

! Estiamte eddy viscosity from WEbel and Schwartman and courant criteria
!evt = 3.75E-5*(DX/DT)**3.
!evb = evt


open(103,file=TRIM(filews),form='formatted',status='unknown')
open(104,file=TRIM(fileflx),form='formatted',status='unknown')
open(105,file=TRIM(filegeo), &
         form='formatted',status='unknown')
!*********************************************************************
elvi = 0.0
if (irst .eq. 0) then
  Read(103,*)
  Read(103,*) telva, elva
  Read(103,*) telvb, elvb
  telva = telva * 86400.
  telvb = telvb * 86400.
  elva = elva * .3048
  elvb = elvb * .3048
  elvi = elva
end if

!*********************************************************************
depthT = elvi
!*********************************************************************

! Read in hypsometry for each cross-section
! RM-19 is the max X.  0 is however far upstream neede to make that happen
 XCMAX = FLOAT(N-1)*DX / (5280.0 * 0.3048) - 19.0
 READ(105,*) CHARDUM
 READ(105,*) IDUM, RDUM, RDUM2
 NCC = 0
 DO I = 1, NCCM
   READ(105,*,END=2112) CHARDUM
   READ(105,*,END=2112) II,XCC(I),ZCC(I),IDUM
   XCC(I) = (XCMAX - XCC(I)) * 5280.0 * 0.3048
   READ(105,*,END=2112) CHARDUM
   DO J = 1, NLVLS
     READ(105,*,END=2112) DH(I,J), ACSEC(I,J), WCSEC(I,J)
   END DO 
   NCC = NCC + 1
 END DO 
2112  CONTINUE
 CLOSE(105)

! SET Z and set index for inerpolating to model cross-sections from opbserved 

DO I = 1,N
  XMILE(I) = -19.0 + FLOAT(N-I)*DX / (5280.0 * 0.3048)
  X(I) = FLOAT(I-1)*DX
  DO J = 1,NCC-1
    IF (X(I) .GE. XCC(J) .AND. X(I) .LT. XCC(J+1)) THEN
      LP0 = (XCC(J+1) - X(I)) / (XCC(J+1) - XCC(J))
      LP1 = (X(I) - XCC(J)) / (XCC(J+1) - XCC(J))      
      Z(I) = ZCC(J) * LP0 + ZCC(J+1) * LP1 
      ILC(I) = J
    END IF
  END DO
END DO

IF (ISDFLG .EQ. 1) THEN
  DO I = 1,N
    LB(I) = MAX(DRELV - Z(I),0.0)
    LT(I) = LB(I) + 1.0
  END DO
  CALL COMPUTE_CSA_FROM_THICK(NCCM,NM,1,N,NLVLSM,NLVLS,ILC,X,XCC,DH,ACSEC,WCSEC,LB,LT,ANDRG,AT,WIF,WWS)
END IF    

! define the diversions
 open(106,file=TRIM(fildiv),form='formatted',status='old')
 read(106,*)
 read(106,*) NUMDIV
 IF (NUMDIV .GT. 0) THEN
  DO I = 1, NUMDIV
    read(106,*)
!     I am disabling the ability to read in invert elevations.  If we want to reenable
!     just enable reading in the DIVIE variable    
    READ(106,*) DIVRM(I), DIVFT(I), DIVW(I), DIVIE(I), DIVSAL(I)
    DIVW(I) = DIVW(I) / m2ft
    DIVIE(I) = DIVIE(I) / m2ft
! this converts continuous rate constant so that the diversion volume is correct    
    DIVFT(I) = -LOG(1-DIVFT(I))
  END DO
  close(106)
  DO I = 1, NUMDIV
    rdum = 1.E8
    DO II = 1, N
      IF(abs(XMILE(II) - DIVRM(I)) .lt. rdum) THEN
        IDIV(I) = II
        rdum = abs(XMILE(II) - DIVRM(I))
      END IF
    END DO
  END DO
  DO I = 1, NUMDIV
    IDIVB(I,1) = IDIV(I) - INT(DIVW(I)/(2.*DX))
    IDIVB(I,2) = IDIV(I) + INT(DIVW(I)/(2.*DX))
  END DO
  DO II = 1, N
    STU(II) = 0.0
    SBU(II) = 0.0
    SDTU(II) = 0.0
    SDBU(II) = 0.0
    DIVF(II) = 0.0
    DIVIA(II) = 0.0
    DO I = 1, NUMDIV
      IF (II .GT. IDIVB(I,1) .AND. II .LE. IDIVB(I,2)) THEN
        DIVF(II) = DIVF(II) + DIVFT(I)/FLOAT(IDIVB(I,2) - IDIVB(I,1))
        LB(II) = MAX((DIVIE(I)-Z(II)),0.0)
        SBU(II) = DIVSAL(I)
        LT(II) =  0.0
        CALL COMPUTE_CSA_FROM_THICK(NCCM,NM,II,II,NLVLSM,NLVLS,ILC,X,XCC,DH,ACSEC,WCSEC,LB,LT,DIVIA,AT,WIF,WWS)
      END IF
    END DO
  END DO
END IF

DO I = 1,N
   UT(I) = 0.0
   UB(I) = 0.0
   LT(I) = depthT - Z(I)
   LB(I) = 0.0
   ST(I) = 0.0
   SB(I) = 0.0
   SDT(I) = 0.0
   SDB(I) = 0.0
   IF (I .GE. (N-1)) THEN
     LT(I) = 0.1*(depthT-Z(I))
     LB(I) = 0.9*(depthT-Z(I))
     ST(I) = 0.0
     SB(I) = SWSAL
     SDT(I) = 0.0
     SDB(I) = SWSED
   END IF
ENDDO

CALL COMPUTE_CSA_FROM_THICK(NCCM,NM,1,N,NLVLSM,NLVLS,ILC,X,XCC,DH,ACSEC,WCSEC,LB,LT,AB,AT,WIF,WWS)

UBAF = 1.0
DO II = 1, N
  DO I = 1, NUMDIV
    IF (IDIV(I) .EQ. II) THEN
       UBAF = UBAF * (1. - DIVFT(I))
    END IF
  END DO
END DO  
! MR. SPECIFIC - USE RM 17 as DS reference
!UBAF = UBAF * (AB(1) + AT(1))/(AB(N) + AT(N) + 1.E-8)
UBAF = 1.43 * UBAF




!*********************************************************************
! Specify downstream Boundary
!*********************************************************************
  open(106,file=TRIM(filop),form='formatted',status='old')
  read(106,*) NOBSP
  DO I = 1, NOBSP
    READ(106,*) ROBRM(I)
  END DO
  close(106)
  DO I = 1, NOBSP
    rdum = 1.E8
    DO II = 1, N
      IF(abs(XMILE(II) - ROBRM(I)) .lt. rdum) THEN
        IOBRM(I) = II
        rdum = abs(XMILE(II) - ROBRM(I)) 
      END IF
    END DO
  END DO
if (irst .eq. 1) then
  open(106,file=TRIM(filere),form='formatted',status='old')
  Read(106,*)
  DO I = 1,N
     II = N-I+1
     Read(106,*) rdum,rdum,zws,zt,zb,ST(II),SB(II),ute,ube, &
                 rdum,rdum 
     LB(II) = zt/ m2ft - Z(II)
     IF (LB(II) .LT. 0.0) THEN
       zt = zt - (LB(II) * m2ft)
       LB(II) = 0.0
       SB(II) = 0.0 
     END IF
     LT(II) = (zws - zt) / m2ft
     UB(II) = ube / m2ft
     UT(II) = ute / m2ft 
  END DO

  CALL COMPUTE_CSA_FROM_THICK(NCCM,NM,1,N,NLVLSM,NLVLS,ILC,X,XCC,DH,ACSEC,WCSEC,LB,LT,AB,AT,WIF,WWS)

  close(106)
  IF (ISDFLG .EQ. 1) THEN
    open(106,file=TRIM(filsre),form='formatted',status='old')
    Read(106,*)
    DO I = 1,N
       II = N-I+1
       Read(106,*) rdum,rdum,rdum,rdum,rdum,SDT(II),SDB(II), &
                   rdum,rdum,ute,ube,rdum
       SSDIB(II) = ube / (m2ft*m2ft)
       SSDIT(II) = ute / (m2ft*m2ft)
    END DO
  close(106)  
  END IF
end if

!*********************************************************************
open(107,file=TRIM(filic),form='formatted',status='unknown')
open(110,file=TRIM(filtsol),form='formatted',status='unknown')
open(112,file=TRIM(filsn),form='formatted',status='unknown')
open(113,file=TRIM(filtloc),form='formatted',status='unknown')
open(114,file='diagnostics.txt',form='formatted',status='unknown')

write(110,'(a)') 'Time River-mile Wsel(ft) Layer-Interface-el(ft) Bed-el(ft) ', &
             'Surface-salt(ppt) Bottom-salt(ppt) ', &
             'Surface-vel(ft/sec) ', &
             'Bottom-vel(ft/sec) Surface-Q(cfs) Bottom-Q(cfs) '

write(112,'(a)') 'Time River-mile Wsel(ft) Layer-Interface-el(ft) Bed-el(ft) ', &
             'Surface-salt(ppt) Bottom-salt(ppt) ', &
             'Surface-vel(ft/sec) ', &
             'Bottom-vel(ft/sec) Surface-Q(cfs) Bottom-Q(cfs) '

write(114,'(a)') 'Time River-mile Horizontal-Disperson-Surface-Layer(sqft/sec) ', &
             'Horizontal-Disperson-Bottom-Layer(sqft/sec) ', &
             'Richardson-Number ', &
             'Vertical-Eddy-Visc(sqft/sec) Vertical-Diff-Coeff(sqft/sec) ', &
             'Entrainment-Cofficient ', &
             'Bed-Drag-Coefficient-Surface-Layer ', &
             'Bed-Drag-Coefficient-Bottom-Layer ', &
             'Water-Surface-Width(ft) Layer-Interface-Width(ft) '
IF (ISDFLG .EQ. 1) THEN
  open(115,file=TRIM(filsds),form='formatted',status='unknown')
  open(116,file=TRIM(filsdn),form='formatted',status='unknown')
  open(117,file=TRIM(filsdr),form='formatted',status='unknown')
  write(115,'(a)') 'Time River-mile Wsel(ft) Layer-Interface-el(ft) ', &
             'Bed-el(ft) ', &
             'Surface-sed(ppm) Bottom-sed(ppm) ', &
             'Surface-grain-shear(Pa) Bottom-grain-shear(Pa) ', &
             'Surface-Layer-Deposition(cubic-feet/foot) ', &
             'Bottom-Layer-Deposition(cubic-feet/foot) ', &
             'Dredging-Requirement(cubic-feet/foot) ' 
  write(116,'(a)') 'Time River-mile Wsel(ft) Layer-Interface-el(ft) ', &
             'Bed-el(ft) ', &
             'Surface-sed(ppm) Bottom-sed(ppm) ', &
             'Surface-grain-shear(Pa) Bottom-grain-shear(Pa) ', &
             'Surface-Layer-Deposition(cubic-feet/foot) ', &
             'Bottom-Layer-Deposition(cubic-feet/foot) ', &
             'Dredging-Requirement(cubic-feet/foot) '
END IF
!*********************************************************************

ttime = trst
!*********************************************************************
if (irst .eq. 0) then 
  Write(110,'(a,f10.4)') 'Time = ',ttime/86400.
  If (ISDFLG .EQ. 1) Write(115,'(a,f10.4)') 'Time = ',ttime/86400.
end if

Write(113,'(a)') 'Time(Days) Toe-Location(River-Mile) 0.25-ppt-at-surface(River-Mile)'
IF (ISDFLG .EQ. 1) Write(117,'(a)') 'Time(Days) Dredge-Volume(cubic-yards)'
!*********************************************************************
DO I = N,1,-1     !!! flipping X and convert to mile
!*********************************************************************
  Write(107,'(I8,1x,10(f12.4,1x),2(f10.1,1x))') I, ttime/86400.,XMILE(I),X(I), &
                                  (LT(I)+LB(I)+Z(I))*m2ft, &
                                  (LB(I)+Z(I))*m2ft,Z(I)*m2ft, &
                            ST(I),SB(I),UT(I)*m2ft,UB(I)*m2ft, &
                            AT(I)*m2ft*m2ft,AB(I)*m2ft*m2ft
  if (irst .eq. 0) then 
     RQT = UT(I)*AT(I)*m2ft**3
     RQB = UB(I)*AB(I)*m2ft**3
     Write(110,'(9(f10.4,1x),2(f14.4,1x))') &
                 ttime/86400., XMILE(I), &
                 (LT(I)+LB(I)+Z(I))*m2ft, &
             (LB(I)+Z(I))*m2ft,Z(I)*m2ft, &
       ST(I),SB(I),UT(I)*m2ft,UB(I)*m2ft, &
       RQT,RQB
     IF (ISDFLG .EQ. 1) THEN
       DREDGE(I) = MAX((SSDIB(I) - ANDRG(I)),0.0)      
       Write(115,'(5(f10.4,1x),7(f12.3,1x))') &
                 ttime/86400., XMILE(I), &
                 (LT(I)+LB(I)+Z(I))*m2ft, &
             (LB(I)+Z(I))*m2ft,Z(I)*m2ft, &
       SDT(I),SDB(I),TAUGT(I),TAUGB(I), &
       SSDIT(I)*m2ft*m2ft,SSDIB(I)*m2ft*m2ft, &
       DREDGE(I)*m2ft*m2ft
     END IF
  end if
ENDDO


!*********************************************************************
! MacCormack's Method solver (explicit)
!*********************************************************************

mu = DT/DX

if (irst .eq. 1) then
  Read(104,*)
  Read(104,*) tflxa, fluxa
  Read(104,*,END=1001) tflxb, fluxb
  tflxa = tflxa * 86400.
  tflxb = tflxb * 86400.
  fluxa = fluxa * 0.3048**3
  fluxb = fluxb * 0.3048**3
  fluxi = fluxa
    do while ((ttime .lt. tflxa) .or. (ttime .ge. tflxb)) 
      tflxa = tflxb
      fluxa = fluxb
      fluxi = fluxa
      Read(104,*,END=1001) tflxb, fluxb
      tflxb = tflxb * 86400.
      fluxb = fluxb * 0.3048**3      
    end do

  print *,'     '
  print 1100,'Restart Interpolation from time-series: '
  print *,'     '
  print 1110,'Restart time in days  ',ttime/86400.0
  print 1110,'Lower interpolation time on discharge series  ',tflxa/86400.0
  print 1110,'Lower interpolation discharge  ',fluxa/0.3048**3
  print 1110,'Upper interpolation time on discharge series  ',tflxb/86400.0
  print 1110,'Upper interpolation discharge  ',fluxb/0.3048**3

  Read(103,*)
  Read(103,*) telva, elva 
  Read(103,*,END=1001) telvb, elvb
  telva = telva * 86400.
  telvb = telvb * 86400.
  elva = elva * .3048
  elvb = elvb * .3048
  elvi = elva
    do while ((ttime .lt. telva) .or. (ttime .ge. telvb))
      telva = telvb
      elva = elvb
      elvi = elva
      Read(103,*,END=1001) telvb, elvb
      telvb = telvb * 86400.
      elvb = elvb * .3048
    end do

  print *,'     '
  print 1110,'Lower interpolation time on wsel series  ',telva/86400.0
  print 1110,'Lower interpolation wsel  ',elva/0.3048
  print 1110,'Upper interpolation time on wsel series  ',telvb/86400.0
  print 1110,'Upper interpolation wsel  ',elvb/0.3048
  print *,'     '

  IF (ISDFLG .EQ. 1) THEN

     open(106,file=TRIM(filsdb),form='formatted',status='old')
     Read(106,*)
     Read(106,*) tsdca, sdcona
     Read(106,*,END=1001) tsdcb, sdconb
     tsdca = tsdca * 86400.
     tsdcb = tsdcb * 86400.
     sdconi = sdcona
     do while ((ttime .lt. tsdca) .or. (ttime .ge. tsdcb))
        tsdca = tsdcb
        sdcona = sdconb
        sdconi = sdcona
        Read(106,*,END=1001) tsdcb, sdconb
        tsdcb = tsdcb * 86400.
      end do

    print *,'     '
    print 1110,'Lower interpolation time on sediment series  ',tsdca/86400.0
    print 1110,'Lower interpolation concentration  ',sdcona
    print 1110,'Upper interpolation time on sediment series  ',tsdcb/86400.0
    print 1110,'Upper interpolation concentration  ',sdconb
    print *,'     '    

  END IF        

end if

!*********************************************************************
! Read in time series flux values
!*********************************************************************
if (irst .eq. 0) then
  Read(104,*)
  Read(104,*) tflxa, fluxa
  Read(104,*) tflxb, fluxb
  tflxa = tflxa * 86400.
  fluxa = fluxa * 0.3048**3
  tflxb = tflxb * 86400.
  fluxb = fluxb * 0.3048**3

  fluxi = fluxa
  IF (ISDFLG .EQ. 1) THEN
     open(106,file=TRIM(filsdb),form='formatted',status='old')
     Read(106,*)
     Read(106,*) tsdca, sdcona
     Read(106,*) tsdcb, sdconb
     tsdca = tsdca * 86400.
     tsdcb = tsdcb * 86400.
     sdconi = sdcona
  END IF
end if 
!*********************************************************************

DO II = 1, nstep

    UB(1) = 0.0
    UT(N+1) = UT(N)
    UB(N+1) = UB(N)
    LT(N+1) = LT(N)
    LB(N+1) = LB(N)
    AT(N+1) = AT(N)
    AB(N+1) = AB(N)
    WIF(N+1) = WIF(N) 
    WWS(N+1) = WWS(N)
    ST(N+1) = ST(N)
    SB(N+1) = SB(N)
    SDT(N+1) = SDT(N)
    SDB(N+1) = SDB(N)
    UT(0) = UT(1)
    UB(0) = UB(1)
    LT(0) = LT(1)
    LB(0) = LB(1)
    AT(0) = AT(1)
    AB(0) = AB(1)
    WIF(0) = WIF(1)
    WWS(0) = WWS(1)
    ST(0) = ST(1)
    SB(0) = SB(1)
    SDT(0) = SDT(1)
    SDB(0) = SDB(1)

!*********************************************************************
! Updating time series data
!*********************************************************************
  flux = fluxi
depthT = elvi
IF (ISDFLG .EQ. 1) sdcon = sdconi

!*********************************************************************
! Compute bed shear stress and interface physics
!*********************************************************************

  CALL INTERFACE_PHYSICS(NM,N,GRAV,SWSAL,SWSED,SPGRAV,DELH,RGHT,RMAF,EVT,UB,UT,SB,ST,SDB,SDT,  &
                         AB,AT,WIF,WWS,PHI,CDBED,CDSRF,EDIF,DDIF,ENTR,DVB,DVT,RICH,0)

  IF (ISDFLG .EQ. 1) THEN
    CALL SEDIMENT_PHYSICS(NM,N,DT,WSD,CSD,CSE,ERC,SPGRAV,POR,UB,UT,  &
                             SDB,SDT,AB,AT,WIF,WWS,SDEPFB,SDEPFT, &
                             SEROFB,SEROFT,TAUGB,TAUGT,SSDIB,SSDIT,0)
  END IF

!*********************************************************************
! predictor step in conservative form (1) top and bottom
!*********************************************************************
  DO I = 1,N-1

!*********************************************************************
!Conservative scheme 
!*********************************************************************

    ATS(I) = AT(I) - (mu)*(UT(I+1)*AT(I+1)) &
                   + (mu)*(UT(I)*AT(I))
    ABBS(I) = AB(I) - (mu)*(UB(I+1)*AB(I+1)) &
                   + (mu)*(UB(I)*AB(I))
! the following logic adds diversions.  Note we onyl divert
! from the top layer, and we ratio the diversion to the relative
! thickness of the top layer
    ABAIV = MAX((AB(I) - DIVIA(I)),0.0)
    ATS(I) = ATS(I) - DIVF(I)*UT(I)*AT(I)*AT(I)*mu &
             /(AT(I) + ABAIV + 1.E-8)
    ABBS(I) = ABBS(I) - DIVF(I)*UB(I)*AB(I)*ABAIV*mu &
             /(AT(I) + ABAIV + 1.E-8)
 ! now add entrainment 
    ENTFAC = 1.0
    ATS(I)  = ATS(I)  + ENTR(I)*ABS(UT(I))*DT*WIF(I)*ENTFAC
    ABBS(I) = ABBS(I) - ENTR(I)*ABS(UT(I))*DT*WIF(I)*ENTFAC
    
!*********************************************************************

  ENDDO


!*********************************************************************
! predictor step in conservative (2) top and bottom
!*********************************************************************
  DO I = 2,N-1

!*********************************************************************
!Conservative scheme 
!*********************************************************************
    VDIFF = UT(I)-UB(I)
    WDCSA = (AT(I) + AB(I))/(WWS(I) + 1.E-8) + 1.E-8
    UTL(I) = (UT(I)*AT(I)) &
           - mu*((UT(I+1)**2.)*AT(I+1)-(UT(I)**2.)*AT(I)) &
           - mu*grav*(AT(I))*(LT(I+1)-LT(I)) &
           - mu*grav*(AT(I))*(LB(I+1)-LB(I)) &
           - mu*grav*(AT(I))*(Z(I+1)-Z(I))   &

!*********************************************************************
           + (mu*evt/DX)*(AT(I)*UT(I-1) &
                     -2.0*AT(I)*UT(I) + AT(I)*UT(I+1) &
           + (UT(I) - UT(I+1)) * (AT(I) - AT(I+1)))
           IF(LB(I).GT.delh .AND. LT(I) .GT. delh) THEN
              UTL(I) = UTL(I) - DT*EDIF(I)*2.*VDIFF*WIF(I)/WDCSA
              UTL(I) = UTL(I) - DT*(CDSRF(I)/2.)*(WWS(I) - WIF(I))*UT(I)*abs(UT(I))
           ELSE IF (LB(I) .LE. delh .AND. LT(I) .GT. delh) THEN
              UTL(I) = UTL(I) - DT*(CDSRF(I)/2.)*WWS(I)*UT(I)*abs(UT(I))
           END IF


    STL(I) = ST(I)*AT(I) - mu*(UT(I+1)*AT(I+1)*ST(I+1)) &
           + mu*(UT(I)*AT(I)*ST(I))       &
           + (mu*DVT(I)/DX)*(AT(I)*ST(I-1) &
                     -2.0*AT(I)*ST(I) + AT(I)*ST(I+1) &
           + (ST(I) - ST(I+1)) * (AT(I) - AT(I+1))) &
           + (mu*AT(I)/DX)*((ST(I) - ST(I+1)) * (DVT(I) - DVT(I+1)))
           IF(LB(I).GT.delh .and. LT(I) .gt. delh) THEN
              STL(I) = STL(I) - DT*DDIF(I)*2.*(ST(I)-SB(I))*WIF(I)/WDCSA
           END IF
           ENTFAC = 1.0
           STL(I) = STL(I) +  &
           SB(I)*ENTR(I)*ABS(UT(I))*DT*WIF(I)*ENTFAC 
    IF (ISDFLG .EQ. 1) THEN
      SDTL(I) = SDT(I)*AT(I) - mu*(UT(I+1)*AT(I+1)*SDT(I+1)) &
           + mu*(UT(I)*AT(I)*SDT(I))       &
           + (mu*DVT(I)/DX)*(AT(I)*SDT(I-1) &
                     -2.0*AT(I)*SDT(I) + AT(I)*SDT(I+1) &
           + (SDT(I) - SDT(I+1)) * (AT(I) - AT(I+1))) &
           + (mu*AT(I)/DX)*((SDT(I) - SDT(I+1)) * (DVT(I) - DVT(I+1)))
           IF(LB(I).GT. 100.*WSD*DT .and. LT(I) .gt. 100.*WSD*DT) THEN
              SDTL(I) = SDTL(I) - DT*DDIF(I)*2.*(SDT(I)-SDB(I))*WIF(I)/WDCSA
              SDTL(I) = SDTL(I) - DT*WSD*SDT(I)*WIF(I)
           END IF   
           IF (LT(I) .gt. 100.*WSD*DT) THEN
              SDTL(I) = SDTL(I) + DT*(SEROFT(I) - SDEPFT(I))*(WWS(I) - WIF(I))
           END IF
           ENTFAC = 1.0
           SDTL(I) = SDTL(I) + &
           SDB(I)*ENTR(I)*ABS(UT(I))*DT*WIF(I)*ENTFAC
    END IF
! if the diversion is out, divert at the current concentration
! if the diversion is in (i.e. water is flowing upstream) divert ar
! the concentration of the last time step where it was out (STU)

    ABAIV = MAX((AB(I) - DIVIA(I)),0.0)
    IF (UT(I) .GT. 0.0) THEN
      STL(I) = STL(I) - DIVF(I)*UT(I)*AT(I)*AT(I)*ST(I)*mu & 
      /(AT(I) + ABAIV + 1.E-8)
    ELSE
      STL(I) = STL(I) - DIVF(I)*UT(I)*AT(I)*AT(I)*STU(I)*mu &
      /(AT(I) + ABAIV + 1.E-8)
    END IF

    IF (ISDFLG .EQ. 1) THEN
      IF (UT(I) .GT. 0.0) THEN
        SDTL(I) = SDTL(I) - DIVF(I)*UT(I)*AT(I)*AT(I)*SDT(I)*mu &
        /(AT(I) + ABAIV + 1.E-8)
      ELSE
        SDTL(I) = SDTL(I) - DIVF(I)*UT(I)*AT(I)*AT(I)*SDTU(I)*mu &
       /(AT(I) + ABAIV + 1.E-8)
      END IF
    END IF

!*********************************************************************
    UBL(I) = (UB(I)*AB(I)) &
           - mu*((UB(I+1)**2.)*AB(I+1)-(UB(I)**2.)*AB(I)) &
           - mu*grav*(AB(I))*(LB(I+1)-LB(I)) &
           - mu*grav*PHI(I)*(AB(I))*(LT(I+1)-LT(I)) &
           - mu*grav*(AB(I))*(Z(I+1)-Z(I))   &

!*********************************************************************
           + (mu*evb/DX)*(AB(I)*UB(I-1) &
                     -2.0*AB(I)*UB(I) + AB(I)*UB(I+1) &
           + (UB(I) - UB(I+1)) * (AB(I) - AB(I+1))) &

           - DT*(CDBED(I)/2.)*WIF(I)*UB(I)*abs(UB(I))  
             IF(LT(I).GT. delh .AND. LB(I).GT.delh) THEN
              UBL(I) = UBL(I) + DT*EDIF(I)*2.*VDIFF*WIF(I)/WDCSA
             END IF


    SBL(I) = SB(I)*AB(I) - mu*(UB(I+1)*AB(I+1)*SB(I+1)) &
           + mu*(UB(I)*AB(I)*SB(I))       &
           + (mu*DVB(I)/DX)*(AB(I)*SB(I-1) &
                     -2.0*AB(I)*SB(I) + AB(I)*SB(I+1) &
           + (SB(I) - SB(I+1)) * (AB(I) - AB(I+1))) &
           + (mu*AB(I)/DX)*((SB(I) - SB(I+1)) * (DVB(I) - DVB(I+1)))
           IF(LB(I).GT.delh .and. LT(I) .gt. delh) THEN
              SBL(I) = SBL(I) - DT*DDIF(I)*2.*WIF(I)*(SB(I)-ST(I))/WDCSA
           END IF
           ENTFAC = 1.0
           SBL(I) = SBL(I) -  &
           SB(I)*ENTR(I)*ABS(UT(I))*DT*WIF(I)*ENTFAC          
    IF (ISDFLG .EQ. 1) THEN
      SDBL(I) = SDB(I)*AB(I) - mu*(UB(I+1)*AB(I+1)*SDB(I+1)) &
           + mu*(UB(I)*AB(I)*SDB(I))       &
           + (mu*DVB(I)/DX)*(AB(I)*SDB(I-1) &
                     -2.0*AB(I)*SDB(I) + AB(I)*SDB(I+1) &
           + (SDB(I) - SDB(I+1)) * (AB(I) - AB(I+1))) &
           + (mu*AB(I)/DX)*((SDB(I) - SDB(I+1)) * (DVB(I) - DVB(I+1)))
           IF(LB(I).GT. 100.*WSD*DT .and. LT(I) .gt. 100.*WSD*DT) THEN
              SDBL(I) = SDBL(I) - DT*DDIF(I)*2.*WIF(I)*(SDB(I)-SDT(I))/WDCSA
              SDBL(I) = SDBL(I) + DT*WSD*SDT(I)*WIF(I)
            END IF
           IF(LB(I).GT. 100.*WSD*DT) THEN    
              SDBL(I) = SDBL(I) + DT*(SEROFB(I) - SDEPFB(I))* WIF(I)
           END IF
           ENTFAC = 1.0
           SDBL(I) = SDBL(I) - &
           SDB(I)*ENTR(I)*ABS(UT(I))*DT*WIF(I)*ENTFAC           
    END IF

! if the diversion is out, divert at the current concentration
! if the diversion is in (i.e. water is flowing upstream) divert ar
! the concentration of seawater 

    ABAIV = MAX((AB(I) - DIVIA(I)),0.0)
    IF (UB(I) .GT. 0.0) THEN
      SBL(I) = SBL(I) - DIVF(I)*UB(I)*AB(I)*ABAIV*SB(I)*mu &
      /(AT(I) + ABAIV + 1.E-8)
    ELSE
      SBL(I) = SBL(I) - DIVF(I)*UB(I)*AB(I)*ABAIV*SBU(I)*mu &
      /(AT(I) + ABAIV + 1.E-8)
    END IF

    IF (ISDFLG .EQ. 1) THEN
      IF (UB(I) .GT. 0.0) THEN
        SDBL(I) = SDBL(I) - DIVF(I)*UB(I)*AB(I)*ABAIV*SDB(I)*mu &
        /(AT(I) + ABAIV + 1.E-8)
      ELSE
        SDBL(I) = SDBL(I) - DIVF(I)*UB(I)*AB(I)*ABAIV*SDBU(I)*mu &
        /(AT(I) + ABAIV + 1.E-8)
      END IF
    END IF
!*********************************************************************
  ENDDO
!*********************************************************************
  DO I = 1,N
     delcsa = delh*(1. + WWS(I)) 
     IF(ATS(I).GT.delcsa) THEN
        UTS(I) = UTL(I)/ATS(I) 
        STS(I) = STL(I)/ATS(I)
        SDTS(I) = SDTL(I)/ATS(I)
     ELSE
        UTS(I) = 0.0 
        SBS(I) = (STL(I) + SBL(I)) / (ATS(I) + ABBS(I))
        SDBS(I) = (SDTL(I) + SDBL(I)) / (ATS(I) + ABBS(I))
        ABBS(I) = ABBS(I) + ATS(I)
        ATS(I) = 0.0 
        STS(I) = 0.0
        SDTS(I) = 0.0 
     ENDIF
     delcsa = delh*(1. + WIF(I))
     IF(ABBS(I).GT.delcsa) THEN
        UBS(I) = UBL(I)/ABBS(I) 
        SBS(I) = SBL(I)/ABBS(I)
        SDBS(I) = SDBL(I)/ABBS(I)
     ELSE
        UBS(I) = 0.0
        STS(I) = (STL(I) + SBL(I)) / (ATS(I) + ABBS(I))
        SDTS(I) = (SDTL(I) + SDBL(I)) / (ATS(I) + ABBS(I))
        ATS(I) = ABBS(I) + ATS(I)
        ABBS(I) = 0.0
        SBS(I) = 0.0
        SDBS(I) = 0.0
     ENDIF
  ENDDO

  CALL COMPUTE_THICK_FROM_CSA(NCCM,NM,1,N,NLVLSM,NLVLS,ILC,X,XCC,DH,ACSEC,WCSEC,ABBS,ATS,LBS,LTS,WIF,WWS)

! set the thinkness of the bottom layer at the downstream boudnary
! equal to densimetric critical depth (Keuelegan)

  CALL CRITICAL_DEPTH_BC(N,NM,GRAV,PHI,AB,AT,WWS,UB,UT,UBAF,TOPFRAC)

  LTS(N) = TOPFRAC*(depthT-Z(N))
  LBS(N) = (1. - TOPFRAC)*(depthT-Z(N))
  
  CALL COMPUTE_CSA_FROM_THICK(NCCM,NM,N,N,NLVLSM,NLVLS,ILC,X,XCC,DH,ACSEC,WCSEC,LBS,LTS,ABBS,ATS,WIF,WWS)



!*********************************************************************
  delcsa = delh*(1. + WWS(1))
  if(ATS(1).GT.delcsa) then
     UTS(1) = flux/ATS(1)
  else
     UTS(1) = 0.0
     ATS(1) = 0.0
     LTS(1) = 0.0
  endif
!*********************************************************************

  UTS(N) = 0.0
  UBS(N) = 0.0
  STS(N) = ST(N)
  SBS(N) = SB(N)
  SDTS(N) = SDT(N)
  SDBS(N) = SDB(N)


!*********************************************************************
  UBS(1) = 0.0
  UTS(0) = UTS(1)
  UBS(0) = UBS(1)
  LTS(0) = LTS(1)
  LBS(0) = LBS(1)
  ATS(N+1) = ATS(N)
  ABBS(N+1) = ABBS(N)  
  STS(1) = 0.0
  SBS(1) = SBS(2)
  STS(0) = STS(1)
  SBS(0) = SBS(1)
  SDTS(1) = sdcon
  SDBS(1) = SDBS(2)
  SDTS(0) = SDTS(1)
  SDBS(0) = SDBS(1)
  UTS(N+1) = UTS(N)
  UBS(N+1) = UBS(N)
  LTS(N+1) = LTS(N)
  LBS(N+1) = LBS(N)
  ATS(0) = ATS(1)
  ABBS(0) = ABBS(1)
  STS(N+1) = STS(N)
  SBS(N+1) = SBS(N)
  SDTS(N+1) = SDTS(N)
  SDBS(N+1) = SDBS(N)


  CALL INTERFACE_PHYSICS(NM,N,GRAV,SWSAL,SWSED,SPGRAV,DELH,RGHT,RMAF,EVT,UBS,UTS,SBS,STS,SDBS,SDTS,  &
                         ABBS,ATS,WIF,WWS,PHI,CDBED,CDSRF,EDIF,DDIF,ENTR,DVB,DVT,RICH,1)

  IF (ISDFLG .EQ. 1) THEN
    CALL SEDIMENT_PHYSICS(NM,N,DT,WSD,CSD,CSE,ERC,SPGRAV,POR,UBS,UTS,  &
                             SDBS,SDTS,ABBS,ATS,WIF,WWS,SDEPFB,SDEPFT, &
                             SEROFB,SEROFT,TAUGB,TAUGT,SSDIB,SSDIT,0)
  END IF

!*********************************************************************
! corrector step in conservative form (1) top and bottom
!*********************************************************************
!  DO I = 2,N
  DO I = N,2,-1
!*********************************************************************
!Conservative scheme 
!*********************************************************************
    ATSS(I) = ATS(I) - mu*(UTS(I)*ATS(I)) &
                     + mu*(UTS(I-1)*ATS(I-1)) 
!*********************************************************************
    ABSS(I) = ABBS(I) - mu*(UBS(I)*ABBS(I)) &
                     + mu*(UBS(I-1)*ABBS(I-1)) 
! the following logic adds diversions.  Note we onyl divert
! from the top layer, and we ratio the diversion to the relative
! thickness of the top layer
    ABAIV = MAX((ABBS(I) - DIVIA(I)),0.0)
    ATSS(I) = ATSS(I) - DIVF(I)*UTS(I)*ATS(I)*ATS(I)*mu &
      /(ATS(I) + ABAIV + 1.E-8)
    ABSS(I) = ABSS(I) - DIVF(I)*UBS(I)*ABBS(I)*ABAIV*mu &
      /(ATS(I) + ABAIV + 1.E-8)
! now add entrainment
    IF (I .LT. N) THEN
      ENTFAC = 1.0   
      ATSS(I) = ATSS(I) + ENTR(I)*ABS(UTS(I))*DT*WIF(I)*ENTFAC
      ABSS(I) = ABSS(I) - ENTR(I)*ABS(UTS(I))*DT*WIF(I)*ENTFAC
    END IF
!*********************************************************************

  ENDDO

    LBSS(1) = LBSS(2)
    LTSS(1) = Z(2) + LTSS(2) + LBSS(2) - Z(1) -  LBSS(1)
    LTSS(1) = MAX(LTSS(1),0.0)
    CALL COMPUTE_CSA_FROM_THICK(NCCM,NM,1,1,NLVLSM,NLVLS,ILC,X,XCC,DH,ACSEC,WCSEC,LBSS,LTSS,ABSS,ATSS,WIF,WWS)

!*********************************************************************
! corrector step (2) top and bottom
!*********************************************************************
!  DO I = 2,N-1
  DO I = N-1,2,-1

    VDIFF = UTS(I)-UBS(I)
    WDCSA = (ATS(I) + ABBS(I))/(WWS(I) + 1.E-8) + 1.E-8
!*********************************************************************
!Conservative scheme 
!*********************************************************************
    UTL(I) = (UTS(I)*ATS(I)) &
           - mu*((UTS(I)**2.)*ATS(I)-(UTS(I-1)**2.)*ATS(I-1)) &
           - mu*grav*(ATS(I))*(LTS(I)-LTS(I-1)) &
           - mu*grav*(ATS(I))*(LBS(I)-LBS(I-1)) &
           - mu*grav*(ATS(I))*(Z(I)-Z(I-1))     &

!*********************************************************************
           + (mu*evt/DX)*(ATS(I)*UTS(I-1) &
                   -2.0*ATS(I)*UTS(I) + ATS(I)*UTS(I+1) &
           + (UTS(I-1) - UTS(I)) * (ATS(I-1) - ATS(I)))

           IF(LBS(I).GT.delh .AND. LTS(I) .GT. delh) THEN
              UTL(I) = UTL(I) - DT*EDIF(I)*2.*VDIFF*WIF(I)/WDCSA
              UTL(I) = UTL(I) - DT*(CDSRF(I)/2.)*(WWS(I) - &
              WIF(I))*UTS(I)*abs(UTS(I))
           ELSE IF (LBS(I) .LE. delh .AND. LTS(I) .GT. delh) THEN
              UTL(I) = UTL(I) - DT*(CDSRF(I)/2.)*WWS(I)*  &
              UTS(I)*abs(UTS(I))
           END IF

    STL(I) = STS(I)*ATS(I) - mu*(UTS(I)*ATS(I)*STS(I)) &
           + mu*(UTS(I-1)*ATS(I-1)*STS(I-1))       &
           + (mu*DVT(I)/DX)*(ATS(I)*STS(I-1) &
                     -2.0*ATS(I)*STS(I) + ATS(I)*STS(I+1) &
             + (STS(I-1) - STS(I)) * (ATS(I-1) - ATS(I))) &
             + (mu*ATS(I)/DX)*((STS(I-1) - STS(I)) * (DVT(I-1) - DVT(I)))
           IF(LBS(I).GT.delh .and. LTS(I) .gt. delh) THEN
              STL(I) = STL(I) - DT*DDIF(I)*2.*(STS(I)-SBS(I))*WIF(I)/WDCSA
           END IF
           ENTFAC = 1.0
           STL(I) = STL(I) +  &
           SBS(I)*ENTR(I)*ABS(UTS(I))*DT*WIF(I)*ENTFAC   

    IF (ISDFLG .EQ. 1) THEN
      SDTL(I) = SDTS(I)*ATS(I) - mu*(UTS(I)*ATS(I)*SDTS(I)) &
           + mu*(UTS(I-1)*ATS(I-1)*SDTS(I-1))       &
           + (mu*DVT(I)/DX)*(ATS(I)*SDTS(I-1) &
                     -2.0*ATS(I)*SDTS(I) + ATS(I)*SDTS(I+1) &
             + (SDTS(I-1) - SDTS(I)) * (ATS(I-1) - ATS(I))) &
             + (mu*ATS(I)/DX)*((SDTS(I-1) - SDTS(I)) * (DVT(I-1) - DVT(I)))
           IF(LBS(I).GT. 100.*WSD*DT .and. LTS(I) .gt. 100.*WSD*DT) THEN
              SDTL(I) = SDTL(I) - DT*DDIF(I)*2.*(SDTS(I)-SDBS(I))*WIF(I)/WDCSA
              SDTL(I) = SDTL(I) - DT*WSD*SDTS(I)*WIF(I)
           END IF   
           IF(LTS(I) .gt. 100.*WSD*DT) THEN
              SDTL(I) = SDTL(I) + DT*(SEROFT(I) - SDEPFT(I))*(WWS(I) - WIF(I))
           END IF
           ENTFAC = 1.0
           SDTL(I) = SDTL(I) +  &
           SDBS(I)*ENTR(I)*ABS(UTS(I))*DT*WIF(I)*ENTFAC  
    END IF
! if the diversion is out, divert at the current concentration
! if the diversion is in (i.e. water is flowing upstream) divert ar
! the concentration of the last time step where it was out (STU)
    ABAIV = MAX((ABBS(I) - DIVIA(I)),0.0)
    IF (UTS(I) .GT. 0.0) THEN
      STL(I) = STL(I) - DIVF(I)*UTS(I)*ATS(I)*ATS(I)*STS(I)*mu &
      /(ATS(I) + ABAIV + 1.E-8)
    ELSE
      STL(I) = STL(I) - DIVF(I)*UTS(I)*ATS(I)*ATS(I)*STU(I)*mu &
      /(ATS(I) + ABAIV + 1.E-8)
    END IF

    IF (ISDFLG .EQ. 1) THEN
      IF (UTS(I) .GT. 0.0) THEN
        SDTL(I) = SDTL(I) - DIVF(I)*UTS(I)*ATS(I)*ATS(I)*SDTS(I)*mu &
        /(ATS(I) + ABAIV + 1.E-8)
      ELSE
        SDTL(I) = SDTL(I) - DIVF(I)*UTS(I)*ATS(I)*ATS(I)*SDTU(I)*mu &
        /(ATS(I) + ABAIV + 1.E-8)
      END IF
    END IF

!*********************************************************************
    UBL(I) = (UBS(I)*ABBS(I)) &
           - mu*((UBS(I)**2.)*ABBS(I)-(UBS(I-1)**2.)*ABBS(I-1)) &
           - mu*grav*(ABBS(I))*(LBS(I)-LBS(I-1)) &
           - mu*grav*PHI(I)*(ABBS(I))*(LTS(I)-LTS(I-1)) &
           - mu*grav*(ABBS(I))*(Z(I)-Z(I-1))   &

!*********************************************************************
          + (mu*evb/DX)*(ABBS(I)*UBS(I-1) &
                     -2.0*ABBS(I)*UBS(I) + ABBS(I)*UBS(I+1) &
             + (UBS(I-1) - UBS(I)) * (ABBS(I-1) - ABBS(I))) &
           - DT*(CDBED(I)/2.)*WIF(I)*UBS(I)*abs(UBS(I))  
             IF(LTS(I).GT. delh .AND. LBS(I) .GT. delh) THEN
              UBL(I) = UBL(I) + DT*EDIF(I)*WIF(I)*2.*VDIFF/WDCSA
             END IF

    SBL(I) = SBS(I)*ABBS(I) - mu*(UBS(I)*ABBS(I)*SBS(I)) &
           + mu*(UBS(I-1)*ABBS(I-1)*SBS(I-1))       &
           + (mu*DVB(I)/DX)*(ABBS(I)*SBS(I-1) &
                     -2.0*ABBS(I)*SBS(I) + ABBS(I)*SBS(I+1) &
             + (SBS(I-1) - SBS(I)) * (ABBS(I-1) - ABBS(I))) &
             + (mu*ABBS(I)/DX)*((SBS(I-1) - SBS(I)) * (DVB(I-1) - DVB(I)))
           IF(LBS(I).GT.delh .and. LTS(I) .gt. delh) THEN
              SBL(I) = SBL(I) - DT*DDIF(I)*2.*WIF(I)*(SBS(I)-STS(I))/WDCSA
           END IF
           ENTFAC = 1.0
           SBL(I) = SBL(I) -  &
           SBS(I)*ENTR(I)*ABS(UTS(I))*DT*WIF(I)*ENTFAC  

    IF (ISDFLG .EQ. 1) THEN
      SDBL(I) = SDBS(I)*ABBS(I) - mu*(UBS(I)*ABBS(I)*SDBS(I)) &
           + mu*(UBS(I-1)*ABBS(I-1)*SDBS(I-1))       &
           + (mu*DVB(I)/DX)*(ABBS(I)*SDBS(I-1) &
                     -2.0*ABBS(I)*SDBS(I) + ABBS(I)*SDBS(I+1) &
             + (SDBS(I-1) - SDBS(I)) * (ABBS(I-1) - ABBS(I))) &
             + (mu*ABBS(I)/DX)*((SDBS(I-1) - SDBS(I)) * (DVB(I-1) - DVB(I)))
           IF(LBS(I).GT. 100.*WSD*DT .and. LTS(I) .gt. 100.*WSD*DT) THEN
              SDBL(I) = SDBL(I) - DT*DDIF(I)*2.*WIF(I)*(SDBS(I)-SDTS(I))/WDCSA
              SDBL(I) = SDBL(I) + DT*WSD*SDTS(I)*WIF(I)
           END IF  
           IF(LBS(I).GT. 100.*WSD*DT) THEN 
              SDBL(I) = SDBL(I) + DT*(SEROFB(I) - SDEPFB(I))*WIF(I)
           END IF
           ENTFAC = 1.0
           SDBL(I) = SDBL(I) -  &
           SDBS(I)*ENTR(I)*ABS(UTS(I))*DT*WIF(I)*ENTFAC    
    END IF 
! if the diversion is out, divert at the current concentration
! if the diversion is in (i.e. water is flowing upstream) divert ar
! the concentration of seawater 
    ABAIV = MAX((ABBS(I) - DIVIA(I)),0.0)
    IF (UBS(I) .GT. 0.0) THEN
      SBL(I) = SBL(I) - DIVF(I)*UBS(I)*ABBS(I)*ABAIV*SBS(I)*mu &
      /(ATS(I) + ABAIV + 1.E-8)
    ELSE
      SBL(I) = SBL(I) - DIVF(I)*UBS(I)*ABBS(I)*ABAIV*SBU(I)*mu &
      /(ATS(I) + ABAIV + 1.E-8)
    END IF

    IF (ISDFLG .EQ. 1) THEN
      IF (UBS(I) .GT. 0.0) THEN
        SDBL(I) = SDBL(I) - DIVF(I)*UBS(I)*ABBS(I)*ABAIV*SDBS(I)*mu &
       /(ATS(I) + ABAIV + 1.E-8)
       ELSE
        SDBL(I) = SDBL(I) - DIVF(I)*UBS(I)*ABBS(I)*ABAIV*SDBU(I)*mu &
        /(ATS(I) + ABAIV + 1.E-8)
      END IF
    END IF
!*********************************************************************
  ENDDO
!*********************************************************************
111 Format(a,4f15.8)
!*********************************************************************
  DO I = 1,N
     delcsa = delh*(1. + WWS(I))
     IF(ATSS(I).GT.delcsa) THEN
        UTSS(I) = UTL(I)/ATSS(I) 
        STSS(I) = STL(I)/ATSS(I)
        SDTSS(I) = SDTL(I)/ATSS(I)
     ELSE
        UTSS(I) = 0.0
        SBSS(I) = (STL(I) + SBL(I)) / (ATSS(I) + ABSS(I))
        SDBSS(I) = (SDTL(I) + SDBL(I)) / (ATSS(I) + ABSS(I))
        ABSS(I) = ABSS(I) + ATSS(I)
        ATSS(I) = 0.0
        STSS(I) = 0.0
        SDTSS(I) = 0.0
     ENDIF
     delcsa = delh*(1. + WIF(I))
     IF(ABSS(I).GT.delcsa) THEN
        UBSS(I) = UBL(I)/ABSS(I) 
        SBSS(I) = SBL(I)/ABSS(I)
        SDBSS(I) = SDBL(I)/ABSS(I)
     ELSE
        UBSS(I) = 0.0
        STSS(I) = (STL(I) + SBL(I)) / (ATSS(I) + ABSS(I))
        SDTSS(I) = (SDTL(I) + SDBL(I)) / (ATSS(I) + ABSS(I))
        ATSS(I) = ABSS(I) + ATSS(I)
        ABSS(I) = 0.0
        SBSS(I) = 0.0 
        SDBSS(I) = 0.0
     ENDIF
  ENDDO

  CALL COMPUTE_THICK_FROM_CSA(NCCM,NM,1,N,NLVLSM,NLVLS,ILC,X,XCC,DH,ACSEC,WCSEC,ABSS,ATSS,LBSS,LTSS,WIF,WWS)

!*********************************************************************
  delcsa = delh*(1. + WWS(1))
  if(ATSS(1).GT.delcsa) then
     UTSS(1) = flux/ATSS(1)
  else
     UTSS(1) = 0.0
     ATSS(1) = 0.0
     LTSS(1) = 0.0
  endif
!*********************************************************************
  UBSS(1) = 0.0 
  STSS(1) = 0.0
  SDTSS(1) = sdcon
  SBSS(1) = SBSS(2)
  SDBSS(1) = SDBSS(2)

  IF (ATSS(N) .gt. 1.E-8) THEN
    UTSS(N) = UTSS(N-1) * ATSS(N-1) / ATSS(N)
  ELSE
    UTSS(N) = 0.0
  END IF
  IF (ABSS(N) .gt. 1.E-8) THEN
    UBSS(N) = UBSS(N-1) * ABSS(N-1) / ABSS(N)
  ELSE
    UBSS(N) = 0.0
  END IF



  IF (UTSS(N) .GE. 0 .and. abs(UTSS(N) * ATSS(N)) .gt. 1.E-8) THEN
    STSS(N) = STSS(N-1) * UTSS(N-1) * ATSS(N-1) / (UTSS(N) * ATSS(N))
  ELSE
    STSS(N) = ST(N)
  END IF
  IF (UBSS(N) .GE. 0 .and. abs(UBSS(N) * LBSS(N)) .gt. 1.E-8) THEN
    SBSS(N) = SBSS(N-1) * UBSS(N-1) * ABSS(N-1) / (UBSS(N) * ABSS(N))
  ELSE
    SBSS(N) = SWSAL
  END IF

  IF (UTSS(N) .GE. 0 .and. abs(UTSS(N) * ATSS(N)) .gt. 1.E-8) THEN
    SDTSS(N) = SDTSS(N-1) * UTSS(N-1) * ATSS(N-1) / (UTSS(N) * ATSS(N))
  ELSE
    SDTSS(N) = SDT(N)
  END IF
  IF (UBSS(N) .GE. 0 .and. abs(UBSS(N) * LBSS(N)) .gt. 1.E-8) THEN
    SDBSS(N) = SDBSS(N-1) * UBSS(N-1) * ABSS(N-1) / (UBSS(N) * ABSS(N))
  ELSE
    SDBSS(N) = SWSED
  END IF


!*********************************************************************
! update solution
!*********************************************************************

  DO I = 1,N
    LTOLD(I) = LT(I)
    LBOLD(I) = LB(I)
    AT(I) = 0.5*(AT(I) + ATSS(I))
    AB(I) = 0.5*(AB(I) + ABSS(I))
  END DO

  CALL COMPUTE_THICK_FROM_CSA(NCCM,NM,1,N,NLVLSM,NLVLS,ILC,X,XCC,DH,ACSEC,WCSEC,AB,AT,LB,LT,WIF,WWS)

  DO I = 1, N
    IF (LTOLD(I) .LE. 1.E-8) ST(I) = STSS(I)
    IF (LBOLD(I) .LE. 1.E-8) SB(I) = SBSS(I)
    IF (LTOLD(I) .LE. 1.E-8) SDT(I) = SDTSS(I)
    IF (LBOLD(I) .LE. 1.E-8) SDB(I) = SDBSS(I)
    IF (LT(I) .GT. 1.E-8 .AND.  &
       ABS(0.5*(UT(I)+UTSS(I))) .LT. (1./(2.*mu)) .AND. &
       (PHI(I) .LT. 0.9970 .OR. (LT(I)/(LB(I) + 1.E-8)) .GT. 0.001)) THEN 
      UT(I) = 0.5*(UT(I) + UTSS(I))
      ST(I) = 0.5*(ST(I) + STSS(I))
      SDT(I) = 0.5*(SDT(I) + SDTSS(I))
    ELSE
      AT(I) = MAX(AT(I),0.0)      
      ST(I) = 0.5*(ST(I) + STSS(I))
      SB(I) = 0.5*(SB(I) + SBSS(I))
      SDT(I) = 0.5*(SDT(I) + SDTSS(I))
      SDB(I) = 0.5*(SDB(I) + SDBSS(I))
      SB(I) = (ST(I)*AT(I) + SB(I)*AB(I)) / (AT(I) + AB(I))
      SDB(I) = (SDT(I)*AT(I) + SDB(I)*AB(I)) / (AT(I) + AB(I))
      LB(I) = LB(I) + MAX(LT(I),0.0)
      AB(I) = AB(I) + AT(I)
      LT(I) = 0.0
      AT(I) = 0.0
      UT(I) = 0.0
      ST(I) = 0.0
      SDT(I) = 0.0
    END IF
    IF (LB(I) .GT. 1.E-8 .AND.  &
       ABS(0.5*(UB(I)+UBSS(I))) .LT. (1./(2.*mu)) .AND.  &
       (PHI(I) .LT. 0.9970 .OR. (LB(I)/(LT(I) + 1.E-8)) .GT. 0.001)) THEN 
      UB(I) = 0.5*(UB(I) + UBSS(I))
      SB(I) = 0.5*(SB(I) + SBSS(I))
      SDB(I) = 0.5*(SDB(I) + SDBSS(I))
    ELSE
      AB(I) = MAX(AB(I),0.0)      
      ST(I) = 0.5*(ST(I) + STSS(I))
      SB(I) = 0.5*(SB(I) + SBSS(I))
      SDT(I) = 0.5*(SDT(I) + SDTSS(I))
      SDB(I) = 0.5*(SDB(I) + SDBSS(I))
      ST(I) = (ST(I)*AT(I) + SB(I)*AB(I)) / (AT(I) + AB(I))
      SDT(I) = (SDT(I)*AT(I) + SDB(I)*AB(I)) / (AT(I) + AB(I))
      LT(I) = LT(I) + MAX(LB(I),0.0)
      AT(I) = AT(I) + AB(I)
      LB(I) = 0.0
      AB(I) = 0.0
      UB(I) = 0.0
      SB(I) = 0.0
      SDB(I) = 0.0
    END IF
    IF (ST(I) .LT. 0.0) ST(I) = 0.0
    IF (ST(I) .GT. (SWSAL + 2.0)) ST(I) = SWSAL + 2.0
    IF (SB(I) .LT. 0.0) SB(I) = 0.0
    IF (SB(I) .GT. (SWSAL + 2.0)) SB(I) = SWSAL + 2.0
    IF (SDT(I) .LT. 0.0) SDT(I) = 0.0
    IF (SDB(I) .LT. 0.0) SDB(I) = 0.0

! don't allow any of the flow to go supercritical
!    IF (UT(I) .GT. SQRT(LT(I)*GRAV)) THEN
!      Write(*,'(a,I8,4f12.6)') 'TOP = '    &
!       ,I,ttime/86400.,XMILE(I),UT(I),SQRT(LT(I)*GRAV)
!      UT(I) = SQRT(LT(I)*GRAV)
!    ELSE IF (UT(I) .LT. -SQRT(LT(I)*GRAV)) THEN
!      Write(*,'(a,I8,4f12.6)') 'TOP = '    &
!       ,I,ttime/86400.,XMILE(I),UT(I),SQRT(LT(I)*GRAV)
!      UT(I) = -SQRT(LT(I)*GRAV)
!    END IF 
!    IF (UB(I) .GT. SQRT(LB(I)*GRAV)) THEN
!      Write(*,'(a,I8,4f12.6)') 'BOT = '    &
!       ,I,ttime/86400.,XMILE(I),UB(I),SQRT(LB(I)*GRAV)
!      UB(I) = SQRT(LB(I)*GRAV)
!    ELSE IF (UB(I) .LT. -SQRT(LB(I)*GRAV)) THEN
!      Write(*,'(a,I8,4f12.6)') 'BOT = '    &
!       ,I,ttime/86400.,XMILE(I),UB(I),SQRT(LB(I)*GRAV)
!      UB(I) = -SQRT(LB(I)*GRAV)
!    END IF
    IF (UT(I) .GT. 0.0) SDTU(I) = SDT(I)
    IF (UB(I) .GT. 0.0) SDBU(I) = SDB(I)
  ENDDO

  IF (ISDFLG .EQ. 1) THEN
    CALL SEDIMENT_PHYSICS(NM,N,DT,WSD,CSD,CSE,ERC,SPGRAV,POR,UB,UT,  &
                             SDB,SDT,AB,AT,WIF,WWS,SDEPFB,SDEPFT, &
                             SEROFB,SEROFT,TAUGB,TAUGT,SSDIB,SSDIT,1)
  END IF

!*********************************************************************
!Writing out time series solution
!*********************************************************************
    ttime = II*DT + trst
  IF(MOD(II,ifreq).EQ.0) THEN
    print 1110,'Writing output for time(days) = ',ttime/86400.
    Write(110,'(a,f10.4)') 'Time = ',ttime/86400.
    Write(114,'(a,f10.4)') 'Time = ',ttime/86400.
    IF (ISDFLG .EQ. 1) THEN
      Write(115,'(a,f10.4)') 'Time = ',ttime/86400.
      Write(116,'(a,f10.4)') 'Time = ',ttime/86400.
    END IF
    toeloc = XMILE(N)
    dwvloc = XMILE(N)
    DRVOL = 0.0
    RGAP = 0.0
    DO I=N,1,-1
     RQT = UT(I)*AT(I)*m2ft**3
     RQB = UB(I)*AB(I)*m2ft**3
     Write(110,'(9(f10.4,1x),2(f14.4,1x))') &
                 ttime/86400., XMILE(I), &
                 (LT(I)+LB(I)+Z(I))*m2ft, &
             (LB(I)+Z(I))*m2ft,Z(I)*m2ft, &
       ST(I),SB(I),UT(I)*m2ft,UB(I)*m2ft, &
       RQT,RQB
       IF(LB(I) .GT. 1.0 .AND. SB(I) .GT. 9.0 .AND. RGAP .LT. 2.0) THEN
          toeloc = XMILE(I) 
          RGAP = 0.0
       ELSE 
          IF (I .LT. N-1) THEN    
            RGAP = RGAP + (XMILE(I) - XMILE(I+1))    
          END IF 
       END IF   
       IF(ST(I) .GT. 0.25) THEN
          dwvloc = XMILE(I)
       END IF       

     Write(114,'(5(f10.4,1x),5(f10.6,1x),2(f10.2,1x))') &
                 ttime/86400., XMILE(I), &
                 DVT(I)*m2ft**2,DVB(I)*m2ft**2, &
             RICH(I),EDIF(I)*m2ft**2,DDIF(I)*m2ft**2, &
       ENTR(I),CDSRF(I),CDBED(I),WWS(I)*m2ft,WIF(I)*m2ft
       IF (ISDFLG .EQ. 1) THEN
         DREDGE(I) = MAX((SSDIB(I) - ANDRG(I)),0.0)    
         DRVOL = DRVOL + DREDGE(I)*DX  
         Write(115,'(5(f10.4,1x),7(f12.3,1x))') &
                   ttime/86400., XMILE(I), &
                   (LT(I)+LB(I)+Z(I))*m2ft, &
               (LB(I)+Z(I))*m2ft,Z(I)*m2ft, &
         SDT(I),SDB(I),TAUGT(I),TAUGB(I), &
         SSDIT(I)*m2ft*m2ft,SSDIB(I)*m2ft*m2ft, &
         DREDGE(I)*m2ft*m2ft
       END IF
    ENDDO

    DO J = 1, NOBSP
     I = IOBRM(J)
     RQT = UT(I)*AT(I)*m2ft**3
     RQB = UB(I)*AB(I)*m2ft**3
     Write(112,'(9(f10.4,1x),2(f14.4,1x))') &
                 ttime/86400., XMILE(I), &
                 (LT(I)+LB(I)+Z(I))*m2ft, &
             (LB(I)+Z(I))*m2ft,Z(I)*m2ft, &
       ST(I),SB(I),UT(I)*m2ft,UB(I)*m2ft, &
       RQT,RQB

     IF (ISDFLG .EQ. 1) THEN
       DREDGE(I) = MAX((SSDIB(I) - ANDRG(I)),0.0)
       Write(116,'(5(f10.4,1x),7(f12.3,1x))') &
                 ttime/86400., XMILE(I), &
                 (LT(I)+LB(I)+Z(I))*m2ft, &
             (LB(I)+Z(I))*m2ft,Z(I)*m2ft, &
       SDT(I),SDB(I),TAUGT(I),TAUGB(I), &
       SSDIT(I)*m2ft*m2ft,SSDIB(I)*m2ft*m2ft, &
       DREDGE(I)*m2ft*m2ft
     END IF
    END DO
    Write(113,'(3f10.3)') ttime/86400., toeloc, dwvloc
    IF (ISDFLG .EQ. 1) Write(117,'(f10.3,1x,f12.1)') ttime/86400., &
                       DRVOL*m2ft*m2ft*m2ft/27.0
  ENDIF

!*********************************************************************
! Handle flux BC upstream time series
!*********************************************************************
if((ttime.ge.tflxa).and.(ttime.le.tflxb)) then
delflux = (fluxb-fluxa)*(ttime-tflxa)/(tflxb-tflxa)
fluxi = fluxa+delflux
endif
!*********************************************************************
if((ttime - tflxb) .gt. 1.E-8) then
tflxa = tflxb
fluxa = fluxb
fluxi = fluxa
Read(104,*,END=999) tflxb, fluxb
tflxb = tflxb * 86400.
fluxb = fluxb * 0.3048**3
endif

!*********************************************************************
! Handle sedimet concentration upstream time series
!*********************************************************************
IF (ISDFLG .EQ. 1) THEN
  if((ttime.ge.tsdca).and.(ttime.le.tsdcb)) then
    delsdc = (sdconb-sdcona)*(ttime-tsdca)/(tsdcb-tsdca)
    sdconi = sdcona+delsdc
  endif
!*********************************************************************
  if((ttime - tsdcb) .gt. 1.E-8) then
    tsdca = tsdcb
    sdcona = sdconb
    sdconi = sdcona
    Read(106,*,END=999) tsdcb, sdconb
    tsdcb = tsdcb * 86400.
  endif
END IF

!*********************************************************************
! Handle elevation downstream time series
!*********************************************************************
if((ttime.ge.telva).and.(ttime.le.telvb)) then
delv = (elvb-elva)*(ttime-telva)/(telvb-telva)
elvi = elva+delv
endif
!*********************************************************************
if((ttime - telvb) .gt. 1.E-8) then
telva = telvb
 elva = elvb
 elvi = elva
Read(103,*,END=999) telvb, elvb
telvb = telvb * 86400.
elvb = elvb * .3048
endif
999 Continue
!*********************************************************************
END DO   ! End of time loop
1001 Continue
!*********************************************************************
! Writing out solution at the last time step
!*********************************************************************
close(113)
!*********************************************************************
close(110)
close(112)
close(114)

print *,'     '
print 1100,'Simulation complete '

END PROGRAM MRSWAT


SUBROUTINE COMPUTE_THICK_FROM_CSA(NCCM,NX,ILS,ILE,NLVLSM,NLVLS,ILC,X,XCC,DH,ACSEC,WCSEC,CSAB,CSAS,LAMB,LAMS,WSB,WSS)
IMPLICIT NONE
!*********************************************************************
REAL :: LTOT, CSAT, XLC, LP0, LP1
INTEGER :: I, J, II, NX, NLVLSM, NLVLS, NCCM, ICOUNT, ILS,ILE
REAL ::  CSAB(0:NX),CSAS(0:NX),LAMB(0:NX),LAMS(0:NX),WSB(0:NX),WSS(0:NX)
REAL ::  X(0:NX)
REAL ::  XCC(1:NCCM),DH(1:NCCM,1:NLVLSM),ACSEC(1:NCCM,1:NLVLSM),WCSEC(1:NCCM,1:NLVLSM)
REAL ::  DHL(1:NLVLSM),ACSECL(1:NLVLSM),WCSECL(1:NLVLSM)
INTEGER :: ILC(0:NX)


    DO II = ILS, ILE
      XLC = X(II)
      J = 1
      LP0 = (XCC(ILC(II)+1) - XLC) / (XCC(ILC(II)+1) - XCC(ILC(II)))
      LP1 = (XLC - XCC(ILC(II))) / (XCC(ILC(II)+1) - XCC(ILC(II)))
      DHL(J) =    DH(ILC(II),J)*LP0 + DH(ILC(II)+1,J)*LP1 
      ACSECL(J) = ACSEC(ILC(II),J)*LP0 + ACSEC(ILC(II)+1,J)*LP1 
      WCSECL(J) = WCSEC(ILC(II),J)*LP0 + WCSEC(ILC(II)+1,J)*LP1 

      ICOUNT = 0
      DO I = 1, NLVLS-1
        J = I + 1
        CSAB(II) = MAX(CSAB(II),0.0)
        DHL(J) =    DH(ILC(II),J)*LP0 + DH(ILC(II)+1,J)*LP1
        ACSECL(J) = ACSEC(ILC(II),J)*LP0 + ACSEC(ILC(II)+1,J)*LP1
        WCSECL(J) = WCSEC(ILC(II),J)*LP0 + WCSEC(ILC(II)+1,J)*LP1        
        IF (CSAB(II) .GE. ACSECL(I) .AND. CSAB(II) .LT. ACSECL(I+1)) THEN
          LAMB(II) = DHL(I) + (DHL(I+1) - DHL(I)) &
         *(CSAB(II) - ACSECL(I)) / (ACSECL(I+1) - ACSECL(I)) 
          WSB(II)  = WCSECL(I) + (WCSECL(I+1) - WCSECL(I)) &
        *(CSAB(II) - ACSECL(I)) / (ACSECL(I+1) - ACSECL(I))
         ICOUNT = ICOUNT + 1
        END IF
        CSAT = CSAB(II) + CSAS(II)
        CSAT = MAX(CSAT,0.0)
        IF (CSAT .GE. ACSECL(I) .AND. CSAT .LT. ACSECL(I+1)) THEN
          LTOT = DHL(I) + (DHL(I+1) - DHL(I)) &
        *(CSAT - ACSECL(I)) / (ACSECL(I+1) - ACSECL(I))
          WSS(II)  = WCSECL(I) + (WCSECL(I+1) - WCSECL(I)) &
        *(CSAT - ACSECL(I)) / (ACSECL(I+1) - ACSECL(I))
          LAMS(II) =  LTOT - LAMB(II)
          IF (LAMS(II) .LE. 0.0) LAMS(II) = 0.0
          IF (WSS(II) .LT. WSB(II)) WSS(II) = WSB(II)
          ICOUNT = ICOUNT + 1
        END IF
        IF (ICOUNT .EQ. 2) GO TO 100
      END DO 
 100  CONTINUE
    END DO

END SUBROUTINE COMPUTE_THICK_FROM_CSA

SUBROUTINE COMPUTE_CSA_FROM_THICK(NCCM,NX,ILS,ILE,NLVLSM,NLVLS,ILC,X,XCC,DH,ACSEC,WCSEC,LAMB,LAMS,CSAB,CSAS,WSB,WSS)
IMPLICIT NONE
!*********************************************************************
REAL :: LTOT, CSAT, XLC, LP0, LP1
INTEGER :: I, J, II, NX, NLVLSM, NLVLS, NCCM, ICOUNT, ILS, ILE
REAL ::  CSAB(0:NX),CSAS(0:NX),LAMB(0:NX),LAMS(0:NX),WSB(0:NX),WSS(0:NX)
REAL ::  X(0:NX)
REAL ::  XCC(1:NCCM),DH(1:NCCM,1:NLVLSM),ACSEC(1:NCCM,1:NLVLSM),WCSEC(1:NCCM,1:NLVLSM)
REAL ::  DHL(1:NLVLSM),ACSECL(1:NLVLSM),WCSECL(1:NLVLSM)
INTEGER :: ILC(0:NX)

   DO II = ILS, ILE
      XLC = X(II) 
      J = 1
      LP0 = (XCC(ILC(II)+1) - XLC) / (XCC(ILC(II)+1) - XCC(ILC(II)))
      LP1 = (XLC - XCC(ILC(II))) / (XCC(ILC(II)+1) - XCC(ILC(II)))
      DHL(J) =    DH(ILC(II),J)*LP0 + DH(ILC(II)+1,J)*LP1
      ACSECL(J) = ACSEC(ILC(II),J)*LP0 + ACSEC(ILC(II)+1,J)*LP1
      WCSECL(J) = WCSEC(ILC(II),J)*LP0 + WCSEC(ILC(II)+1,J)*LP1      

      ICOUNT = 0
      DO I = 1, NLVLS-1
        J = I + 1
        LAMB(II) = MAX(LAMB(II),0.0)
        DHL(J) =    DH(ILC(II),J)*LP0 + DH(ILC(II)+1,J)*LP1
        ACSECL(J) = ACSEC(ILC(II),J)*LP0 + ACSEC(ILC(II)+1,J)*LP1
        WCSECL(J) = WCSEC(ILC(II),J)*LP0 + WCSEC(ILC(II)+1,J)*LP1        
        IF (LAMB(II) .GE. DHL(I) .AND. LAMB(II) .LT. DHL(I+1)) THEN
          CSAB(II) = ACSECL(I) + (ACSECL(I+1) - ACSECL(I)) &
         *(LAMB(II) - DHL(I)) / (DHL(I+1) - DHL(I))
          WSB(II)  = WCSECL(I) + (WCSECL(I+1) - WCSECL(I)) &
        *(CSAB(II) - ACSECL(I)) / (ACSECL(I+1) - ACSECL(I))
         ICOUNT = ICOUNT + 1
        END IF
        LTOT = LAMB(II) + LAMS(II)
        LTOT = MAX(LTOT,0.0)
        IF (LTOT .GE. DHL(I) .AND. LTOT .LT. DHL(I+1)) THEN
          CSAT = ACSECL(I) + (ACSECL(I+1) - ACSECL(I)) &
         *(LTOT - DHL(I)) / (DHL(I+1) - DHL(I))
          WSS(II)  = WCSECL(I) + (WCSECL(I+1) - WCSECL(I)) &
        *(CSAT - ACSECL(I)) / (ACSECL(I+1) - ACSECL(I))
          CSAS(II) =  CSAT - CSAB(II)
          IF (CSAS(II) .LE. 0.0) CSAS(II) = 0.0
          IF (WSS(II) .LT. WSB(II)) WSS(II) = WSB(II)
          ICOUNT = ICOUNT + 1
        END IF
        IF (ICOUNT .EQ. 2) GO TO 100
      END DO
 100  CONTINUE
    END DO

END SUBROUTINE COMPUTE_CSA_FROM_THICK  
!*********************************************************************
!*********************************************************************
!*********************************************************************
!*********************************************************************

SUBROUTINE INTERFACE_PHYSICS(NM,N,GRAV,SWSAL,SWSED,SPGRAV,DELH,RGHT,RMAF,EVT,UB,UT,SB,ST,SDB,SDT,  &
                             AB,AT,WIF,WWS,PHI,CDBED,CDSRF,EDIF,DDIF,ENTR,DVB,DVT,RICH,ICALL)
IMPLICIT NONE
!*********************************************************************
REAL :: STTMP, SBTMP, DENT, DENB, DEPFR, ROUGHT, DRH, BTA, BLTFR, ZHEV, VFUS
REAL :: MIXL, DENGR, LMIXL, UDIFF, VDIFF, PHRAT, RMAF, EVT, VVF 
REAL :: GRAV,DELH, RGHT, SDFACT, LBCSA, LTCSA, USTB, USTT, SWSAL, SWSED
REAL :: SDTTMP, SDBTMP, SPGRAV, SPGT, SPGB, KINVISC, CDINT, RETL, USRE, DUOZ
INTEGER :: I, ICALL, NM, N, J
REAL ::  UB(0:NM),UT(0:NM),SB(0:NM),ST(0:NM),SDB(0:NM),SDT(0:NM)
REAL ::  AB(0:NM),AT(0:NM),WIF(0:NM),WWS(0:NM)
REAL ::  CDBED(0:NM),CDSRF(0:NM)
REAL ::  EDIF(0:NM),DDIF(0:NM),ENTR(0:NM),PHI(0:NM)
REAL ::  DVT(0:NM),DVB(0:NM),RICH(0:NM)

!*********************************************************************
! Compute bed shear stress and interface physics
!*********************************************************************

  KINVISC = 8.6E-7
  USRE = UT(1)*(AT(1)/WWS(1))/KINVISC
  IF (USRE .LT. 500.0) USRE = 500.0

  DO I = 1, N

! Limit acceptable range of salinity
    STTMP = ST(I)
    IF (STTMP .GT. SWSAL) STTMP = SWSAL
    IF (STTMP .LT. 0.0) STTMP = 0.0
    SBTMP = SB(I)
    IF (SBTMP .GT. SWSAL) SBTMP = SWSAL
    IF (SBTMP .LT. 0.0) SBTMP = 0.0
    SDTTMP = SDT(I) / 1000000.0
    IF (SDTTMP .GT. 0.3) SDTTMP = 0.3
    IF (SDTTMP .LT. 0.0) SDTTMP = 0.0
    SDBTMP = SDB(I) / 1000000.0
    IF (SDBTMP .GT. 0.3) SDBTMP = 0.3
    IF (SDBTMP .LT. 0.0) SDBTMP = 0.0
! calculate density using linearized equaiton of state
    DENT = 1000. + 0.78*STTMP
    DENB = 1000. + 0.78*SBTMP
    SPGT = SPGRAV*1000./DENT
    SPGB = SPGRAV*1000./DENB
    DENT = DENT * (1. + ((SPGT-1.)*SDTTMP)/(SPGT - (SPGT-1.)*SDTTMP))
    DENB = DENB * (1. + ((SPGB-1.)*SDBTMP)/(SPGB - (SPGB-1.)*SDBTMP))
! calculate density ratio.  Do not allow negative buoyancy
    PHI(I) = DENT/DENB
    IF (PHI(I) .GT. 0.9977) PHI(I) =0.9977 

! calculate bed drag coefficient assuming standard log velocity profile 
! assume flow between plates, so 1/2 bottom thick is the length
! of the boudnary layer

    LBCSA = AB(I)/(WIF(I) + 1.E-8)
    LTCSA = AT(I)/(WWS(I) + 1.E-8)
    DEPFR = MAX(LBCSA/2.,1.E-8)
    ROUGHT = MAX(RGHT,1.E-8)
    IF (ROUGHT .GT. DEPFR) THEN
       DRH = 1. 
    ELSE
      DRH = DEPFR / ROUGHT
    END IF

    SDFACT = 1.0
    IF (ROUGHT .GT. DEPFR) SDFACT = MAX((2.*(DEPFR/ROUGHT) - 1.),0.0)

    BTA = 29.7 * DRH
    CDBED(I) = SDFACT * 2.*((0.4*BTA)/((BTA+1.)*(LOG(BTA+1.)-1.)+1.))**2.
    USTB = SQRT(CDBED(I)/2.) * ABS(UB(I))
! calculate surface drag coefficient assuming standard log velocity profile
! assume open channel flow, so the top thick is the length
! of the boudnary layer

    DEPFR = MAX(LTCSA,1.E-8)
    ROUGHT = MAX(RGHT,1.E-8)
    IF (ROUGHT .GT. DEPFR) THEN
       DRH = 1.
    ELSE
      DRH = DEPFR / ROUGHT
    END IF

    SDFACT = 1.0
    IF (ROUGHT .GT. DEPFR) SDFACT = MAX((2.*(DEPFR/ROUGHT) - 1.),0.0)

    BTA = 29.7 * DRH
    CDSRF(I) = SDFACT * 2.*((0.4*BTA)/((BTA+1.)*(LOG(BTA+1.)-1.)+1.))**2.
    USTT = SQRT(CDSRF(I)/2.) * ABS(UT(I))


    DEPFR = MAX((LTCSA + LBCSA),1.E-8)    
    BLTFR = MAX(LBCSA,1.E-8)
    ZHEV = BLTFR/DEPFR
    VFUS = ABS(UT(I))
    IF (ABS(UB(I)) .GT. VFUS) VFUS = ABS(UB(I))
! calculate eddy viscosity and turbulent diffision from 
! mellor yamada
! limiting mixing length from ratio of vertiical turbulent
! velocity fluxuation to the Brunt Vaisala frequency (Andre et al, 1978).
! the limitation arises from the principle you can't have an amplitude 
! that exceeds the vertical velocity over the buoyant frequency.
! root mean square of the vertical velocity is estimated from 
! Mellor Yamada description of vertical velocity (w = l du/dz) 
! Note, the user is permitted to adjust this, but adjustments shoudl be small because
! salt wedge behavior is very sensitive to this

!    MIXL = 0.4*DEPFR*ZHEV*(1.-ZHEV)
!    DENGR = 4.*GRAV*(1. - PHIP(I))/(DEPFR*(1. + PHIP(I)))
!    IF (DENGR .LT. 0.0) DENGR = 0.0
!    UDIFF = ABS(UT(I) - UB(I))
!    VVF = 2. * MIXL * UDIFF / DEPFR
!    LMIXL = 0.53 * RMAF * VVF/(SQRT(DENGR) + 1.E-8)
!    IF (MIXL .GE. LMIXL) MIXL = LMIXL
!    EDIF(I) = MIXL*MIXL*2.*UDIFF/DEPFR
!    DDIF(I) = 1.26*EDIF(I)

    UDIFF = ABS(UT(I) - UB(I))
    VDIFF = abs(UT(I)-UB(I))
    IF (VDIFF .LT. .01) VDIFF = .01
    RICH(I) = 1.0
    ENTR(I) = 0.0

    SDFACT = 1.0
    DEPFR = MAX(LTCSA,1.E-8)
    IF (ROUGHT .GT. DEPFR) SDFACT = MAX((2.*(DEPFR/ROUGHT) - 1.),0.0)
    DEPFR = MAX(LBCSA/2.,1.E-8)
    IF (ROUGHT .GT. DEPFR) SDFACT = MAX((2.*(DEPFR/ROUGHT) - 1.),0.0)

    IF (LTCSA .GT. DELH .AND. LBCSA .GT. DELH) THEN
      CDINT = 0.001
      RETL = MAX((LTCSA*VDIFF/KINVISC),1.0)
      DO J = 1, 10
        CDINT = 2.*(2.5*LOG(4.5*RETL*SQRT(CDINT/2.)))**(-2.0)
      END DO
      PHRAT = 2.*(DENT-DENB)/(DENT+DENB)
      IF (PHRAT .GT. -0.002337) PHRAT = -0.002337
      RICH(I) = -GRAV*(PHRAT/2.)*(LTCSA+LBCSA)/(VDIFF**2.)
      IF (RICH(I) .LT. 1.0) RICH(I) = 1.0
      EDIF(I) = (1./4.)*CDINT*DEPFR*UDIFF
      DDIF(I) = (1./4.)*CDINT*DEPFR*UDIFF
      EDIF(I) = EDIF(I)/(1.+0.74*RICH(I))
      DDIF(I) = DDIF(I)/(1.+37.*RICH(I)*RICH(I))
! Try Jirka, use complete formula    
      DUOZ = SQRT(500./USRE) + (0.25/(SQRT(RICH(I)**2.+0.0625)))*(1.-SQRT(500./USRE))
      ENTR(I) = 0.038*(1.-RICH(I)/(SQRT(RICH(I)**2.+0.0625))) + (2./USRE)*(1./DUOZ)
      ENTR(I) = ENTR(I) / RMAF
!      ENTR(I) = (0.038/RMAF)*(1.-RICH(I)/(SQRT(RICH(I)**2.+0.0625)))
    ELSE
      EDIF(I) = KINVISC
      DDIF(I) = KINVISC
      ENTR(I) = 0.0
    END IF

    EDIF(I) = MAX(EDIF(I) * SDFACT, KINVISC)
    DDIF(I) = MAX(DDIF(I) * SDFACT, KINVISC)
    ENTR(I) = ENTR(I) * SDFACT


! calculate henderson-sellers adjustment for stratificaiton.  This is 
! based on the Richardson number (RICH).  
!   IF (LTCSA .GT. DELH .AND. LBCSA .GT. DELH) THEN
!     RICH(I) = -GRAV*((DENT-DENB)/(DENT+DENB))*(LTCSA+LBCSA)/(VDIFF**2.)
!     IF (RICH(I) .LT. 0.0) RICH(I) = 0.0 
!     RICH(I) = RICH(I) 
!     EDIF(I) = EDIF(I)/(1.+0.74*RICH(I))
!     DDIF(I) = DDIF(I)/(1.+37.*RICH(I)*RICH(I))
!     EDIF(I) = MAX(EDIF(I),KINEV)
!     DDIF(I) = MAX(DDIF(I),KINEV)
!   ELSE
!     EDIF(I) = KINEV
!     DDIF(I) = KINEV
!   END IF

! estimate horizontal dispersion.  Top layer assumes open channel (Elder)
! bottom layer assumes pipe flow (Taylor)

  IF (ICALL .EQ. 0) THEN
    DVT(I) = MAX((5.93 * LTCSA * USTT),(EVT/10.))
    DVB(I) = MAX((10.06 * LBCSA * USTB),(EVT/10.))
  END IF

  END DO

END SUBROUTINE INTERFACE_PHYSICS

SUBROUTINE SEDIMENT_PHYSICS(NM,N,DT,WSD,CSD,CSE,ERC,SPGRAV,POR,UB,UT,  &
                             SDB,SDT,AB,AT,WIF,WWS,SDEPFB,SDEPFT, &
                             SEROFB,SEROFT,TAUGB,TAUGT,SSDIB,SSDIT,IETS)
IMPLICIT NONE
!*********************************************************************
REAL :: DEPFR, ROUGHT, DRH, BTA, LBCSA, LTCSA
REAL :: WSD,CSD,CSE,ERC,SPGRAV,POR,SDCONV 
REAL :: CDBED, CDSRF, USTB, USTT, DT
INTEGER :: I, NM, N, IETS
REAL ::  UB(0:NM),UT(0:NM),SDB(0:NM),SDT(0:NM)
REAL ::  AB(0:NM),AT(0:NM),WIF(0:NM),WWS(0:NM)
REAL ::  SDEPFT(0:NM),SDEPFB(0:NM),SEROFT(0:NM),SEROFB(0:NM)
REAL ::  SSDIT(0:NM),SSDIB(0:NM),TAUGB(0:NM),TAUGT(0:NM)

!*********************************************************************
! Compute grain shear stress 
!*********************************************************************

! calculate bed drag coefficient assuming standard log velocity profile
! assume flow between plates, so 1/2 bottom thick is the length
! of the boudnary layer

  SDCONV = 1000000.
! sed grain roughness to 1 mm    
  ROUGHT = 0.001

  DO I = 1, N

    LBCSA = AB(I)/(WIF(I) + 1.E-8)
    LTCSA = AT(I)/(WWS(I) + 1.E-8)
    DEPFR = MAX(LBCSA/2.,1.E-8)
    IF (ROUGHT .GT. 29.7 * DEPFR) THEN
       DRH = 1. / 29.7
    ELSE
      DRH = DEPFR / ROUGHT
    END IF

    BTA = 29.7 * DRH

    CDBED = 2.*((0.4*BTA)/((BTA+1.)*(LOG(BTA+1.)-1.)+1.))**2.
    USTB = SQRT(CDBED/2.) * ABS(UB(I))
    TAUGB(I) = 1000.*USTB*USTB
! calculate surface drag coefficient assuming standard log velocity profile
! assume open channel flow, so the top thick is the length
! of the boudnary layer

    DEPFR = MAX(LTCSA,1.E-8)
    IF (ROUGHT .GT. 29.7 * DEPFR) THEN
       DRH = 1. / 29.7
    ELSE
      DRH = DEPFR / ROUGHT
    END IF

    BTA = 29.7 * DRH
    CDSRF = 2.*((0.4*BTA)/((BTA+1.)*(LOG(BTA+1.)-1.)+1.))**2.
    USTT = SQRT(CDSRF/2.) * ABS(UT(I))
    TAUGT(I) = 1000.*USTT*USTT

! calcualte rates of deposition and erosion

    IF (CSD .LT. 1.E-8) THEN
      SDEPFT(I) = 0.0
      SDEPFB(I) = 0.0
    ELSE    
      SDEPFT(I) = WSD*SDT(I)*(1. - MIN((TAUGT(I)/CSD),1.0))
      SDEPFB(I) = WSD*SDB(I)*(1. - MIN((TAUGB(I)/CSD),1.0))
    END IF  

    IF (CSE .LT. 1.E-8) THEN
      SEROFT(I) = 0.0
      SEROFB(I) = 0.0
    ELSE    
      SEROFT(I) = SDCONV*ERC*(MAX((TAUGT(I)/CSE),1.0) - 1.)
      SEROFB(I) = SDCONV*ERC*(MAX((TAUGB(I)/CSE),1.0) - 1.)
    END IF  

    IF (SSDIT(I) .LT. 1.E-3) SEROFT(I) = 0.0
    IF (SSDIB(I) .LT. 1.E-3) SEROFB(I) = 0.0 

  END DO

! calcualte cumulate bed change (volume per unit channel length)

    IF (IETS .EQ. 1) THEN
      DO I = 1, N
        SSDIT(I) = SSDIT(I) + DT*(SDEPFT(I) - SEROFT(I)) * (WWS(I) - WIF(I)) / &
                   (SDCONV * SPGRAV * (1. - POR))
        SSDIB(I) = SSDIB(I) + DT*(SDEPFB(I) - SEROFB(I)) * WIF(I) / &
                   (SDCONV * SPGRAV * (1. - POR))
        SSDIT(I) = MAX(SSDIT(I),0.0)
        SSDIB(I) = MAX(SSDIB(I),0.0) 
      END DO
    END IF


       

END SUBROUTINE SEDIMENT_PHYSICS

SUBROUTINE CRITICAL_DEPTH_BC(I,NM,GRAV,PHI,AB,AT,WWS,UB,UT,UBAF,TOPFRAC)
IMPLICIT NONE
!*********************************************************************
REAL :: PHRAT, DFN, TOPFRAC, UBAF
REAL :: GRAV,WDCSA
INTEGER :: I, NM, N
REAL ::  WWS(0:NM),UB(0:NM),UT(0:NM),SB(0:NM),ST(0:NM)
REAL ::  AB(0:NM),AT(0:NM),PHI(0:NM)

!*********************************************************************
! Use Keuelegan criteria of critical depth to set top later thick at
! OCean boudnary
!*********************************************************************
    PHRAT = 1. - PHI(I)

    WDCSA = (AT(I) + AB(I)) / (WWS(I) + 1.E-8)

    DFN = UBAF*(ABS(UT(1))) / &
          ((PHRAT*GRAV*WDCSA)**(0.5))
    TOPFRAC = (DFN)**(2./3.)
    IF (TOPFRAC .LT. 0.1) TOPFRAC = 0.1
    IF (TOPFRAC .GT. 0.9) TOPFRAC = 0.9

END SUBROUTINE CRITICAL_DEPTH_BC

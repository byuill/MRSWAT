import numpy as np
import sys
import math
import os

class MRSWAT:
    def __init__(self):
        # Parameters
        self.NM = 10000
        self.NLVLSM = 30
        self.NCCM = 1000
        self.NUMDVM = 1000
        self.grav = 9.806  # acceleration of gravity
        self.DELH = 1.E-8  # specify wet/dry limit in meter
        self.m2ft = 1. / 0.3048  # metric conversion

        # Arrays (0:NM) -> size NM + 1
        self.X = np.zeros(self.NM + 1)
        self.Z = np.zeros(self.NM + 1)
        self.UT = np.zeros(self.NM + 1)
        self.UTS = np.zeros(self.NM + 1)
        self.UTSS = np.zeros(self.NM + 1)
        self.LT = np.zeros(self.NM + 1)
        self.LTS = np.zeros(self.NM + 1)
        self.LTSS = np.zeros(self.NM + 1)
        self.AT = np.zeros(self.NM + 1)
        self.ATS = np.zeros(self.NM + 1)
        self.ATSS = np.zeros(self.NM + 1)
        self.WIF = np.zeros(self.NM + 1)
        self.WWS = np.zeros(self.NM + 1)
        self.UB = np.zeros(self.NM + 1)
        self.UBS = np.zeros(self.NM + 1)
        self.UBSS = np.zeros(self.NM + 1)
        self.LB = np.zeros(self.NM + 1)
        self.LBS = np.zeros(self.NM + 1)
        self.LBSS = np.zeros(self.NM + 1)
        self.AB = np.zeros(self.NM + 1)
        self.ABBS = np.zeros(self.NM + 1)
        self.ABSS = np.zeros(self.NM + 1)
        self.UTL = np.zeros(self.NM + 1)
        self.UBL = np.zeros(self.NM + 1)
        self.ST = np.zeros(self.NM + 1)
        self.STS = np.zeros(self.NM + 1)
        self.STSS = np.zeros(self.NM + 1)
        self.SB = np.zeros(self.NM + 1)
        self.SBS = np.zeros(self.NM + 1)
        self.SBSS = np.zeros(self.NM + 1)
        self.STL = np.zeros(self.NM + 1)
        self.SBL = np.zeros(self.NM + 1)
        self.SDT = np.zeros(self.NM + 1)
        self.SDTS = np.zeros(self.NM + 1)
        self.SDTSS = np.zeros(self.NM + 1)
        self.SDB = np.zeros(self.NM + 1)
        self.SDBS = np.zeros(self.NM + 1)
        self.SDBSS = np.zeros(self.NM + 1)
        self.SDTL = np.zeros(self.NM + 1)
        self.SDBL = np.zeros(self.NM + 1)
        self.SSDIT = np.zeros(self.NM + 1)
        self.SSDIB = np.zeros(self.NM + 1)
        self.ANDRG = np.zeros(self.NM + 1)
        self.SDEPFT = np.zeros(self.NM + 1)
        self.SDEPFB = np.zeros(self.NM + 1)
        self.SEROFT = np.zeros(self.NM + 1)
        self.SEROFB = np.zeros(self.NM + 1)
        self.TAUGB = np.zeros(self.NM + 1)
        self.TAUGT = np.zeros(self.NM + 1)
        self.DREDGE = np.zeros(self.NM + 1)
        self.DVT = np.zeros(self.NM + 1)
        self.DVB = np.zeros(self.NM + 1)
        self.RICH = np.zeros(self.NM + 1)
        self.CDBED = np.zeros(self.NM + 1)
        self.CDSRF = np.zeros(self.NM + 1)
        self.EDIF = np.zeros(self.NM + 1)
        self.DDIF = np.zeros(self.NM + 1)
        self.ENTR = np.zeros(self.NM + 1)
        self.PHI = np.zeros(self.NM + 1)
        self.LTOLD = np.zeros(self.NM + 1)
        self.LBOLD = np.zeros(self.NM + 1)

        # Additional arrays
        self.XMILE = np.zeros(self.NM + 1) # Note: Fortran declared XMILE(NM) which is 1:NM usually, but let's make it consistent
        self.DIVF = np.zeros(self.NM + 1)
        self.STU = np.zeros(self.NM + 1)
        self.SBU = np.zeros(self.NM + 1)
        self.DIVIA = np.zeros(self.NM + 1)
        self.SDTU = np.zeros(self.NM + 1)
        self.SDBU = np.zeros(self.NM + 1)

        # Geometry arrays (1:NCCM) -> size NCCM + 1
        self.XCC = np.zeros(self.NCCM + 1)
        self.ZCC = np.zeros(self.NCCM + 1)
        self.DH = np.zeros((self.NCCM + 1, self.NLVLSM + 1))
        self.ACSEC = np.zeros((self.NCCM + 1, self.NLVLSM + 1))
        self.WCSEC = np.zeros((self.NCCM + 1, self.NLVLSM + 1))

        # Diversion arrays (NUMDVM)
        self.DIVRM = np.zeros(self.NUMDVM + 1)
        self.DIVFT = np.zeros(self.NUMDVM + 1)
        self.DIVW = np.zeros(self.NUMDVM + 1)
        self.DIVIE = np.zeros(self.NUMDVM + 1)
        self.DIVSAL = np.zeros(self.NUMDVM + 1)
        self.ROBRM = np.zeros(self.NUMDVM + 1)

        self.ILC = np.zeros(self.NM + 1, dtype=int)
        self.IOBRM = np.zeros(self.NUMDVM + 1, dtype=int)
        self.IDIVB = np.zeros((self.NUMDVM + 1, 3), dtype=int) # Fortran: IDIVB(NUMDVM, 2)
        self.IDIV = np.zeros(self.NUMDVM + 1, dtype=int)

    def compute_thick_from_csa(self, nccm, nx, ils, ile, nlvlsm, nlvls, ilc, x, xcc, dh, acsec, wcsec, csab, csas, lamb, lams, wsb, wss):
        dhl = np.zeros(nlvlsm + 1)
        acsecl = np.zeros(nlvlsm + 1)
        wcsecl = np.zeros(nlvlsm + 1)

        for ii in range(ils, ile + 1):
            xlc = x[ii]
            j = 1
            # Interpolation factors
            # Ensure ilc[ii] is within bounds. Fortran assumes valid.
            lp0 = (xcc[ilc[ii]+1] - xlc) / (xcc[ilc[ii]+1] - xcc[ilc[ii]])
            lp1 = (xlc - xcc[ilc[ii]]) / (xcc[ilc[ii]+1] - xcc[ilc[ii]])

            dhl[j] = dh[ilc[ii], j]*lp0 + dh[ilc[ii]+1, j]*lp1
            acsecl[j] = acsec[ilc[ii], j]*lp0 + acsec[ilc[ii]+1, j]*lp1
            wcsecl[j] = wcsec[ilc[ii], j]*lp0 + wcsec[ilc[ii]+1, j]*lp1

            icount = 0
            for i in range(1, nlvls): # NLVLS-1 loop, i goes 1 to NLVLS-1
                j = i + 1
                csab[ii] = max(csab[ii], 0.0)
                dhl[j] = dh[ilc[ii], j]*lp0 + dh[ilc[ii]+1, j]*lp1
                acsecl[j] = acsec[ilc[ii], j]*lp0 + acsec[ilc[ii]+1, j]*lp1
                wcsecl[j] = wcsec[ilc[ii], j]*lp0 + wcsec[ilc[ii]+1, j]*lp1

                if csab[ii] >= acsecl[i] and csab[ii] < acsecl[i+1]:
                    lamb[ii] = dhl[i] + (dhl[i+1] - dhl[i]) * \
                               (csab[ii] - acsecl[i]) / (acsecl[i+1] - acsecl[i])
                    wsb[ii] = wcsecl[i] + (wcsecl[i+1] - wcsecl[i]) * \
                              (csab[ii] - acsecl[i]) / (acsecl[i+1] - acsecl[i])
                    icount += 1

                csat = csab[ii] + csas[ii]
                csat = max(csat, 0.0)

                if csat >= acsecl[i] and csat < acsecl[i+1]:
                    ltot = dhl[i] + (dhl[i+1] - dhl[i]) * \
                           (csat - acsecl[i]) / (acsecl[i+1] - acsecl[i])
                    wss[ii] = wcsecl[i] + (wcsecl[i+1] - wcsecl[i]) * \
                              (csat - acsecl[i]) / (acsecl[i+1] - acsecl[i])
                    lams[ii] = ltot - lamb[ii]
                    if lams[ii] <= 0.0: lams[ii] = 0.0
                    if wss[ii] < wsb[ii]: wss[ii] = wsb[ii]
                    icount += 1

                if icount == 2:
                    break

    def compute_csa_from_thick(self, nccm, nx, ils, ile, nlvlsm, nlvls, ilc, x, xcc, dh, acsec, wcsec, lamb, lams, csab, csas, wsb, wss):
        dhl = np.zeros(nlvlsm + 1)
        acsecl = np.zeros(nlvlsm + 1)
        wcsecl = np.zeros(nlvlsm + 1)

        for ii in range(ils, ile + 1):
            xlc = x[ii]
            j = 1
            lp0 = (xcc[ilc[ii]+1] - xlc) / (xcc[ilc[ii]+1] - xcc[ilc[ii]])
            lp1 = (xlc - xcc[ilc[ii]]) / (xcc[ilc[ii]+1] - xcc[ilc[ii]])

            dhl[j] = dh[ilc[ii], j]*lp0 + dh[ilc[ii]+1, j]*lp1
            acsecl[j] = acsec[ilc[ii], j]*lp0 + acsec[ilc[ii]+1, j]*lp1
            wcsecl[j] = wcsec[ilc[ii], j]*lp0 + wcsec[ilc[ii]+1, j]*lp1

            icount = 0
            for i in range(1, nlvls):
                j = i + 1
                lamb[ii] = max(lamb[ii], 0.0)
                dhl[j] = dh[ilc[ii], j]*lp0 + dh[ilc[ii]+1, j]*lp1
                acsecl[j] = acsec[ilc[ii], j]*lp0 + acsec[ilc[ii]+1, j]*lp1
                wcsecl[j] = wcsec[ilc[ii], j]*lp0 + wcsec[ilc[ii]+1, j]*lp1

                if lamb[ii] >= dhl[i] and lamb[ii] < dhl[i+1]:
                    csab[ii] = acsecl[i] + (acsecl[i+1] - acsecl[i]) * \
                               (lamb[ii] - dhl[i]) / (dhl[i+1] - dhl[i])
                    wsb[ii] = wcsecl[i] + (wcsecl[i+1] - wcsecl[i]) * \
                              (csab[ii] - acsecl[i]) / (acsecl[i+1] - acsecl[i])
                    icount += 1

                ltot = lamb[ii] + lams[ii]
                ltot = max(ltot, 0.0)

                if ltot >= dhl[i] and ltot < dhl[i+1]:
                    csat = acsecl[i] + (acsecl[i+1] - acsecl[i]) * \
                           (ltot - dhl[i]) / (dhl[i+1] - dhl[i])
                    wss[ii] = wcsecl[i] + (wcsecl[i+1] - wcsecl[i]) * \
                              (csat - acsecl[i]) / (acsecl[i+1] - acsecl[i])
                    csas[ii] = csat - csab[ii]
                    if csas[ii] <= 0.0: csas[ii] = 0.0
                    if wss[ii] < wsb[ii]: wss[ii] = wsb[ii]
                    icount += 1

                if icount == 2:
                    break

    def interface_physics(self, nm, n, grav, swsal, swsed, spgrav, delh, rght, rmaf, evt, ub, ut, sb, st, sdb, sdt, ab, at, wif, wws, phi, cdbed, cdsrf, edif, ddif, entr, dvb, dvt, rich, icall):
        kinvisc = 8.6E-7
        usre = ut[1] * (at[1] / wws[1]) / kinvisc
        if usre < 500.0: usre = 500.0

        for i in range(1, n + 1):
            sttmp = st[i]
            if sttmp > swsal: sttmp = swsal
            if sttmp < 0.0: sttmp = 0.0
            sbtmp = sb[i]
            if sbtmp > swsal: sbtmp = swsal
            if sbtmp < 0.0: sbtmp = 0.0
            sdttmp = sdt[i] / 1000000.0
            if sdttmp > 0.3: sdttmp = 0.3
            if sdttmp < 0.0: sdttmp = 0.0
            sdbtmp = sdb[i] / 1000000.0
            if sdbtmp > 0.3: sdbtmp = 0.3
            if sdbtmp < 0.0: sdbtmp = 0.0

            dent = 1000. + 0.78 * sttmp
            denb = 1000. + 0.78 * sbtmp
            spgt = spgrav * 1000. / dent
            spgb = spgrav * 1000. / denb
            dent = dent * (1. + ((spgt-1.)*sdttmp)/(spgt - (spgt-1.)*sdttmp))
            denb = denb * (1. + ((spgb-1.)*sdbtmp)/(spgb - (spgb-1.)*sdbtmp))

            phi[i] = dent / denb
            if phi[i] > 0.9977: phi[i] = 0.9977

            lbcsa = ab[i] / (wif[i] + 1.E-8)
            ltcsa = at[i] / (wws[i] + 1.E-8)
            depfr = max(lbcsa/2., 1.E-8)
            rought = max(rght, 1.E-8)

            if rought > depfr:
                drh = 1.
            else:
                drh = depfr / rought

            sdfact = 1.0
            if rought > depfr:
                sdfact = max((2.*(depfr/rought) - 1.), 0.0)

            bta = 29.7 * drh
            cdbed[i] = sdfact * 2.*((0.4*bta)/((bta+1.)*(math.log(bta+1.)-1.)+1.))**2.
            ustb = math.sqrt(cdbed[i]/2.) * abs(ub[i])

            depfr = max(ltcsa, 1.E-8)
            rought = max(rght, 1.E-8)
            if rought > depfr:
                drh = 1.
            else:
                drh = depfr / rought

            sdfact = 1.0
            if rought > depfr:
                sdfact = max((2.*(depfr/rought) - 1.), 0.0)

            bta = 29.7 * drh
            cdsrf[i] = sdfact * 2.*((0.4*bta)/((bta+1.)*(math.log(bta+1.)-1.)+1.))**2.
            ustt = math.sqrt(cdsrf[i]/2.) * abs(ut[i])

            depfr = max((ltcsa + lbcsa), 1.E-8)
            bltfr = max(lbcsa, 1.E-8)
            zhev = bltfr / depfr
            vfus = abs(ut[i])
            if abs(ub[i]) > vfus: vfus = abs(ub[i])

            udiff = abs(ut[i] - ub[i])
            vdiff = abs(ut[i] - ub[i])
            if vdiff < .01: vdiff = .01
            rich[i] = 1.0
            entr[i] = 0.0

            sdfact = 1.0
            depfr = max(ltcsa, 1.E-8)
            if rought > depfr:
                sdfact = max((2.*(depfr/rought) - 1.), 0.0)
            depfr = max(lbcsa/2., 1.E-8)
            if rought > depfr:
                sdfact = max((2.*(depfr/rought) - 1.), 0.0)

            if ltcsa > delh and lbcsa > delh:
                cdint = 0.001
                retl = max((ltcsa*vdiff/kinvisc), 1.0)
                for j in range(10):
                    cdint = 2.*(2.5*math.log(4.5*retl*math.sqrt(cdint/2.)))**(-2.0)

                phrat = 2.*(dent-denb)/(dent+denb)
                if phrat > -0.002337: phrat = -0.002337
                rich[i] = -grav*(phrat/2.)*(ltcsa+lbcsa)/(vdiff**2.)
                if rich[i] < 1.0: rich[i] = 1.0
                edif[i] = (1./4.)*cdint*depfr*udiff
                ddif[i] = (1./4.)*cdint*depfr*udiff
                edif[i] = edif[i]/(1.+0.74*rich[i])
                ddif[i] = ddif[i]/(1.+37.*rich[i]*rich[i])

                duoz = math.sqrt(500./usre) + (0.25/(math.sqrt(rich[i]**2.+0.0625)))*(1.-math.sqrt(500./usre))
                entr[i] = 0.038*(1.-rich[i]/(math.sqrt(rich[i]**2.+0.0625))) + (2./usre)*(1./duoz)
                entr[i] = entr[i] / rmaf
            else:
                edif[i] = kinvisc
                ddif[i] = kinvisc
                entr[i] = 0.0

            edif[i] = max(edif[i] * sdfact, kinvisc)
            ddif[i] = max(ddif[i] * sdfact, kinvisc)
            entr[i] = entr[i] * sdfact

            if icall == 0:
                dvt[i] = max((5.93 * ltcsa * ustt), (evt/10.))
                dvb[i] = max((10.06 * lbcsa * ustb), (evt/10.))

    def sediment_physics(self, nm, n, dt, wsd, csd, cse, erc, spgrav, por, ub, ut, sdb, sdt, ab, at, wif, wws, sdepfb, sdepft, serofb, seroft, taugb, taugt, ssdib, ssdit, iets):
        sdconv = 1000000.
        rought = 0.001

        for i in range(1, n + 1):
            lbcsa = ab[i]/(wif[i] + 1.E-8)
            ltcsa = at[i]/(wws[i] + 1.E-8)
            depfr = max(lbcsa/2., 1.E-8)

            if rought > 29.7 * depfr:
                drh = 1. / 29.7
            else:
                drh = depfr / rought

            bta = 29.7 * drh
            cdbed = 2.*((0.4*bta)/((bta+1.)*(math.log(bta+1.)-1.)+1.))**2.
            ustb = math.sqrt(cdbed/2.) * abs(ub[i])
            taugb[i] = 1000.*ustb*ustb

            depfr = max(ltcsa, 1.E-8)
            if rought > 29.7 * depfr:
                drh = 1. / 29.7
            else:
                drh = depfr / rought

            bta = 29.7 * drh
            cdsrf = 2.*((0.4*bta)/((bta+1.)*(math.log(bta+1.)-1.)+1.))**2.
            ustt = math.sqrt(cdsrf/2.) * abs(ut[i])
            taugt[i] = 1000.*ustt*ustt

            if csd < 1.E-8:
                sdepft[i] = 0.0
                sdepfb[i] = 0.0
            else:
                sdepft[i] = wsd*sdt[i]*(1. - min((taugt[i]/csd), 1.0))
                sdepfb[i] = wsd*sdb[i]*(1. - min((taugb[i]/csd), 1.0))

            if cse < 1.E-8:
                seroft[i] = 0.0
                serofb[i] = 0.0
            else:
                seroft[i] = sdconv*erc*(max((taugt[i]/cse), 1.0) - 1.)
                serofb[i] = sdconv*erc*(max((taugb[i]/cse), 1.0) - 1.)

            if ssdit[i] < 1.E-3: seroft[i] = 0.0
            if ssdib[i] < 1.E-3: serofb[i] = 0.0

        if iets == 1:
            for i in range(1, n + 1):
                ssdit[i] = ssdit[i] + dt*(sdepft[i] - seroft[i]) * (wws[i] - wif[i]) / \
                           (sdconv * spgrav * (1. - por))
                ssdib[i] = ssdib[i] + dt*(sdepfb[i] - serofb[i]) * wif[i] / \
                           (sdconv * spgrav * (1. - por))
                ssdit[i] = max(ssdit[i], 0.0)
                ssdib[i] = max(ssdib[i], 0.0)

    def critical_depth_bc(self, i, nm, grav, phi, ab, at, wws, ub, ut, ubaf):
        phrat = 1. - phi[i]
        wdcsa = (at[i] + ab[i]) / (wws[i] + 1.E-8)

        # Ensure non-negative argument for sqrt
        denom_sq = phrat * grav * wdcsa
        if denom_sq < 0: denom_sq = 0

        if denom_sq > 0:
            dfn = ubaf * abs(ut[1]) / (denom_sq**0.5)
        else:
             dfn = 0.0 # Or handle appropriately

        topfrac = (dfn)**(2./3.)

        if topfrac < 0.1: topfrac = 0.1
        if topfrac > 0.9: topfrac = 0.9

        return topfrac

    def run(self):
        print('***************************************************')
        print(' MRSWAT - Mississippi River Salt Wedge Analysis Tool')
        print(' version 3  - 12-24 ')
        print(' written by Gary L. Brown and Phu V. Luong')
        print(' Coastal and Hydraulics Laboratory')
        print(' Engineer Research and Development Center')
        print('***************************************************')
        print()
        print('enter the input file containing the requested run parameters')

        if len(sys.argv) > 1:
            filein = sys.argv[1]
        else:
            filein = input()

        try:
            with open(filein.strip(), 'r') as f:
                lines = [line.strip() for line in f if line.strip()]

            lines_iter = iter(lines)

            def get_val(iter_lines, type_func=float):
                try:
                    next(iter_lines) # Skip header
                    val_line = next(iter_lines)
                    return type_func(val_line.split()[0])
                except StopIteration:
                    return None

            def get_str(iter_lines):
                try:
                    next(iter_lines)
                    return next(iter_lines)
                except StopIteration:
                    return None

            self.DT = get_val(lines_iter)
            self.DX = get_val(lines_iter)
            self.DX = self.DX / self.m2ft
            self.N = get_val(lines_iter, int)
            self.NLVLS = get_val(lines_iter, int)
            self.time = get_val(lines_iter)
            self.time = self.time * 86400.
            freq = get_val(lines_iter)
            freq = freq * 86400.
            self.ifreq = int(freq / self.DT)
            self.RGHT = get_val(lines_iter)
            self.RGHT = self.RGHT / 100.0
            self.evt = get_val(lines_iter)
            self.evt = self.evt / (self.m2ft * self.m2ft)
            self.evb = self.evt
            self.SWSAL = get_val(lines_iter)
            self.RMAF = get_val(lines_iter)
            self.irst = get_val(lines_iter, int)

            self.trst = 0.0
            filere = ""
            if self.irst == 1:
                filere = get_str(lines_iter)
                self.trst = get_val(lines_iter)
                self.trst = self.trst * 86400.

            filews = get_str(lines_iter)
            fileflx = get_str(lines_iter)
            fildiv = get_str(lines_iter)
            filegeo = get_str(lines_iter)
            filop = get_str(lines_iter)
            filic = get_str(lines_iter)
            filtsol = get_str(lines_iter)
            filsn = get_str(lines_iter)
            filtloc = get_str(lines_iter)
            self.ISDFLG = get_val(lines_iter, int)

            filsdi = ""
            filsdb = ""
            filsds = ""
            filsdn = ""
            filsdr = ""
            filsre = ""

            if self.ISDFLG == 1:
                filsdi = get_str(lines_iter)
                filsdb = get_str(lines_iter)
                filsds = get_str(lines_iter)
                filsdn = get_str(lines_iter)
                filsdr = get_str(lines_iter)
                if self.irst == 1:
                    filsre = get_str(lines_iter)

        except Exception as e:
            print(f"Error reading input file: {e}")
            return

        spgrav = 2.65
        wsd = 0.0
        cse = 0.0
        erc = 0.0
        csd = 0.0
        por = 0.3
        drelv = 0.0
        swsed = 0.0

        if self.ISDFLG == 1:
            try:
                with open(filsdi.strip(), 'r') as f:
                    lines_sdi = [line.strip() for line in f if line.strip()]
                iter_sdi = iter(lines_sdi)
                spgrav = get_val(iter_sdi)
                wsd = get_val(iter_sdi)
                wsd = wsd / 1000.0
                cse = get_val(iter_sdi)
                erc = get_val(iter_sdi)
                erc = erc / 1000.0
                csd = get_val(iter_sdi)
                por = get_val(iter_sdi)
                drelv = get_val(iter_sdi)
                drelv = drelv / self.m2ft
                swsed = get_val(iter_sdi)
            except Exception as e:
                print(f"Error reading sediment parameter file: {e}")
                return

        print('     ')
        print(f'You requested the following: ')
        print('     ')
        print(f'Total Number of Days of Simulation {self.time/86400.:14.4f}')
        print(f'Frequency of output in days {freq/86400.:14.4f}')
        print(f'Roughness Height in Centimeters  {self.RGHT*100.:14.4f}')
        print('     ')
        print(f'Filename for downstream wsel info -- {filews}')
        print(f'Filename for upstream discharge info -- {fileflx}')
        print(f'Filename for diversion info -- {fildiv}')
        print(f'Filename for geometry/hypsometry data -- {filegeo}')
        print(f'Filename for observation locations -- {filop}')
        print('     ')
        if self.irst == 1:
            print('Restart Conditions: ')
            print(f'Filename for restart conditions -- {filere}')
            print(f'Restart time in days  {self.trst/86400.0:14.4f}')
            print('     ')

        print(f'Filename for initial conditions output -- {filic}')
        print(f'Filename for full output -- {filtsol}')
        print(f'Filename for output at selected locations -- {filsn}')
        print(f'Filename for toe location output -- {filtloc}')
        print('     ')

        if self.ISDFLG == 1:
            print('SEDIMENT TRANSPORT IS ACTIVE: ')
            print('     ')
            print(f'Filename for sediment parameters -- {filsdi}')
            print(f'Filename for upstream sediment concentrations  -- {filsdb}')
            print(f'Filename for sediment output at selected locations -- {filsdn}')
            print(f'Filename for full sediment output -- {filsds}')
            print(f'Filename for dredge volume output -- {filsdr}')
            if self.irst == 1:
                print(f'Filename for sediment restart conditions -- {filsre}')
            print('     ')

        print('Beginning simulation now ')
        print('     ')

        nstep = int((self.time - self.trst) / self.DT)

        elvi = 0.0
        telva = 0.0
        elva = 0.0
        telvb = 0.0
        elvb = 0.0

        # Read wsel (103) for initial elvi
        lines_103 = []
        try:
             with open(filews.strip(), 'r') as f103:
                lines_103 = [line.strip() for line in f103 if line.strip()]

             if self.irst == 0:
                 iter_103 = iter(lines_103)
                 next(iter_103) # Header
                 line = next(iter_103)
                 parts = line.split()
                 telva = float(parts[0])
                 elva = float(parts[1])
                 line = next(iter_103)
                 parts = line.split()
                 telvb = float(parts[0])
                 elvb = float(parts[1])

                 telva = telva * 86400.
                 telvb = telvb * 86400.
                 elva = elva * 0.3048
                 elvb = elvb * 0.3048
                 elvi = elva
        except Exception as e:
            print(f"Error reading wsel file: {e}")
            return

        depthT = elvi

        # Read hypsometry (105)
        try:
            with open(filegeo.strip(), 'r') as f105:
                lines_geo = [line.strip() for line in f105 if line.strip()]
            iter_geo = iter(lines_geo)

            next(iter_geo) # CHARDUM
            next(iter_geo) # IDUM, RDUM, RDUM2

            xcmax = float(self.N - 1) * self.DX / (5280.0 * 0.3048) - 19.0

            ncc = 0
            for i in range(1, self.NCCM + 1):
                try:
                    next(iter_geo) # CHARDUM

                    line = next(iter_geo)
                    parts = line.split()
                    self.XCC[i] = float(parts[1])
                    self.ZCC[i] = float(parts[2])

                    self.XCC[i] = (xcmax - self.XCC[i]) * 5280.0 * 0.3048

                    next(iter_geo) # CHARDUM

                    for j in range(1, self.NLVLS + 1):
                        line = next(iter_geo)
                        parts = line.split()
                        self.DH[i, j] = float(parts[0])
                        self.ACSEC[i, j] = float(parts[1])
                        self.WCSEC[i, j] = float(parts[2])

                    ncc += 1
                except StopIteration:
                    break

            for i in range(1, self.N + 1):
                self.XMILE[i] = -19.0 + float(self.N - i) * self.DX / (5280.0 * 0.3048)
                self.X[i] = float(i - 1) * self.DX

                for j in range(1, ncc):
                    if self.X[i] >= self.XCC[j] and self.X[i] < self.XCC[j+1]:
                        lp0 = (self.XCC[j+1] - self.X[i]) / (self.XCC[j+1] - self.XCC[j])
                        lp1 = (self.X[i] - self.XCC[j]) / (self.XCC[j+1] - self.XCC[j])
                        self.Z[i] = self.ZCC[j] * lp0 + self.ZCC[j+1] * lp1
                        self.ILC[i] = j

            if self.ISDFLG == 1:
                for i in range(1, self.N + 1):
                    self.LB[i] = max(drelv - self.Z[i], 0.0)
                    self.LT[i] = self.LB[i] + 1.0

                self.compute_csa_from_thick(self.NCCM, self.NM, 1, self.N, self.NLVLSM, self.NLVLS, self.ILC, self.X, self.XCC, self.DH, self.ACSEC, self.WCSEC, self.LB, self.LT, self.ANDRG, self.AT, self.WIF, self.WWS)

        except Exception as e:
            print(f"Error reading geometry file: {e}")
            return

        # Diversions
        numdiv = 0
        try:
             with open(fildiv.strip(), 'r') as f106:
                 lines_div = [line.strip() for line in f106 if line.strip()]
             iter_div = iter(lines_div)

             next(iter_div)
             numdiv = int(next(iter_div).split()[0])

             if numdiv > 0:
                 for i in range(1, numdiv + 1):
                     next(iter_div)
                     line = next(iter_div)
                     parts = line.split()
                     self.DIVRM[i] = float(parts[0])
                     self.DIVFT[i] = float(parts[1])
                     self.DIVW[i] = float(parts[2])
                     self.DIVIE[i] = float(parts[3])
                     self.DIVSAL[i] = float(parts[4])

                     self.DIVW[i] = self.DIVW[i] / self.m2ft
                     self.DIVIE[i] = self.DIVIE[i] / self.m2ft
                     self.DIVFT[i] = -math.log(1.0 - self.DIVFT[i])

                 for i in range(1, numdiv + 1):
                     rdum = 1.E8
                     for ii in range(1, self.N + 1):
                         if abs(self.XMILE[ii] - self.DIVRM[i]) < rdum:
                             self.IDIV[i] = ii
                             rdum = abs(self.XMILE[ii] - self.DIVRM[i])

                 for i in range(1, numdiv + 1):
                     self.IDIVB[i, 1] = self.IDIV[i] - int(self.DIVW[i] / (2. * self.DX))
                     self.IDIVB[i, 2] = self.IDIV[i] + int(self.DIVW[i] / (2. * self.DX))

                 for ii in range(1, self.N + 1):
                     self.STU[ii] = 0.0
                     self.SBU[ii] = 0.0
                     self.SDTU[ii] = 0.0
                     self.SDBU[ii] = 0.0
                     self.DIVF[ii] = 0.0
                     self.DIVIA[ii] = 0.0

                     for i in range(1, numdiv + 1):
                         if ii > self.IDIVB[i, 1] and ii <= self.IDIVB[i, 2]:
                             self.DIVF[ii] = self.DIVF[ii] + self.DIVFT[i] / float(self.IDIVB[i, 2] - self.IDIVB[i, 1])
                             self.LB[ii] = max((self.DIVIE[i] - self.Z[ii]), 0.0)
                             self.SBU[ii] = self.DIVSAL[i]
                             self.LT[ii] = 0.0
                             self.compute_csa_from_thick(self.NCCM, self.NM, ii, ii, self.NLVLSM, self.NLVLS, self.ILC, self.X, self.XCC, self.DH, self.ACSEC, self.WCSEC, self.LB, self.LT, self.DIVIA, self.AT, self.WIF, self.WWS)
        except Exception as e:
             # It's possible fildiv doesn't exist or is empty if no diversions, but code assumes it exists.
             # Actually Fortran opens `status='old'`, so it must exist.
             # But loop `IF (NUMDIV .GT. 0)` handles 0 diversions.
             pass

        for i in range(1, self.N + 1):
            self.UT[i] = 0.0
            self.UB[i] = 0.0
            self.LT[i] = depthT - self.Z[i]
            self.LB[i] = 0.0
            self.ST[i] = 0.0
            self.SB[i] = 0.0
            self.SDT[i] = 0.0
            self.SDB[i] = 0.0
            if i >= (self.N - 1):
                self.LT[i] = 0.1 * (depthT - self.Z[i])
                self.LB[i] = 0.9 * (depthT - self.Z[i])
                self.ST[i] = 0.0
                self.SB[i] = self.SWSAL
                self.SDT[i] = 0.0
                self.SDB[i] = swsed

        self.compute_csa_from_thick(self.NCCM, self.NM, 1, self.N, self.NLVLSM, self.NLVLS, self.ILC, self.X, self.XCC, self.DH, self.ACSEC, self.WCSEC, self.LB, self.LT, self.AB, self.AT, self.WIF, self.WWS)

        ubaf = 1.0
        if numdiv > 0:
            for ii in range(1, self.N + 1):
                for i in range(1, numdiv + 1):
                    if self.IDIV[i] == ii:
                        ubaf = ubaf * (1. - self.DIVFT[i])

        ubaf = 1.43 * ubaf

        # Observation points
        nobsp = 0
        try:
             with open(filop.strip(), 'r') as f106:
                 lines_op = [line.strip() for line in f106 if line.strip()]
             iter_op = iter(lines_op)
             nobsp = int(next(iter_op).split()[0])
             for i in range(1, nobsp + 1):
                 self.ROBRM[i] = float(next(iter_op).split()[0])

             for i in range(1, nobsp + 1):
                 rdum = 1.E8
                 for ii in range(1, self.N + 1):
                     if abs(self.XMILE[ii] - self.ROBRM[i]) < rdum:
                         self.IOBRM[i] = ii
                         rdum = abs(self.XMILE[ii] - self.ROBRM[i])
        except Exception as e:
             pass

        # Restart
        if self.irst == 1:
            try:
                with open(filere.strip(), 'r') as f106:
                    lines_re = [line.strip() for line in f106 if line.strip()]
                iter_re = iter(lines_re)
                next(iter_re) # header

                for i in range(1, self.N + 1):
                    ii = self.N - i + 1
                    line = next(iter_re)
                    parts = line.split()
                    zws = float(parts[2])
                    zt = float(parts[3])
                    zb = float(parts[4])
                    self.ST[ii] = float(parts[5])
                    self.SB[ii] = float(parts[6])
                    ute = float(parts[7])
                    ube = float(parts[8])

                    self.LB[ii] = zt / self.m2ft - self.Z[ii]
                    if self.LB[ii] < 0.0:
                        zt = zt - (self.LB[ii] * self.m2ft)
                        self.LB[ii] = 0.0
                        self.SB[ii] = 0.0
                    self.LT[ii] = (zws - zt) / self.m2ft
                    self.UB[ii] = ube / self.m2ft
                    self.UT[ii] = ute / self.m2ft

                self.compute_csa_from_thick(self.NCCM, self.NM, 1, self.N, self.NLVLSM, self.NLVLS, self.ILC, self.X, self.XCC, self.DH, self.ACSEC, self.WCSEC, self.LB, self.LT, self.AB, self.AT, self.WIF, self.WWS)

                if self.ISDFLG == 1:
                    with open(filsre.strip(), 'r') as f106:
                        lines_sre = [line.strip() for line in f106 if line.strip()]
                    iter_sre = iter(lines_sre)
                    next(iter_sre)
                    for i in range(1, self.N + 1):
                         ii = self.N - i + 1
                         line = next(iter_sre)
                         parts = line.split()
                         self.SDT[ii] = float(parts[5])
                         self.SDB[ii] = float(parts[6])
                         ute = float(parts[9])
                         ube = float(parts[10])

                         self.SSDIB[ii] = ube / (self.m2ft * self.m2ft)
                         self.SSDIT[ii] = ute / (self.m2ft * self.m2ft)
            except Exception as e:
                print(f"Error reading restart file: {e}")
                return

        # Open output files
        f107 = open(filic.strip(), 'w')
        f110 = open(filtsol.strip(), 'w')
        f112 = open(filsn.strip(), 'w')
        f113 = open(filtloc.strip(), 'w')
        f114 = open('diagnostics.txt', 'w')

        f110.write('Time River-mile Wsel(ft) Layer-Interface-el(ft) Bed-el(ft) Surface-salt(ppt) Bottom-salt(ppt) Surface-vel(ft/sec) Bottom-vel(ft/sec) Surface-Q(cfs) Bottom-Q(cfs) \\n')
        f112.write('Time River-mile Wsel(ft) Layer-Interface-el(ft) Bed-el(ft) Surface-salt(ppt) Bottom-salt(ppt) Surface-vel(ft/sec) Bottom-vel(ft/sec) Surface-Q(cfs) Bottom-Q(cfs) \\n')
        f114.write('Time River-mile Horizontal-Disperson-Surface-Layer(sqft/sec) Horizontal-Disperson-Bottom-Layer(sqft/sec) Richardson-Number Vertical-Eddy-Visc(sqft/sec) Vertical-Diff-Coeff(sqft/sec) Entrainment-Cofficient Bed-Drag-Coefficient-Surface-Layer Bed-Drag-Coefficient-Bottom-Layer Water-Surface-Width(ft) Layer-Interface-Width(ft) \\n')

        f115 = None
        f116 = None
        f117 = None

        if self.ISDFLG == 1:
            f115 = open(filsds.strip(), 'w')
            f116 = open(filsdn.strip(), 'w')
            f117 = open(filsdr.strip(), 'w')
            f115.write('Time River-mile Wsel(ft) Layer-Interface-el(ft) Bed-el(ft) Surface-sed(ppm) Bottom-sed(ppm) Surface-grain-shear(Pa) Bottom-grain-shear(Pa) Surface-Layer-Deposition(cubic-feet/foot) Bottom-Layer-Deposition(cubic-feet/foot) Dredging-Requirement(cubic-feet/foot) \\n')
            f116.write('Time River-mile Wsel(ft) Layer-Interface-el(ft) Bed-el(ft) Surface-sed(ppm) Bottom-sed(ppm) Surface-grain-shear(Pa) Bottom-grain-shear(Pa) Surface-Layer-Deposition(cubic-feet/foot) Bottom-Layer-Deposition(cubic-feet/foot) Dredging-Requirement(cubic-feet/foot) \\n')

        ttime = self.trst

        if self.irst == 0:
            f110.write(f'Time = {ttime/86400.:10.4f}\\n')
            if self.ISDFLG == 1:
                f115.write(f'Time = {ttime/86400.:10.4f}\\n')

        f113.write('Time(Days) Toe-Location(River-Mile) 0.25-ppt-at-surface(River-Mile)\\n')
        if self.ISDFLG == 1:
            f117.write('Time(Days) Dredge-Volume(cubic-yards)\\n')

        for i in range(self.N, 0, -1):
            f107.write(f'{i:8d} {ttime/86400.:12.4f} {self.XMILE[i]:12.4f} {self.X[i]:12.4f} {(self.LT[i]+self.LB[i]+self.Z[i])*self.m2ft:12.4f} {(self.LB[i]+self.Z[i])*self.m2ft:12.4f} {self.Z[i]*self.m2ft:12.4f} {self.ST[i]:12.4f} {self.SB[i]:12.4f} {self.UT[i]*self.m2ft:12.4f} {self.UB[i]*self.m2ft:12.4f} {self.AT[i]*self.m2ft*self.m2ft:12.4f} {self.AB[i]*self.m2ft*self.m2ft:10.1f}\\n')

            if self.irst == 0:
                rqt = self.UT[i]*self.AT[i]*self.m2ft**3
                rqb = self.UB[i]*self.AB[i]*self.m2ft**3
                f110.write(f'{ttime/86400.:10.4f} {self.XMILE[i]:10.4f} {(self.LT[i]+self.LB[i]+self.Z[i])*self.m2ft:10.4f} {(self.LB[i]+self.Z[i])*self.m2ft:10.4f} {self.Z[i]*self.m2ft:10.4f} {self.ST[i]:10.4f} {self.SB[i]:10.4f} {self.UT[i]*self.m2ft:10.4f} {self.UB[i]*self.m2ft:10.4f} {rqt:14.4f} {rqb:14.4f}\\n')

                if self.ISDFLG == 1:
                     self.DREDGE[i] = max((self.SSDIB[i] - self.ANDRG[i]), 0.0)
                     f115.write(f'{ttime/86400.:10.4f} {self.XMILE[i]:10.4f} {(self.LT[i]+self.LB[i]+self.Z[i])*self.m2ft:10.4f} {(self.LB[i]+self.Z[i])*self.m2ft:10.4f} {self.Z[i]*self.m2ft:10.4f} {self.SDT[i]:10.4f} {self.SDB[i]:10.4f} {self.TAUGT[i]:12.3f} {self.TAUGB[i]:12.3f} {self.SSDIT[i]*self.m2ft*self.m2ft:12.3f} {self.SSDIB[i]*self.m2ft*self.m2ft:12.3f} {self.DREDGE[i]*self.m2ft*self.m2ft:12.3f}\\n')

        mu = self.DT / self.DX

        try:
            with open(fileflx.strip(), 'r') as f104:
                lines_104 = [line.strip() for line in f104 if line.strip()]
            iter_104 = iter(lines_104)
            next(iter_104)

            with open(filews.strip(), 'r') as f103:
                lines_103 = [line.strip() for line in f103 if line.strip()]
            iter_103 = iter(lines_103)
            next(iter_103)

            lines_sdb = []
            iter_sdb = None
            if self.ISDFLG == 1:
                with open(filsdb.strip(), 'r') as f106_sed:
                    lines_sdb = [line.strip() for line in f106_sed if line.strip()]
                iter_sdb = iter(lines_sdb)
                next(iter_sdb)
        except Exception as e:
            print(f"Error reading time series files: {e}")
            return

        tflxa, fluxa = 0.0, 0.0
        tflxb, fluxb = 0.0, 0.0
        telva, elva = 0.0, 0.0
        telvb, elvb = 0.0, 0.0
        tsdca, sdcona = 0.0, 0.0
        tsdcb, sdconb = 0.0, 0.0

        def get_pair(iter_obj):
            try:
                line = next(iter_obj)
                parts = line.split()
                return float(parts[0]), float(parts[1])
            except StopIteration:
                return None, None

        if self.irst == 1:
            tflxa, fluxa = get_pair(iter_104)
            tflxb, fluxb = get_pair(iter_104)
            tflxa *= 86400.
            tflxb *= 86400.
            fluxa *= 0.3048**3
            fluxb *= 0.3048**3
            fluxi = fluxa

            while (ttime < tflxa) or (ttime >= tflxb):
                tflxa = tflxb
                fluxa = fluxb
                fluxi = fluxa
                tflxb, fluxb = get_pair(iter_104)
                if tflxb is None: break
                tflxb *= 86400.
                fluxb *= 0.3048**3

            telva, elva = get_pair(iter_103)
            telvb, elvb = get_pair(iter_103)
            telva *= 86400.
            telvb *= 86400.
            elva *= 0.3048
            elvb *= 0.3048
            elvi = elva

            while (ttime < telva) or (ttime >= telvb):
                telva = telvb
                elva = elvb
                elvi = elva
                telvb, elvb = get_pair(iter_103)
                if telvb is None: break
                telvb *= 86400.
                elvb *= 0.3048

            if self.ISDFLG == 1:
                tsdca, sdcona = get_pair(iter_sdb)
                tsdcb, sdconb = get_pair(iter_sdb)
                tsdca *= 86400.
                tsdcb *= 86400.
                sdconi = sdcona

                while (ttime < tsdca) or (ttime >= tsdcb):
                    tsdca = tsdcb
                    sdcona = sdconb
                    sdconi = sdcona
                    tsdcb, sdconb = get_pair(iter_sdb)
                    if tsdcb is None: break
                    tsdcb *= 86400.
        else:
             tflxa, fluxa = get_pair(iter_104)
             tflxb, fluxb = get_pair(iter_104)
             tflxa *= 86400.
             fluxa *= 0.3048**3
             tflxb *= 86400.
             fluxb *= 0.3048**3
             fluxi = fluxa

             telva, elva = get_pair(iter_103)
             telvb, elvb = get_pair(iter_103)
             telva *= 86400.
             elva *= 0.3048
             telvb *= 86400.
             elvb *= 0.3048

             if self.ISDFLG == 1:
                 tsdca, sdcona = get_pair(iter_sdb)
                 tsdcb, sdconb = get_pair(iter_sdb)
                 tsdca *= 86400.
                 tsdcb *= 86400.
                 sdconi = sdcona

        flux = fluxi
        sdcon = 0.0

        for ii in range(1, nstep + 1):
            self.UB[1] = 0.0
            self.UT[self.N+1] = self.UT[self.N]
            self.UB[self.N+1] = self.UB[self.N]
            self.LT[self.N+1] = self.LT[self.N]
            self.LB[self.N+1] = self.LB[self.N]
            self.AT[self.N+1] = self.AT[self.N]
            self.AB[self.N+1] = self.AB[self.N]
            self.WIF[self.N+1] = self.WIF[self.N]
            self.WWS[self.N+1] = self.WWS[self.N]
            self.ST[self.N+1] = self.ST[self.N]
            self.SB[self.N+1] = self.SB[self.N]
            self.SDT[self.N+1] = self.SDT[self.N]
            self.SDB[self.N+1] = self.SDB[self.N]

            self.UT[0] = self.UT[1]
            self.UB[0] = self.UB[1]
            self.LT[0] = self.LT[1]
            self.LB[0] = self.LB[1]
            self.AT[0] = self.AT[1]
            self.AB[0] = self.AB[1]
            self.WIF[0] = self.WIF[1]
            self.WWS[0] = self.WWS[1]
            self.ST[0] = self.ST[1]
            self.SB[0] = self.SB[1]
            self.SDT[0] = self.SDT[1]
            self.SDB[0] = self.SDB[1]

            flux = fluxi
            depthT = elvi
            if self.ISDFLG == 1: sdcon = sdconi

            self.interface_physics(self.NM, self.N, self.grav, self.SWSAL, swsed, spgrav, self.DELH, self.RGHT, self.RMAF, self.evt, self.UB, self.UT, self.SB, self.ST, self.SDB, self.SDT, self.AB, self.AT, self.WIF, self.WWS, self.PHI, self.CDBED, self.CDSRF, self.EDIF, self.DDIF, self.ENTR, self.DVB, self.DVT, self.RICH, 0)

            if self.ISDFLG == 1:
                self.sediment_physics(self.NM, self.N, self.DT, wsd, csd, cse, erc, spgrav, por, self.UB, self.UT, self.SDB, self.SDT, self.AB, self.AT, self.WIF, self.WWS, self.SDEPFB, self.SDEPFT, self.SEROFB, self.SEROFT, self.TAUGB, self.TAUGT, self.SSDIB, self.SSDIT, 0)

            # Predictor step
            for i in range(1, self.N):
                self.ATS[i] = self.AT[i] - mu*(self.UT[i+1]*self.AT[i+1]) + mu*(self.UT[i]*self.AT[i])
                self.ABBS[i] = self.AB[i] - mu*(self.UB[i+1]*self.AB[i+1]) + mu*(self.UB[i]*self.AB[i])

                abaiv = max((self.AB[i] - self.DIVIA[i]), 0.0)
                self.ATS[i] = self.ATS[i] - self.DIVF[i]*self.UT[i]*self.AT[i]*self.AT[i]*mu / (self.AT[i] + abaiv + 1.E-8)
                self.ABBS[i] = self.ABBS[i] - self.DIVF[i]*self.UB[i]*self.AB[i]*abaiv*mu / (self.AT[i] + abaiv + 1.E-8)

                entfac = 1.0
                self.ATS[i] = self.ATS[i] + self.ENTR[i]*abs(self.UT[i])*self.DT*self.WIF[i]*entfac
                self.ABBS[i] = self.ABBS[i] - self.ENTR[i]*abs(self.UT[i])*self.DT*self.WIF[i]*entfac

            for i in range(2, self.N):
                vdiff = self.UT[i] - self.UB[i]
                wdcsa = (self.AT[i] + self.AB[i]) / (self.WWS[i] + 1.E-8) + 1.E-8

                self.UTL[i] = (self.UT[i]*self.AT[i]) \
                    - mu*((self.UT[i+1]**2.)*self.AT[i+1] - (self.UT[i]**2.)*self.AT[i]) \
                    - mu*self.grav*(self.AT[i])*(self.LT[i+1]-self.LT[i]) \
                    - mu*self.grav*(self.AT[i])*(self.LB[i+1]-self.LB[i]) \
                    - mu*self.grav*(self.AT[i])*(self.Z[i+1]-self.Z[i]) \
                    + (mu*self.evt/self.DX)*(self.AT[i]*self.UT[i-1] \
                    - 2.0*self.AT[i]*self.UT[i] + self.AT[i]*self.UT[i+1] \
                    + (self.UT[i] - self.UT[i+1]) * (self.AT[i] - self.AT[i+1]))

                if self.LB[i] > self.DELH and self.LT[i] > self.DELH:
                    self.UTL[i] -= self.DT*self.EDIF[i]*2.*vdiff*self.WIF[i]/wdcsa
                    self.UTL[i] -= self.DT*(self.CDSRF[i]/2.)*(self.WWS[i] - self.WIF[i])*self.UT[i]*abs(self.UT[i])
                elif self.LB[i] <= self.DELH and self.LT[i] > self.DELH:
                    self.UTL[i] -= self.DT*(self.CDSRF[i]/2.)*self.WWS[i]*self.UT[i]*abs(self.UT[i])

                self.STL[i] = self.ST[i]*self.AT[i] - mu*(self.UT[i+1]*self.AT[i+1]*self.ST[i+1]) \
                    + mu*(self.UT[i]*self.AT[i]*self.ST[i]) \
                    + (mu*self.DVT[i]/self.DX)*(self.AT[i]*self.ST[i-1] \
                    - 2.0*self.AT[i]*self.ST[i] + self.AT[i]*self.ST[i+1] \
                    + (self.ST[i] - self.ST[i+1]) * (self.AT[i] - self.AT[i+1])) \
                    + (mu*self.AT[i]/self.DX)*((self.ST[i] - self.ST[i+1]) * (self.DVT[i] - self.DVT[i+1]))

                if self.LB[i] > self.DELH and self.LT[i] > self.DELH:
                    self.STL[i] -= self.DT*self.DDIF[i]*2.*(self.ST[i]-self.SB[i])*self.WIF[i]/wdcsa

                entfac = 1.0
                self.STL[i] += self.SB[i]*self.ENTR[i]*abs(self.UT[i])*self.DT*self.WIF[i]*entfac

                if self.ISDFLG == 1:
                    self.SDTL[i] = self.SDT[i]*self.AT[i] - mu*(self.UT[i+1]*self.AT[i+1]*self.SDT[i+1]) \
                        + mu*(self.UT[i]*self.AT[i]*self.SDT[i]) \
                        + (mu*self.DVT[i]/self.DX)*(self.AT[i]*self.SDT[i-1] \
                        - 2.0*self.AT[i]*self.SDT[i] + self.AT[i]*self.SDT[i+1] \
                        + (self.SDT[i] - self.SDT[i+1]) * (self.AT[i] - self.AT[i+1])) \
                        + (mu*self.AT[i]/self.DX)*((self.SDT[i] - self.SDT[i+1]) * (self.DVT[i] - self.DVT[i+1]))

                    if self.LB[i] > 100.*wsd*self.DT and self.LT[i] > 100.*wsd*self.DT:
                        self.SDTL[i] -= self.DT*self.DDIF[i]*2.*(self.SDT[i]-self.SDB[i])*self.WIF[i]/wdcsa
                        self.SDTL[i] -= self.DT*wsd*self.SDT[i]*self.WIF[i]

                    if self.LT[i] > 100.*wsd*self.DT:
                        self.SDTL[i] += self.DT*(self.SEROFT[i] - self.SDEPFT[i])*(self.WWS[i] - self.WIF[i])

                    entfac = 1.0
                    self.SDTL[i] += self.SDB[i]*self.ENTR[i]*abs(self.UT[i])*self.DT*self.WIF[i]*entfac

                abaiv = max((self.AB[i] - self.DIVIA[i]), 0.0)
                if self.UT[i] > 0.0:
                    self.STL[i] -= self.DIVF[i]*self.UT[i]*self.AT[i]*self.AT[i]*self.ST[i]*mu / (self.AT[i] + abaiv + 1.E-8)
                else:
                    self.STL[i] -= self.DIVF[i]*self.UT[i]*self.AT[i]*self.AT[i]*self.STU[i]*mu / (self.AT[i] + abaiv + 1.E-8)

                if self.ISDFLG == 1:
                    if self.UT[i] > 0.0:
                        self.SDTL[i] -= self.DIVF[i]*self.UT[i]*self.AT[i]*self.AT[i]*self.SDT[i]*mu / (self.AT[i] + abaiv + 1.E-8)
                    else:
                        self.SDTL[i] -= self.DIVF[i]*self.UT[i]*self.AT[i]*self.AT[i]*self.SDTU[i]*mu / (self.AT[i] + abaiv + 1.E-8)

                self.UBL[i] = (self.UB[i]*self.AB[i]) \
                    - mu*((self.UB[i+1]**2.)*self.AB[i+1] - (self.UB[i]**2.)*self.AB[i]) \
                    - mu*self.grav*(self.AB[i])*(self.LB[i+1]-self.LB[i]) \
                    - mu*self.grav*self.PHI[i]*(self.AB[i])*(self.LT[i+1]-self.LT[i]) \
                    - mu*self.grav*(self.AB[i])*(self.Z[i+1]-self.Z[i]) \
                    + (mu*self.evb/self.DX)*(self.AB[i]*self.UB[i-1] \
                    - 2.0*self.AB[i]*self.UB[i] + self.AB[i]*self.UB[i+1] \
                    + (self.UB[i] - self.UB[i+1]) * (self.AB[i] - self.AB[i+1])) \
                    - self.DT*(self.CDBED[i]/2.)*self.WIF[i]*self.UB[i]*abs(self.UB[i])

                if self.LT[i] > self.DELH and self.LB[i] > self.DELH:
                    self.UBL[i] += self.DT*self.EDIF[i]*2.*vdiff*self.WIF[i]/wdcsa

                self.SBL[i] = self.SB[i]*self.AB[i] - mu*(self.UB[i+1]*self.AB[i+1]*self.SB[i+1]) \
                    + mu*(self.UB[i]*self.AB[i]*self.SB[i]) \
                    + (mu*self.DVB[i]/self.DX)*(self.AB[i]*self.SB[i-1] \
                    - 2.0*self.AB[i]*self.SB[i] + self.AB[i]*self.SB[i+1] \
                    + (self.SB[i] - self.SB[i+1]) * (self.AB[i] - self.AB[i+1])) \
                    + (mu*self.AB[i]/self.DX)*((self.SB[i] - self.SB[i+1]) * (self.DVB[i] - self.DVB[i+1]))

                if self.LB[i] > self.DELH and self.LT[i] > self.DELH:
                    self.SBL[i] -= self.DT*self.DDIF[i]*2.*self.WIF[i]*(self.SB[i]-self.ST[i])/wdcsa

                entfac = 1.0
                self.SBL[i] -= self.SB[i]*self.ENTR[i]*abs(self.UT[i])*self.DT*self.WIF[i]*entfac

                if self.ISDFLG == 1:
                    self.SDBL[i] = self.SDB[i]*self.AB[i] - mu*(self.UB[i+1]*self.AB[i+1]*self.SDB[i+1]) \
                        + mu*(self.UB[i]*self.AB[i]*self.SDB[i]) \
                        + (mu*self.DVB[i]/self.DX)*(self.AB[i]*self.SDB[i-1] \
                        - 2.0*self.AB[i]*self.SDB[i] + self.AB[i]*self.SDB[i+1] \
                        + (self.SDB[i] - self.SDB[i+1]) * (self.AB[i] - self.AB[i+1])) \
                        + (mu*self.AB[i]/self.DX)*((self.SDB[i] - self.SDB[i+1]) * (self.DVB[i] - self.DVB[i+1]))

                    if self.LB[i] > 100.*wsd*self.DT and self.LT[i] > 100.*wsd*self.DT:
                        self.SDBL[i] -= self.DT*self.DDIF[i]*2.*self.WIF[i]*(self.SDB[i]-self.SDT[i])/wdcsa
                        self.SDBL[i] += self.DT*wsd*self.SDT[i]*self.WIF[i]

                    if self.LB[i] > 100.*wsd*self.DT:
                        self.SDBL[i] += self.DT*(self.SEROFB[i] - self.SDEPFB[i])*self.WIF[i]

                    entfac = 1.0
                    self.SDBL[i] -= self.SDBS[i]*self.ENTR[i]*abs(self.UTS[i])*self.DT*self.WIF[i]*entfac

                abaiv = max((self.ABBS[i] - self.DIVIA[i]), 0.0)
                if self.UBS[i] > 0.0:
                    self.SBL[i] -= self.DIVF[i]*self.UBS[i]*self.ABBS[i]*abaiv*self.SBS[i]*mu / (self.ATS[i] + abaiv + 1.E-8)
                else:
                    self.SBL[i] -= self.DIVF[i]*self.UBS[i]*self.ABBS[i]*abaiv*self.SBU[i]*mu / (self.ATS[i] + abaiv + 1.E-8)

                if self.ISDFLG == 1:
                    if self.UBS[i] > 0.0:
                        self.SDBL[i] -= self.DIVF[i]*self.UBS[i]*self.ABBS[i]*abaiv*self.SDBS[i]*mu / (self.ATS[i] + abaiv + 1.E-8)
                    else:
                        self.SDBL[i] -= self.DIVF[i]*self.UBS[i]*self.ABBS[i]*abaiv*self.SDBU[i]*mu / (self.ATS[i] + abaiv + 1.E-8)

            for i in range(1, self.N + 1):
                delcsa = self.DELH*(1. + self.WWS[i])
                if self.ATSS[i] > delcsa:
                    self.UTSS[i] = self.UTL[i]/self.ATSS[i]
                    self.STSS[i] = self.STL[i]/self.ATSS[i]
                    self.SDTSS[i] = self.SDTL[i]/self.ATSS[i]
                else:
                    self.UTSS[i] = 0.0
                    self.SBSS[i] = (self.STL[i] + self.SBL[i]) / (self.ATSS[i] + self.ABSS[i])
                    self.SDBSS[i] = (self.SDTL[i] + self.SDBL[i]) / (self.ATSS[i] + self.ABSS[i])
                    self.ABSS[i] = self.ABSS[i] + self.ATSS[i]
                    self.ATSS[i] = 0.0
                    self.STSS[i] = 0.0
                    self.SDTSS[i] = 0.0

                delcsa = self.DELH*(1. + self.WIF[i])
                if self.ABSS[i] > delcsa:
                    self.UBSS[i] = self.UBL[i]/self.ABSS[i]
                    self.SBSS[i] = self.SBL[i]/self.ABSS[i]
                    self.SDBSS[i] = self.SDBL[i]/self.ABSS[i]
                else:
                    self.UBSS[i] = 0.0
                    self.STSS[i] = (self.STL[i] + self.SBL[i]) / (self.ATSS[i] + self.ABSS[i])
                    self.SDTSS[i] = (self.SDTL[i] + self.SDBL[i]) / (self.ATSS[i] + self.ABSS[i])
                    self.ATSS[i] = self.ABSS[i] + self.ATSS[i]
                    self.ABSS[i] = 0.0
                    self.SBSS[i] = 0.0
                    self.SDBSS[i] = 0.0

            self.compute_thick_from_csa(self.NCCM, self.NM, 1, self.N, self.NLVLSM, self.NLVLS, self.ILC, self.X, self.XCC, self.DH, self.ACSEC, self.WCSEC, self.ABSS, self.ATSS, self.LBSS, self.LTSS, self.WIF, self.WWS)

            delcsa = self.DELH*(1. + self.WWS[1])
            if self.ATSS[1] > delcsa:
                self.UTSS[1] = flux / self.ATSS[1]
            else:
                self.UTSS[1] = 0.0
                self.ATSS[1] = 0.0
                self.LTSS[1] = 0.0

            self.UBSS[1] = 0.0
            self.STSS[1] = 0.0
            self.SDTSS[1] = sdcon
            self.SBSS[1] = self.SBSS[2]
            self.SDBSS[1] = self.SDBSS[2]

            if self.ATSS[self.N] > 1.E-8:
                self.UTSS[self.N] = self.UTSS[self.N-1] * self.ATSS[self.N-1] / self.ATSS[self.N]
            else:
                self.UTSS[self.N] = 0.0

            if self.ABSS[self.N] > 1.E-8:
                self.UBSS[self.N] = self.UBSS[self.N-1] * self.ABSS[self.N-1] / self.ABSS[self.N]
            else:
                self.UBSS[self.N] = 0.0

            if self.UTSS[self.N] >= 0 and abs(self.UTSS[self.N]*self.ATSS[self.N]) > 1.E-8:
                self.STSS[self.N] = self.STSS[self.N-1] * self.UTSS[self.N-1] * self.ATSS[self.N-1] / (self.UTSS[self.N] * self.ATSS[self.N])
            else:
                self.STSS[self.N] = self.ST[self.N]

            if self.UBSS[self.N] >= 0 and abs(self.UBSS[self.N]*self.LBSS[self.N]) > 1.E-8:
                self.SBSS[self.N] = self.SBSS[self.N-1] * self.UBSS[self.N-1] * self.ABSS[self.N-1] / (self.UBSS[self.N] * self.ABSS[self.N])
            else:
                self.SBSS[self.N] = self.SWSAL

            if self.UTSS[self.N] >= 0 and abs(self.UTSS[self.N]*self.ATSS[self.N]) > 1.E-8:
                self.SDTSS[self.N] = self.SDTSS[self.N-1] * self.UTSS[self.N-1] * self.ATSS[self.N-1] / (self.UTSS[self.N] * self.ATSS[self.N])
            else:
                self.SDTSS[self.N] = self.SDT[self.N]

            if self.UBSS[self.N] >= 0 and abs(self.UBSS[self.N]*self.LBSS[self.N]) > 1.E-8:
                self.SDBSS[self.N] = self.SDBSS[self.N-1] * self.UBSS[self.N-1] * self.ABSS[self.N-1] / (self.UBSS[self.N] * self.ABSS[self.N])
            else:
                self.SDBSS[self.N] = swsed

            for i in range(1, self.N + 1):
                self.LTOLD[i] = self.LT[i]
                self.LBOLD[i] = self.LB[i]
                self.AT[i] = 0.5*(self.AT[i] + self.ATSS[i])
                self.AB[i] = 0.5*(self.AB[i] + self.ABSS[i])

            self.compute_thick_from_csa(self.NCCM, self.NM, 1, self.N, self.NLVLSM, self.NLVLS, self.ILC, self.X, self.XCC, self.DH, self.ACSEC, self.WCSEC, self.AB, self.AT, self.LB, self.LT, self.WIF, self.WWS)

            for i in range(1, self.N + 1):
                if self.LTOLD[i] <= 1.E-8: self.ST[i] = self.STSS[i]
                if self.LBOLD[i] <= 1.E-8: self.SB[i] = self.SBSS[i]
                if self.LTOLD[i] <= 1.E-8: self.SDT[i] = self.SDTSS[i]
                if self.LBOLD[i] <= 1.E-8: self.SDB[i] = self.SDBSS[i]

                if self.LT[i] > 1.E-8 and abs(0.5*(self.UT[i]+self.UTSS[i])) < (1./(2.*mu)) and \
                   (self.PHI[i] < 0.9970 or (self.LT[i]/(self.LB[i] + 1.E-8)) > 0.001):
                    self.UT[i] = 0.5*(self.UT[i] + self.UTSS[i])
                    self.ST[i] = 0.5*(self.ST[i] + self.STSS[i])
                    self.SDT[i] = 0.5*(self.SDT[i] + self.SDTSS[i])
                else:
                    self.AT[i] = max(self.AT[i], 0.0)
                    self.ST[i] = 0.5*(self.ST[i] + self.STSS[i])
                    self.SB[i] = 0.5*(self.SB[i] + self.SBSS[i])
                    self.SDT[i] = 0.5*(self.SDT[i] + self.SDTSS[i])
                    self.SDB[i] = 0.5*(self.SDB[i] + self.SDBSS[i])
                    self.SB[i] = (self.ST[i]*self.AT[i] + self.SB[i]*self.AB[i]) / (self.AT[i] + self.AB[i])
                    self.SDB[i] = (self.SDT[i]*self.AT[i] + self.SDB[i]*self.AB[i]) / (self.AT[i] + self.AB[i])
                    self.LB[i] = self.LB[i] + max(self.LT[i], 0.0)
                    self.AB[i] = self.AB[i] + self.AT[i]
                    self.LT[i] = 0.0
                    self.AT[i] = 0.0
                    self.UT[i] = 0.0
                    self.ST[i] = 0.0
                    self.SDT[i] = 0.0

                if self.LB[i] > 1.E-8 and abs(0.5*(self.UB[i]+self.UBSS[i])) < (1./(2.*mu)) and \
                   (self.PHI[i] < 0.9970 or (self.LB[i]/(self.LT[i] + 1.E-8)) > 0.001):
                    self.UB[i] = 0.5*(self.UB[i] + self.UBSS[i])
                    self.SB[i] = 0.5*(self.SB[i] + self.SBSS[i])
                    self.SDB[i] = 0.5*(self.SDB[i] + self.SDBSS[i])
                else:
                    self.AB[i] = max(self.AB[i], 0.0)
                    self.ST[i] = 0.5*(self.ST[i] + self.STSS[i])
                    self.SB[i] = 0.5*(self.SB[i] + self.SBSS[i])
                    self.SDT[i] = 0.5*(self.SDT[i] + self.SDTSS[i])
                    self.SDB[i] = 0.5*(self.SDB[i] + self.SDBSS[i])
                    self.ST[i] = (self.ST[i]*self.AT[i] + self.SB[i]*self.AB[i]) / (self.AT[i] + self.AB[i])
                    self.SDT[i] = (self.SDT[i]*self.AT[i] + self.SDB[i]*self.AB[i]) / (self.AT[i] + self.AB[i])
                    self.LT[i] = self.LT[i] + max(self.LB[i], 0.0)
                    self.AT[i] = self.AT[i] + self.AB[i]
                    self.LB[i] = 0.0
                    self.AB[i] = 0.0
                    self.UB[i] = 0.0
                    self.SB[i] = 0.0
                    self.SDB[i] = 0.0

                if self.ST[i] < 0.0: self.ST[i] = 0.0
                if self.ST[i] > (self.SWSAL + 2.0): self.ST[i] = self.SWSAL + 2.0
                if self.SB[i] < 0.0: self.SB[i] = 0.0
                if self.SB[i] > (self.SWSAL + 2.0): self.SB[i] = self.SWSAL + 2.0
                if self.SDT[i] < 0.0: self.SDT[i] = 0.0
                if self.SDB[i] < 0.0: self.SDB[i] = 0.0

                if self.UT[i] > 0.0: self.SDTU[i] = self.SDT[i]
                if self.UB[i] > 0.0: self.SDBU[i] = self.SDB[i]

            if self.ISDFLG == 1:
                self.sediment_physics(self.NM, self.N, self.DT, wsd, csd, cse, erc, spgrav, por, self.UB, self.UT, self.SDB, self.SDT, self.AB, self.AT, self.WIF, self.WWS, self.SDEPFB, self.SDEPFT, self.SEROFB, self.SEROFT, self.TAUGB, self.TAUGT, self.SSDIB, self.SSDIT, 1)

            ttime = ii * self.DT + self.trst
            if ii % self.ifreq == 0:
                print(f"Writing output for time(days) = {ttime/86400.:.4f}")
                f110.write(f'Time = {ttime/86400.:10.4f}\\n')
                f114.write(f'Time = {ttime/86400.:10.4f}\\n')
                if self.ISDFLG == 1:
                    f115.write(f'Time = {ttime/86400.:10.4f}\\n')
                    f116.write(f'Time = {ttime/86400.:10.4f}\\n')

                toeloc = self.XMILE[self.N]
                dwvloc = self.XMILE[self.N]
                drvol = 0.0
                rgap = 0.0

                for i in range(self.N, 0, -1):
                    rqt = self.UT[i]*self.AT[i]*self.m2ft**3
                    rqb = self.UB[i]*self.AB[i]*self.m2ft**3
                    f110.write(f'{ttime/86400.:10.4f} {self.XMILE[i]:10.4f} {(self.LT[i]+self.LB[i]+self.Z[i])*self.m2ft:10.4f} {(self.LB[i]+self.Z[i])*self.m2ft:10.4f} {self.Z[i]*self.m2ft:10.4f} {self.ST[i]:10.4f} {self.SB[i]:10.4f} {self.UT[i]*self.m2ft:14.4f} {self.UB[i]*self.m2ft:14.4f} {rqt:14.4f} {rqb:14.4f}\\n')

                    if self.LB[i] > 1.0 and self.SB[i] > 9.0 and rgap < 2.0:
                        toeloc = self.XMILE[i]
                        rgap = 0.0
                    else:
                        if i < self.N-1:
                            rgap = rgap + (self.XMILE[i] - self.XMILE[i+1])

                    if self.ST[i] > 0.25:
                        dwvloc = self.XMILE[i]

                    f114.write(f'{ttime/86400.:10.4f} {self.XMILE[i]:10.4f} {self.DVT[i]*self.m2ft**2:10.6f} {self.DVB[i]*self.m2ft**2:10.6f} {self.RICH[i]:10.6f} {self.EDIF[i]*self.m2ft**2:10.6f} {self.DDIF[i]*self.m2ft**2:10.6f} {self.ENTR[i]:10.6f} {self.CDSRF[i]:10.6f} {self.CDBED[i]:10.6f} {self.WWS[i]*self.m2ft:10.2f} {self.WIF[i]*self.m2ft:10.2f}\\n')

                    if self.ISDFLG == 1:
                        self.DREDGE[i] = max((self.SSDIB[i] - self.ANDRG[i]), 0.0)
                        drvol = drvol + self.DREDGE[i]*self.DX
                        f115.write(f'{ttime/86400.:10.4f} {self.XMILE[i]:10.4f} {(self.LT[i]+self.LB[i]+self.Z[i])*self.m2ft:10.4f} {(self.LB[i]+self.Z[i])*self.m2ft:10.4f} {self.Z[i]*self.m2ft:10.4f} {self.SDT[i]:12.3f} {self.SDB[i]:12.3f} {self.TAUGT[i]:12.3f} {self.TAUGB[i]:12.3f} {self.SSDIT[i]*self.m2ft*self.m2ft:12.3f} {self.SSDIB[i]*self.m2ft*self.m2ft:12.3f} {self.DREDGE[i]*self.m2ft*self.m2ft:12.3f}\\n')

                if nobsp > 0:
                    for j in range(1, nobsp + 1):
                        i = self.IOBRM[j]
                        rqt = self.UT[i]*self.AT[i]*self.m2ft**3
                        rqb = self.UB[i]*self.AB[i]*self.m2ft**3
                        f112.write(f'{ttime/86400.:10.4f} {self.XMILE[i]:10.4f} {(self.LT[i]+self.LB[i]+self.Z[i])*self.m2ft:10.4f} {(self.LB[i]+self.Z[i])*self.m2ft:10.4f} {self.Z[i]*self.m2ft:10.4f} {self.ST[i]:10.4f} {self.SB[i]:10.4f} {self.UT[i]*self.m2ft:14.4f} {self.UB[i]*self.m2ft:14.4f} {rqt:14.4f} {rqb:14.4f}\\n')

                        if self.ISDFLG == 1:
                            self.DREDGE[i] = max((self.SSDIB[i] - self.ANDRG[i]), 0.0)
                            f116.write(f'{ttime/86400.:10.4f} {self.XMILE[i]:10.4f} {(self.LT[i]+self.LB[i]+self.Z[i])*self.m2ft:10.4f} {(self.LB[i]+self.Z[i])*self.m2ft:10.4f} {self.Z[i]*self.m2ft:10.4f} {self.SDT[i]:12.3f} {self.SDB[i]:12.3f} {self.TAUGT[i]:12.3f} {self.TAUGB[i]:12.3f} {self.SSDIT[i]*self.m2ft*self.m2ft:12.3f} {self.SSDIB[i]*self.m2ft*self.m2ft:12.3f} {self.DREDGE[i]*self.m2ft*self.m2ft:12.3f}\\n')

                f113.write(f'{ttime/86400.:10.3f} {toeloc:10.3f} {dwvloc:10.3f}\\n')
                if self.ISDFLG == 1:
                    f117.write(f'{ttime/86400.:10.3f} {drvol*self.m2ft**3/27.0:12.1f}\\n')

            if (ttime >= tflxa) and (ttime <= tflxb):
                delflux = (fluxb-fluxa)*(ttime-tflxa)/(tflxb-tflxa)
                fluxi = fluxa + delflux

            if tflxb is not None and (ttime - tflxb) > 1.E-8:
                tflxa = tflxb
                fluxa = fluxb
                fluxi = fluxa
                tflxb, fluxb = get_pair(iter_104)
                if tflxb is not None:
                    tflxb *= 86400.
                    fluxb *= 0.3048**3
                else:
                    # Keep using last values if EOF
                    pass

            if self.ISDFLG == 1:
                if (ttime >= tsdca) and (ttime <= tsdcb):
                    delsdc = (sdconb-sdcona)*(ttime-tsdca)/(tsdcb-tsdca)
                    sdconi = sdcona + delsdc

                if tsdcb is not None and (ttime - tsdcb) > 1.E-8:
                    tsdca = tsdcb
                    sdcona = sdconb
                    sdconi = sdcona
                    tsdcb, sdconb = get_pair(iter_sdb)
                    if tsdcb is not None:
                        tsdcb *= 86400.

            if (ttime >= telva) and (ttime <= telvb):
                delv = (elvb-elva)*(ttime-telva)/(telvb-telva)
                elvi = elva + delv

            if telvb is not None and (ttime - telvb) > 1.E-8:
                telva = telvb
                elva = elvb
                elvi = elva
                telvb, elvb = get_pair(iter_103)
                if telvb is not None:
                    telvb *= 86400.
                    elvb *= 0.3048

        f113.close()
        f110.close()
        f112.close()
        f114.close()
        f107.close()

        if self.ISDFLG == 1:
            f115.close()
            f116.close()
            f117.close()

        print('     ')
        print('Simulation complete ')

if __name__ == "__main__":
    model = MRSWAT()
    model.run()

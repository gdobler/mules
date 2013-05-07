import numpy as np
from ml_counters import *


# -------- global utility (note: The purpose of glones was to
#          implement a sum via dot(vector,glones[ind]), which is the
#          optimal method for summing vectors.  Surprisingly, there
#          turned out to be no significant speed compared to a simple
#          for loop for the sums.  As such, glones is depricated
#          here.)
mmrange = counters.mmrange
glones  = counters.glones


# -------- allocate arrays and utilities outside of the recursion
sqrt     = np.sqrt
dot      = np.dot
arctan2  = np.arctan2
cos      = np.cos
sin      = np.sin
dt       = np.float64
nmm      = 20
mm       = mmrange[0:nmm] + 1.0
amp      = np.zeros(11, dtype=dt)
cmom     = np.zeros(20, dtype=dt)
smom     = np.zeros(20, dtype=dt)
oormp1   = np.zeros(20, dtype=dt)
cmmtcell = np.zeros(20, dtype=dt)
smmtcell = np.zeros(20, dtype=dt)
alphar   = np.zeros(20, dtype=dt)
alphat   = np.zeros(20, dtype=dt)

def ml_defang_iter(ximg, yimg, cells, stars, cellnum):

    # -------- utilities
    maxnst  = counters.maxnst
    mi      = np.zeros(maxnst, dtype=dt)
    delxi   = np.zeros(maxnst, dtype=dt) + 1e-6 # avoid divide by 0
    delyi   = np.zeros(maxnst, dtype=dt) + 1e-6 # avoid divide by 0
    mi2ori2 = np.zeros(maxnst, dtype=dt)
    alphax  = 0.0
    alphay  = 0.0
    clist   = [0,1,2,3]


    
    # -------- loop through the cell list (inserting subcells when 
    #          necessary).
    for cellnum in clist:

        # -------- unpack the cell
        tcind   = clist.index(cellnum)
        delx    = ximg - cells[cellnum].xcell
        dely    = yimg - cells[cellnum].ycell
        rcell   = sqrt(delx*delx + dely*dely)
        high    = cells[cellnum].high

        # -------- you're obviously inside the cell, break it up
        if (rcell<1.0e-5) and ~high:
            clist.insert(tcind+1,4*(cellnum+1) + 3)
            clist.insert(tcind+1,4*(cellnum+1) + 2)
            clist.insert(tcind+1,4*(cellnum+1) + 1)
            clist.insert(tcind+1,4*(cellnum+1) + 0)

            continue

        # -------- continue unpacking
        tcell   = arctan2(dely,delx)
        mcell   = cells[cellnum].mcell
        cmom[:] = cells[cellnum].cmom
        smom[:] = cells[cellnum].smom
        nstar   = cells[cellnum].nstar


        # -------- there are no stars in this cell
        if nstar==0: continue


        # -------- calculate the deflection angle up to order 20
        oormp1[:]   = rcell**(-mm-1.0)
        alphar0     = mcell/rcell
        cmmtcell[:] = cos(mm*tcell)
        smmtcell[:] = sin(mm*tcell)
        alphar[:]   = (cmom*cmmtcell+smom*smmtcell)*oormp1
        alphat[:]   = (cmom*smmtcell-smom*cmmtcell)*oormp1


        # -------- large multipole moments are too high, break up cell
        amp[:] = sqrt(alphar*alphar+alphat*alphat)[nmm-11:nmm]
        bcond  = 0.0
        for i in xrange(11): bcond += amp[i]

        if bcond>=0.1:

            # -------- already at the highest depth cell, use all stars
            if high:
                stind = cells[cellnum].stind
                stind = stind[stind > 0]
                stsz  = stind.size

                counters.starcnt += stsz

                mi[:stsz]    = stars.rein[stind]
                delxi[:stsz] = ximg - stars.xstar[stind]
                delyi[:stsz] = yimg - stars.ystar[stind]
                mi2ori2[:]   = (mi*mi)/(delxi*delxi + delyi*delyi)

                for i in xrange(stsz):
                    alphax += (delxi*mi2ori2)[i]
                    alphay += (delyi*mi2ori2)[i]

                counters.clist.append(-cellnum)
            else:
                clist.insert(tcind+1,4*(cellnum+1) + 3)
                clist.insert(tcind+1,4*(cellnum+1) + 2)
                clist.insert(tcind+1,4*(cellnum+1) + 1)
                clist.insert(tcind+1,4*(cellnum+1) + 0)

                continue
        else:
            counters.cellcnt += 1

            talphar   = alphar0
            talphat   = 0.0

            for i in xrange(20):
                talphar += alphar[i]
                talphat += alphat[i]

            alphax += (delx*talphar - dely*talphat)/rcell
            alphay += (dely*talphar + delx*talphat)/rcell

            counters.clist.append(cellnum)

    return alphax, alphay

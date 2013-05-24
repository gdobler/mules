import numpy as np
from ml_counters import *

def ml_defang_iter(ximg, yimg, cells, stars, cellnum):

    # -------- allocate arrays and utilities outside of the iteration
    sqrt     = np.sqrt
    dot      = np.dot
    arctan2  = np.arctan2
    cos      = np.cos
    sin      = np.sin
    dt       = np.float64
    nmm      = 20
    mm       = np.arange(nmm) + 1.0
    amp      = np.zeros(11, dtype=dt)
    cmom     = np.zeros(20, dtype=dt)
    smom     = np.zeros(20, dtype=dt)
    oormp1   = np.zeros(20, dtype=dt)
    cmmtcell = np.zeros(20, dtype=dt)
    smmtcell = np.zeros(20, dtype=dt)
    alphar   = np.zeros(20, dtype=dt)
    alphat   = np.zeros(20, dtype=dt)


    # -------- utilities
    maxnst  = counters.maxnst
    mi      = np.zeros(maxnst, dtype=dt)
    delxi   = np.zeros(maxnst, dtype=dt) + 1e-6 # avoid divide by 0
    delyi   = np.zeros(maxnst, dtype=dt) + 1e-6 # avoid divide by 0
    mi2ori2 = np.zeros(maxnst, dtype=dt)
    alphax  = 0.0
    alphay  = 0.0
    clist   = [0,1,2,3]

    # -------- utility function
    def useallstars(cell, alphax, alphay):
        stind = cell.stind
        stind = stind[stind > 0]
        stsz  = stind.size

        counters.starcnt += stsz

        mi[:stsz]    = stars.rein[stind]
        delxi[:stsz] = stars.xstar[stind] - ximg
        delyi[:stsz] = stars.ystar[stind] - yimg
        mi2ori2[:]   = (mi*mi)/(delxi*delxi + delyi*delyi)
        axvec        = delxi*mi2ori2
        ayvec        = delyi*mi2ori2

        for i in xrange(stsz):
            alphax += axvec[i]
            alphay += ayvec[i]

        counters.clist.append(-cellnum)

        return alphax, alphay

    
    # -------- loop through the cell list (inserting subcells when 
    #          necessary).
    for cellnum in clist:
        # -------- grab the cell
        cell = cells[cellnum]


        # -------- are there stars in this cell?
        nstar = cell.nstar
        if nstar==0: continue


        # -------- unpack the cell
        tcind   = clist.index(cellnum)
        delx    = cell.xcell - ximg
        dely    = cell.ycell - yimg
        rcell   = sqrt(delx*delx + dely*dely)
        high    = cell.high
        cside   = 0.5*(cell.xside+cell.yside)

        # -------- break up the cell
        if (rcell<0.6*1.5*cside):
            if high!=1:
                clist.insert(tcind+1,4*(cellnum+1) + 3)
                clist.insert(tcind+1,4*(cellnum+1) + 2)
                clist.insert(tcind+1,4*(cellnum+1) + 1)
                clist.insert(tcind+1,4*(cellnum+1) + 0)
                continue
            else:
                alphax, alphay = useallstars(cell,alphax,alphay)
                continue

        # -------- continue unpacking
        tcell   = arctan2(dely,delx)
        mcell   = cell.mcell
        cmom[:] = cell.cmom
        smom[:] = cell.smom


        # -------- calculate the deflection angle up to order 20
        oormp1[:]   = rcell**(-mm-1.0)
        alphar0     = mcell/rcell
        cmmtcell[:] = cos(mm*tcell)
        smmtcell[:] = sin(mm*tcell)
        alphar[:]   = (cmom*cmmtcell+smom*smmtcell)*oormp1
        alphat[:]   = (cmom*smmtcell-smom*cmmtcell)*oormp1


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

import numpy as np
from ml_counters import *

# -------- allocate arrays and utilities outside of the recursion
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



# -------- utility function
def useallstars(ximg, yimg, cell, stars, cellnum, alphax, alphay):
    stind = cell.stind
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

    return alphax, alphay



# -------- break up cell
def breakup_cell(ximg, yimg, cells, stars, cellnum):
    basecell = 0 if cellnum==None else 4*(cellnum+1)

    defang0 = ml_defang(ximg,yimg,cells,stars,basecell+0)
    defang1 = ml_defang(ximg,yimg,cells,stars,basecell+1)
    defang2 = ml_defang(ximg,yimg,cells,stars,basecell+2)
    defang3 = ml_defang(ximg,yimg,cells,stars,basecell+3)
    defang  = defang0[0] + defang1[0] + defang2[0] + defang3[0], \
              defang0[1] + defang1[1] + defang2[1] + defang3[1]

    return defang



# -------- main recursion function
def ml_defang(ximg, yimg, cells, stars, cellnum):

    # -------- begin running lowest level cell
    if cellnum==None:
        global mi, delxi, delyi, mi2ori2

        maxnst  = counters.maxnst
        mi      = np.zeros(maxnst, dtype=dt)
        delxi   = np.zeros(maxnst, dtype=dt) + 1e-6 # avoid divide by 0
        delyi   = np.zeros(maxnst, dtype=dt) + 1e-6 # avoid divide by 0
        mi2ori2 = np.zeros(maxnst, dtype=dt)

        return breakup_cell(ximg, yimg, cells, stars, cellnum)


    # -------- grab the cell and intialize the deflection angle for it
    cell   = cells[cellnum]
    alphax = 0.0
    alphay = 0.0
    defang = alphax, alphay


    # -------- are there stars in this cell?
    nstar = cell.nstar
    if nstar==0: return defang


    # -------- utilities
    delx   = ximg - cell.xcell
    dely   = yimg - cell.ycell
    rcell  = sqrt(delx*delx + dely*dely)
    high   = cell.high

    if (rcell<1.0e-5) and high!=1:
        return breakup_cell(ximg, yimg, cells, stars, cellnum)
    elif (rcell<1.0e-5) and high==1:
        return useallstars(ximg,yimg,cell,stars,cellnum,alphax,alphay)

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


    # -------- large multipole moments are too high, break up cell
    amp[:] = sqrt(alphar*alphar+alphat*alphat)[nmm-11:nmm]
    bcond  = 0.0
    for i in xrange(11): bcond += amp[i]

    if bcond>=0.1:
        # -------- already at the highest depth cell, use all stars
        if high:
            defang = useallstars(ximg,yimg,cell,stars,cellnum,alphax,alphay)
        else:
            return breakup_cell(ximg, yimg, cells, stars, cellnum)
    else:
        counters.cellcnt += 1

        talphar   = alphar0
        talphat   = 0.0

        for i in xrange(20):
            talphar += alphar[i]
            talphat += alphat[i]

        alphax = (delx*talphar - dely*talphat)/rcell
        alphay = (dely*talphar + delx*talphat)/rcell
        defang = alphax, alphay

        counters.clist.append(cellnum)


    # -------- return the deflection angle
    return defang

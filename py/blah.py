import numpy as np
from ml_counters import *

# -------- global utility
mmrange = counters.mmrange


def ml_defang(ximg, yimg, cells, stars, cellnum):

# -------- initialize the deflection angle and run lowest level cell
    if cellnum==None:
        defang0 = ml_defang(ximg, yimg, cells, stars, 0)
        defang1 = ml_defang(ximg, yimg, cells, stars, 1)
        defang2 = ml_defang(ximg, yimg, cells, stars, 2)
        defang3 = ml_defang(ximg, yimg, cells, stars, 3)

#        print "ML_DEFANG: MUST SORT CLIST!!!"

        return defang0[0] + defang1[0] + defang2[0] + defang3[0], \
            defang0[1] + defang1[1] + defang2[1] + defang3[1]



# -------- utilities
    delx   = ximg - cells[cellnum].xcell
    dely   = yimg - cells[cellnum].ycell
    rcell  = np.sqrt(delx*delx + dely*dely)
    tcell  = np.arctan2(dely,delx)
    mcell  = cells[cellnum].mcell
    cmom   = cells[cellnum].cmom
    smom   = cells[cellnum].smom
    nstar  = cells[cellnum].nstar
    alphax = 0.0
    alphay = 0.0
    defang = alphax, alphay

    if nstar==0: return defang



# -------- calculate the deflection angle up to order 20
    nmm      = cmom.size
    mm       = mmrange[0:nmm] + 1.0
    oormp1   = 1.0/rcell**(mm+1.0)
    alphar0  = mcell/rcell
    cmmtcell = np.cos(mm*tcell)
    smmtcell = np.sin(mm*tcell)
    alphar   = (cmom*cmmtcell+smom*smmtcell)*oormp1
    alphat   = (cmom*smmtcell-smom*cmmtcell)*oormp1

    foo = (np.sqrt(alphar*alphar+alphat*alphat))[nmm-11:nmm]
    bar = foo[0] + foo[1] +foo[2] + foo[3] + foo[4] + foo[5] + foo[6] + \
        foo[7] + foo[8] + foo[9]

    # -------- large multipole moments are too high, break up cell
#    if np.add.reduce((np.sqrt(alphar*alphar+alphat*alphat))[nmm-11:nmm])>=0.1:
    if bar>=0.1:

        # -------- alread at the highest depth cell, use all stars
        if cells[cellnum].high:
            stind = cells[cellnum].stind
            stind = stind[stind > 0]

            counters.starcnt += stind.size

            mi        = stars.rein[stind]
            mi2       = mi*mi
            delxi     = ximg - stars.xstar[stind]
            delyi     = yimg - stars.ystar[stind]
            ri2       = delxi*delxi + delyi*delyi
#            alphax    = np.add.reduce(mi2*(delxi/ri2))
#            alphay    = np.add.reduce(mi2*(delyi/ri2))
            defang    = alphax, alphay

#            counters.clist.append(-cellnum)

        else:
            defang0 = ml_defang(ximg,yimg,cells,stars,4*(cellnum+1))
            defang1 = ml_defang(ximg,yimg,cells,stars,4*(cellnum+1)+1)
            defang2 = ml_defang(ximg,yimg,cells,stars,4*(cellnum+1)+2)
            defang3 = ml_defang(ximg,yimg,cells,stars,4*(cellnum+1)+3)
            defang  = defang0[0] + defang1[0] + defang2[0] + defang3[0], \
                defang0[1] + defang1[1] + defang2[1] + defang3[1]
    else:
        counters.cellcnt += 1

#        talphar   = alphar0 + np.add.reduce(alphar)
#        talphat   = np.add.reduce(alphat)
        talphar   = 0.0
        talphat   = 0.0
#        alphax    = np.cos(tcell)*talphar - np.sin(tcell)*talphat
#        alphay    = np.sin(tcell)*talphar + np.cos(tcell)*talphat
        alphax    = delx/rcell*talphar - dely/rcell*talphat
        alphay    = dely/rcell*talphar + delx/rcell*talphat
        defang    = alphax, alphay

#        counters.clist.append(cellnum)

    return defang

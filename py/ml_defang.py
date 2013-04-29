import numpy as np
from ml_counters import *

# -------- global utility
mmrange = counters.mmrange
glones  = counters.glones

def ml_defang(ximg, yimg, cells, stars, cellnum):

    # -------- break up cell
    def breakup_cell(cellnum):
        basecell = 0 if cellnum==None else 4*(cellnum+1)
        defang0  = ml_defang(ximg,yimg,cells,stars,basecell+0)
        defang1  = ml_defang(ximg,yimg,cells,stars,basecell+1)
        defang2  = ml_defang(ximg,yimg,cells,stars,basecell+2)
        defang3  = ml_defang(ximg,yimg,cells,stars,basecell+3)
        defang   = defang0[0] + defang1[0] + defang2[0] + defang3[0], \
                   defang0[1] + defang1[1] + defang2[1] + defang3[1]

        return defang



    # -------- begin running lowest level cell
    if cellnum==None: return breakup_cell(cellnum)



    # -------- utilities
    delx   = ximg - cells[cellnum].xcell
    dely   = yimg - cells[cellnum].ycell
    rcell  = np.sqrt(delx*delx + dely*dely)

    if (rcell<1.0e-5) and ~cells[cellnum].high: return breakup_cell(cellnum)

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
    oormp1   = rcell**(-mm-1.0)
    alphar0  = mcell/rcell
    cmmtcell = np.cos(mm*tcell)
    smmtcell = np.sin(mm*tcell)
    alphar   = (cmom*cmmtcell+smom*smmtcell)*oormp1
    alphat   = (cmom*smmtcell-smom*cmmtcell)*oormp1

    # -------- large multipole moments are too high, break up cell
    bcond = np.dot((np.sqrt(alphar*alphar + \
                                alphat*alphat))[nmm-11:nmm],glones[0:11])

    if bcond>=0.1:

        # -------- alread at the highest depth cell, use all stars
        if cells[cellnum].high:
            stind = cells[cellnum].stind
            stind = stind[stind > 0]

            counters.starcnt += stind.size

            mi        = stars.rein[stind]
            mi2       = mi*mi
            delxi     = ximg - stars.xstar[stind]
            delyi     = yimg - stars.ystar[stind]
            oori2     = 1.0/(delxi*delxi + delyi*delyi)
            alphax    = np.dot(mi2*delxi*oori2,glones[:stind.size])
            alphay    = np.dot(mi2*delyi*oori2,glones[:stind.size])
            defang    = alphax, alphay

            counters.clist.append(-cellnum)

        else:
            return breakup_cell(cellnum)
    else:
        counters.cellcnt += 1

        talphar   = np.dot(alphar,glones[0:20]) + alphar0
        talphat   = np.dot(alphat,glones[0:20])
        alphax    = (delx*talphar - dely*talphat)/rcell
        alphay    = (dely*talphar + delx*talphat)/rcell
        defang    = alphax, alphay

        counters.clist.append(cellnum)

    return defang

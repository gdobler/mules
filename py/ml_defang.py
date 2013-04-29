from numpy import sqrt, dot, arctan2, cos, sin
from ml_counters import *

# -------- global utility (note: The purpose of glones was to
#          implement a sum via dot(vector,glones[ind]), which is the
#          optimal method for summing vectors.  Surprisingly, there
#          turned out to be no significant speed compared to a simple
#          for loop for the sums.  As such, glones is depricated
#          here.)
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
    rcell  = sqrt(delx*delx + dely*dely)
    high   = cells[cellnum].high

    if (rcell<1.0e-5) and ~high: return breakup_cell(cellnum)

    tcell  = arctan2(dely,delx)
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
    cmmtcell = cos(mm*tcell)
    smmtcell = sin(mm*tcell)
    alphar   = (cmom*cmmtcell+smom*smmtcell)*oormp1
    alphat   = (cmom*smmtcell-smom*cmmtcell)*oormp1

    # -------- large multipole moments are too high, break up cell
    amp   = sqrt(alphar*alphar+alphat*alphat)[nmm-11:nmm]
    bcond = 0.0
    for i in xrange(11): bcond += amp[i]

    if bcond>=0.1:

        # -------- already at the highest depth cell, use all stars
        if high:
            stind = cells[cellnum].stind
            stind = stind[stind > 0]
            stsz  = stind.size

            counters.starcnt += stsz

            mi        = stars.rein[stind]
            mi2       = mi*mi
            delxi     = ximg - stars.xstar[stind]
            delyi     = yimg - stars.ystar[stind]
            oori2     = 1.0/(delxi*delxi + delyi*delyi)

            for i in xrange(stsz):
                alphax += mi2[i]*delxi[i]*oori2[i]
                alphay += mi2[i]*delyi[i]*oori2[i]

            defang    = alphax, alphay

            counters.clist.append(-cellnum)

        else:
            return breakup_cell(cellnum)
    else:
        counters.cellcnt += 1

        talphar   = alphar0
        talphat   = 0.0

        for i in xrange(20):
            talphar += alphar[i]
            talphat += alphat[i]

        alphax    = (delx*talphar - dely*talphat)/rcell
        alphay    = (dely*talphar + delx*talphat)/rcell
        defang    = alphax, alphay

        counters.clist.append(cellnum)

    return defang

#+
# NAME:
#
# PURPOSE:
#
# CALLING SEQUENCE:
#
# INPUTS:
#
# OPTIONAL INPUTS:
#
# KEYWORDS:
#
# OUTPUTS:
#
# OPTIONAL OUTPUTS:
#
# EXAMPLES:
#
# COMMENTS:
#
# REVISION HISTORY:
#  2012/01/19 - Written by Greg Dobler (KITP/UCSB)
#
#------------------------------------------------------------
def ml_gencells(kappas, xrange, yrange, nside, xseed=xseed, yseed=yseed, 
                rseed=rseed, stars=stars):

# -------- utilities
    from math import *
    import ml_cellpos

    xmin  = xrange[0]
    xmax  = xrange[1]
    dx    = xmax - xmin
    ymin  = yrange[0]
    ymax  = yrange[1]
    dy    = ymax - ymin
    nstar = kappas*dx*dy/pi # Nstar = kappa*Lx*Ly/pi*bbar^2
    cpos  = ml_cellpos(xrange, yrange, nside)
    nlev  = long(log(nside)/log(2))

    print 'Creating cell structure with ' + str(nstar) + ' stars...'




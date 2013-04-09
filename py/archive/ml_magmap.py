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
#  2012/01/18 - Written by Greg Dobler (KITP/UCSB)
#
#------------------------------------------------------------
def ml_magmap(kappa, gamma, fstar, nxpix=300L, nypix=300L, xr=[-15.,15.], 
              yr=[-15.,15.], xseed=3L, yseed=31L, rseed=314L, nside=16L, 
              nximg=1001L, nyimg=1001L, eps=0.05, scheme='scheme', bin='bin', 
              defarr='defarr', xsrc='xsrc', ysrc='ysrc', ixsrc='ixsrc', 
              iysrc='iysrc', stars='stars', cells='cells'):

# -------- set arrays
    from numpy import array
    xr = array(xr)
    yr = array(yr)



# -------- utilities
    kappas   = kappa*fstar
    kappac   = kappa - kappas
    xmin     = xr[0]
    xmax     = xr[1]
    dx       = (xmax-xmin)/float(nxpix)
    ymin     = yr[0]
    ymax     = yr[1]
    dy       = (ymax-ymin)/float(nypix)
    xrimg    = 1.2*xr / abs(1.0 - kappa - gamma) # emperical factor
    yrimg    = 1.2*yr / abs(1.0 - kappa + gamma) # emperical factor
    xminimg  = xrimg[0]
    xmaximg  = xrimg[1]
    yminimg  = yrimg[0]
    ymaximg  = yrimg[1]
    dximg    = (xmaximg - xminimg)/float(nximg)
    dyimg    = (ymaximg - yminimg)/float(nyimg)
    xrst     = [min(xrimg[0],yrimg[0]),max(xrimg[1],yrimg[1])]
    yrst     = xrst
    nraymag1 = dx*dy/(dximg*dyimg)

    print ('Shooting region size is ' + str(dximg*float(nximg)) + ' by ' + 
           str(dyimg*float(nyimg)))

    print ('Number of rays per pixel for unity magnification is ' + 
            str(nraymag1))

    print 'size of stars region: ' + str(xrst) + str(yrst)



# -------- initialize the stars and cells
#    import ml_gencells(kappas, xrst, yrst, nside, xseed=xseed, yseed=yseed,
#                       rseed=rseed, stars=stars)

    return

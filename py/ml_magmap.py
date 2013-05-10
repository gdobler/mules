import numpy as np
from ml_gencells import *
from ml_genstars import *
from ml_counters import *
from ml_defarr import *

class magmap():

    """
      magmap class includes the following data.
        ??? : ???
        ??? : ???
        ??? : ???
    """

    # -------- initialize the magnification map parameters
    def __init__(self, kappa, gamma, fstar, nxpix=None, nypix=None,
                 xr=None, yr=None, seed_pos=None, seed_rein=None, nside=None,
                 nximg=None, nyimg=None, eps=None, scheme=None, bins=None,
                 xsrc=None, ysrc=None, ixsrc=None, iysrc=None):


        """ Initialize the magnification map parameters """
        # --------  defaults
        nxpix     = 300 if nxpix==None else nxpix
        nypix     = 300 if nypix==None else nypix
        xr        = [-15., 15.] if xr==None else xr # in stellar Rein
        yr        = [-15., 15.] if yr==None else yr # in stellar Rein
        nside     = 16 if nside==None else nside # nside of stellar population
        nximg     = 100 if nximg==None else nximg
        nyimg     = 100 if nyimg==None else nyimg
        eps       = 0.05 if eps==None else eps # % of deflection of 6th moment
        seed_pos  = 111 if seed_pos==None else seed_pos
        seed_rein = 222 if seed_rein==None else seed_rein


        # -------- utilities
        xr, yr = np.array([xr,yr])

        kappas = kappa*fstar
        kappac = kappa - kappas

        [xmin, xmax], [ymin, ymax] = xr, yr
        dx, dy = (xmax-xmin)/float(nxpix), (ymax-ymin)/float(nypix)

        xrimg = 1.2*xr / abs(1.0 - kappa - gamma) # emperical factor
        yrimg = 1.2*yr / abs(1.0 - kappa + gamma) # emperical factor
        [xminimg, xmaximg], [yminimg, ymaximg] = xrimg, yrimg
        dximg = (xmaximg-xminimg)/float(nximg)
        dyimg = (ymaximg-yminimg)/float(nyimg)

        xrst = np.array([min(xrimg[0],yrimg[0]),max(xrimg[1],yrimg[1])])*2.2
        yrst = xrst

        nraymag1 = dx*dy/(dximg*dyimg)


        # -------- messages
        print "ML_MAGMAP: Shooting region size: ", dximg*nximg, ' by ', \
            dyimg*nyimg

        print "ML_MAGMAP: Number of rays per pixel for unity " + \
            "magnification: ", nraymag1

        print "ML_MAGMAP: Size of stars region: ", xrst, yrst


        # -------- initialize the stars and cells
        stars  = genstars(kappas, xrst, yrst, nside, seed_pos, seed_rein)
        cells  = gencells(stars)
        nstars = np.array([i.nstar for i in cells])

        counters.maxnst = max(nstars[np.array([i.high for i in cells])==1])

        print "ML_MAGMAP: Mean number of stars per highest cell: ", \
            np.mean(nstars[np.array([i.high for i in cells])==1])


        # -------- initialize the deflection angle map
        defarr = np.zeros([nximg, nyimg, 2])
        xpos   = xmin + dx*np.arange(nxpix)
        ypos   = ymin + dy*np.arange(nypix)


        # -------- shoot the rays in the image plane (calculate the
        #          deflection angle at each point
        ml_defarr(defarr, xrimg, yrimg, nximg, nyimg, cells, stars, \
                      scheme=scheme, bins=bins)

        print "ML_MAGMAP: Total number of cells used: ", counters.cellcnt
        print "ML_MAGMAP: Total number of stars used: ", counters.starcnt

        # -------- gather rays into source plane pixels
        #          note: definition of meshgrid backwards
        ximg, yimg = np.meshgrid(np.linspace(xrimg[0],xrimg[1],nximg), \
                                     np.linspace(yrimg[0],yrimg[1],nyimg))
        ximg, yimg = ximg.T, yimg.T

        xsrc  = ximg*(1.0-kappac-gamma) - defarr[:,:,0]
        ysrc  = yimg*(1.0-kappac+gamma) - defarr[:,:,1]
        ixsrc = (np.floor((xsrc-xmin)/dx).astype(int)).reshape(nximg*nyimg)
        iysrc = (np.floor((ysrc-ymin)/dy).astype(int)).reshape(nximg*nyimg)

        w = np.where((ixsrc >= 0) & (ixsrc < nxpix) & (iysrc >= 0) & \
                         (iysrc < nypix))[0]

        if w.size==0:
            print "ML_MAGMAP: NO RAYS DEFLECTED INTO SOURCE PLANE GRID!!!"
            return

        isrc_ring = np.bincount(ixsrc[w] + iysrc[w]*nxpix)        


        # -------- convert to magnification
        magarr = np.zeros(nxpix*nypix)
        magarr[:len(isrc_ring)] += isrc_ring
        magarr = magarr.reshape(nxpix,nypix)

        print "ML_MAGMAP: Number of rays shot: ", nximg*nyimg
        print "ML_MAGMAP: Number of rays landed: ", w.size
        print "ML_MAGMAP: Double check: ", np.sum(magarr)

        magarr /= nraymag1


        # -------- set attributes
        self.nxpix     = nxpix
        self.nypix     = nypix
        self.xr        = xr
        self.yr        = yr
        self.nside     = nside
        self.nximg     = nximg
        self.nyimg     = nyimg
        self.eps       = eps
        self.seed_pos  = seed_pos
        self.seed_rein = seed_rein
        self.nraymag1  = nraymag1
        self.stars     = stars
        self.cells     = cells
        self.defarr    = defarr
        self.magarr    = magarr



    # -------- convolve the magnification map with a source morphology
#    def convmap(self, type=None):

# type 1 is a Gaussian
#
# set the source size
#
# determine the units of the source size in pixels (compared to Rein).
# this is going to involve ratios of number of pixels.
#
# set a minimum number of pixels for the convolution.  for example
# problems could arise in which the source size is like 0.1, but
# magarr is only sampled at 1 Rein...  what to do in that case?
#
# do the convolution
#
# what to do if we don't want a Gaussian?
#
# uniform source is possible
#
# random little clouds
#
# passing a convolution kernel should be allowed
#
# we should store the postage same of the convolution kernel
#
# we should also store some parameters for the convolution kernel.
# what should those be?  rsrc is an obvious one, but what else?
#
# need some support functions to go from Rein to physical
# size/wavelength

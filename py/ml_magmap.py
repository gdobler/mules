import sys
import numpy as np
from scipy import ndimage, interpolate
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
                 xr=None, yr=None, seed_pos=None, seed_rein=None,
                 nside=None, nximg=None, nyimg=None, eps=None,
                 multi=None, recur=None, bins=None, xsrc=None,
                 ysrc=None, ixsrc=None, iysrc=None, beta=None,
                 mbar=None, mrat=None):


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
        bins      = 1 if bins==None else bins # subgridding in deflection angle
        seed_pos  = 111 if seed_pos==None else seed_pos
        seed_rein = 222 if seed_rein==None else seed_rein
        beta      = 0.0 if beta==None else beta
        mbar      = 1.0 if mbar==None else mbar
        mrat      = 1.0 if mrat==None else mrat


        # -------- make sure source plane is a square
        if (nxpix!=nypix) or (xr[0]!=yr[0]) or (xr[1]!=yr[1]):
            print("ML_MAGMAP: ERROR - source plane must be square...")
            print("ML_MAGMAP:         i.e., nxpix=nypix, xr=yr")
            sys.exit(-1)


        # -------- utilities
        xr, yr = np.array([xr,yr])

        kappas = kappa*fstar
        kappac = kappa - kappas

        [xmin, xmax], [ymin, ymax] = xr, yr
        dx, dy = (xmax-xmin)/float(nxpix), (ymax-ymin)/float(nypix)

#        xrimg = 1.2*xr / abs(1.0 - kappa - gamma) # emperical factor
#        yrimg = 1.2*yr / abs(1.0 - kappa + gamma) # emperical factor
        xrimg = xr / abs(1.0 - kappa - gamma) # emperical factor
        yrimg = yr / abs(1.0 - kappa + gamma) # emperical factor

        xrimg[0] -= 2.0
        xrimg[1] += 2.0
        yrimg[0] -= 2.0
        yrimg[1] += 2.0

        [xminimg, xmaximg], [yminimg, ymaximg] = xrimg, yrimg
        dximg = (xmaximg-xminimg)/float(nximg*bins)
        dyimg = (ymaximg-yminimg)/float(nyimg*bins)

        xrst = np.array([min(xrimg[0],yrimg[0]),max(xrimg[1],yrimg[1])])*2.2
        yrst = xrst

        nraymag1 = dx*dy/(dximg*dyimg)


        # -------- messages
        print "ML_MAGMAP: Shooting region size: ", dximg*nximg*bins, ' by ', \
            dyimg*nyimg*bins

        print "ML_MAGMAP: Number of rays per pixel for unity " + \
            "magnification: ", nraymag1

        print "ML_MAGMAP: Size of stars region: ", xrst, yrst


        # -------- initialize the stars and cells
        stars = genstars(kappas, xrst, yrst, nside, seed_pos,
                         seed_rein, beta=beta, mbar=mbar, mrat=mrat)
        cells  = gencells(stars)
        nstars = np.array([i.nstar for i in cells])

        counters.maxnst = max(nstars[np.array([i.high for i in cells])==1])

        print "ML_MAGMAP: Mean number of stars per highest cell: ", \
            np.mean(nstars[np.array([i.high for i in cells])==1])


        # -------- initialize the deflection angle map
        xpos   = xmin + dx*np.arange(nxpix)
        ypos   = ymin + dy*np.arange(nypix)


        # -------- shoot the rays in the image plane (calculate the
        #          deflection angle at each point
        defarr = ml_defarr(xrimg, yrimg, nximg, nyimg, cells, stars, \
                           multi=multi, recur=recur, bins=bins)

        if (len(defarr)==1) and (defarr==-1):
            print("ML_MAGMAP: ERROR: # of processors (multi kw) ")
            print("ML_MAGMAP:   inappropriately set... aborting.")
            return

        print "ML_MAGMAP: Total number of cells used: ", counters.cellcnt
        print "ML_MAGMAP: Total number of stars used: ", counters.starcnt

        # -------- gather rays into source plane pixels
        ximg, yimg = np.meshgrid(np.linspace(xrimg[0], xrimg[1],
                                             nximg*bins, endpoint=False),
                                 np.linspace(yrimg[0], yrimg[1],
                                             nyimg*bins, endpoint=False),
                                 indexing='ij')

        xsrc  = ximg*(1.0-kappac-gamma) + defarr[:,:,0]
        ysrc  = yimg*(1.0-kappac+gamma) + defarr[:,:,1]
        ixsrc = (np.floor((xsrc-xmin)/dx).astype(int)).reshape(nximg*nyimg* 
                                                               bins*bins)
        iysrc = (np.floor((ysrc-ymin)/dy).astype(int)).reshape(nximg*nyimg* 
                                                               bins*bins)

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

        print "ML_MAGMAP: Number of rays shot: ", nximg*nyimg*bins*bins
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
        self.bins      = bins
        self.eps       = eps
        self.seed_pos  = seed_pos
        self.seed_rein = seed_rein
        self.nraymag1  = nraymag1
        self.stars     = stars
        self.cells     = cells
        self.defarr    = defarr
        self.magarr    = magarr
        self.rsrc      = 0.0
        self.stype     = ''
        self.kernel    = np.zeros([1,1])
        self.magcon    = magarr



    # -------- convolve the magnification map with a source morphology
    def convmap(self, rsrc, stype=None):

        # defaults
        if stype==None: stype='gaussian' # default Gaussian


        # Gaussian type
        if stype.lower()=='gaussian':

            # set the appropriate side length for the kernel
            rmax   = 2.0*rsrc*np.sqrt(np.log(10.)) # 99% containment
            pixsz  = (self.xr[1]-self.xr[0])/float(self.nxpix)
            side   = np.arange(-rmax,rmax,pixsz)

            # create the kernel
            xm, ym = np.meshgrid(side,side,indexing='ij')
            kernel = np.exp(-(xm**2+ym**2)/(2.0*rsrc**2))/(2.0*np.pi*rsrc**2)


        # do the convolution and set attributes
        self.rsrc    = rsrc
        self.stype   = stype
        self.kernel  = kernel
        self.magcon  = ndimage.convolve(self.magarr,kernel)
        self.magcon *= pixsz*pixsz # ndimage.convolve is sum not integral

        return

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



    # -------- generate a light curve for the magnification map
    def lightcurve(self, distance, angle, nsamp, conv=None, origin=None):

        """
          Samples the magnification map (or convolved map) from the 
          origin for a given distance and angle for nsamp points.  To 
          get actual physical values, the distances must be converted 
          to times via some velocity.

          [distance] = map units, [angle] = degrees, [origin] = pixels
        """

        # -------- set the map and origin
        nxpix = self.nxpix
        nypix = self.nypix
        x0    = self.xr[0]
        dx    = self.xr[1]-x0
        mmap  = self.magarr if conv==None else self.magcon

        if origin==None:
            origin = [i//2 for i in mmap.shape]
        else:
            origin = [int(round((float(i)-x0)/dx)) for i in origin]

        


        # -------- determine the x and y pixel values
        rad   = angle*np.pi/180.0
        crad  = np.cos(rad)
        srad  = np.sin(rad)
        pixsz = dx/float(nxpix)
        drpix = distance/dx * float(nxpix)/float(nsamp)
        ix    = [i*drpix*crad+origin[0] for i in range(nsamp)]
        iy    = [i*drpix*srad+origin[1] for i in range(nsamp)]


        # -------- check boundaries
        if (min(ix)<0) or (min(iy)<0) or (max(ix)>nxpix) or (max(iy)>nypix):
            print("ML_MAGMAP: Error - lightcurve too long!!!")
            sys.exit(-1)


        # -------- interpolate onto light curve pixels
        #          (using interp2d!!!)
        x, y   = np.arange(float(nxpix)), np.arange(float(nypix))
        mapint = interpolate.interp2d(x,y,mmap,kind='cubic')

        return np.array([mapint(i,j)[0] for i,j in zip(*[iy,ix])])

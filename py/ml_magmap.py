import sys
import datetime
import numpy as np
from scipy import ndimage, interpolate
import pyfits as fits
from .ml_gencells import *
from .ml_genstars import *
from .ml_counters import *
from .ml_defarr import *

class magmap():

    """
      magmap class includes the following data.
        lens system params:
          kappa - local convergence
          gamma - local shear
          fstar - fraction of surface mass density in stars
        source plane params:
          nxpix  - number of y pixels in the source plane
          nypix  - number of y pixels in the source plane
          xr     - [xmin,xmax] of the source plane
          yr     - [ymin,ymax] of the source plane
          magarr - magnification as a function of source position
          magcon - magarr convolved with kernel (source morphology)
        image plane params:
          nximg  - number of coarse x pixels in the image plane
          nyimg  - number of coarse y pixels in the image plane
          bins   - number of x and y subgrid bins
          defarr - deflection angle array (nximg*bins,nyimg*bins,2)
        star params:
          stars - the star field of type "starfield" class
          cells - list of tree cells each of type "treecell" class
        misc params:
          nraymag1 - number of rays per unity magnification per cell
        source params:
          rsrc   - size of the source
          stype  - source morphology
          kernel - convolution kernel (source morphology)
    """

    # -------- initialize the magnification map parameters
    def __init__(self, kappa, gamma, fstar, nxpix=1024, nypix=1024,
                 xr=None, yr=None, seed_pos=None, seed_rein=None,
                 nside=None, nximg=512, nyimg=512, multi=None,
                 bins=15, xsrc=None, ysrc=None, ixsrc=None,
                 iysrc=None, beta=None, mbar=None, mrat=None,
                 test=None):

        """ Initialize the magnification map parameters """

        # --------  defaults
        xr        = [-15., 15.] if xr==None else xr # in stellar Rein
        yr        = [-15., 15.] if yr==None else yr # in stellar Rein
        nside     = 16 if nside==None else nside # nside of stellar population
        seed_pos  = 111 if seed_pos==None else seed_pos
        seed_rein = 222 if seed_rein==None else seed_rein
        beta      = 0.0 if beta==None else beta
        mbar      = 1.0 if mbar==None else mbar
        mrat      = 1.0 if mrat==None else mrat

        if test!=None:
            nxpix     = 300
            nypix     = 300
            xr        = [-15., 15.]
            yr        = [-15., 15.]
            nside     = 16
            nximg     = 100
            nyimg     = 100
            bins      = 1
            seed_pos  = 111
            seed_rein = 222
            beta      = 0.0
            mbar      = 1.0
            mrat      = 1.0


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

        xrimg     = 1.05*xr / abs(1.0 - kappa - gamma) # emperical factor
        yrimg     = 1.05*yr / abs(1.0 - kappa + gamma) # emperical factor
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
        print("ML_MAGMAP: Shooting region size: {0} by {1}"
              .format(dximg*nximg*bins, dyimg*nyimg*bins))
        print("ML_MAGMAP: Number of rays per pixel for unity " + 
              "magnification: {0}".format(nraymag1))
        print("ML_MAGMAP: Size of stars region: {0} {1}".format(xrst, yrst))


        # -------- initialize the stars and cells
        stars = genstars(kappas, xrst, yrst, nside, seed_pos,
                         seed_rein, beta=beta, mbar=mbar, mrat=mrat)
        cells  = gencells(stars)
        nstars = np.array([i.nstar for i in cells])

        counters.maxnst = max(nstars[np.array([i.high for i in cells])==1])

        print("ML_MAGMAP: Mean number of stars per highest cell: {0}"
              .format(np.mean(nstars[np.array([i.high for i in cells])==1])))


        # -------- initialize the deflection angle map
        xpos   = xmin + dx*np.arange(nxpix)
        ypos   = ymin + dy*np.arange(nypix)


        # -------- shoot the rays in the image plane (calculate the
        #          deflection angle at each point
        defarr = ml_defarr(xrimg, yrimg, nximg, nyimg, cells, stars, \
                           multi=multi, bins=bins)

        if (len(defarr)==1) and (defarr==-1):
            print("ML_MAGMAP: ERROR: # of processors (multi kw) ")
            print("ML_MAGMAP:   inappropriately set... aborting.")
            return


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
            print("ML_MAGMAP: NO RAYS DEFLECTED INTO SOURCE PLANE GRID!!!")
            return

        isrc_ring = np.bincount(ixsrc[w] + iysrc[w]*nxpix)        


        # -------- convert to magnification
        magarr = np.zeros(nxpix*nypix)
        magarr[:len(isrc_ring)] += isrc_ring
        magarr = magarr.reshape(nxpix,nypix)

        print("ML_MAGMAP: Number of rays shot: {0}"
              .format(nximg*nyimg*bins*bins))
        print("ML_MAGMAP: Number of rays landed: {0}"
              .format(w.size))

        magarr /= nraymag1


        # -------- set attributes
        self.names = ['kappa', 'gamma', 'fstar', 'nxpix', 'nypix',
                      'xr', 'yr', 'nximg', 'nyimg', 'bins',
                      'nraymag1', 'stars', 'cells', 'defarr',
                      'magarr', 'rsrc', 'stype', 'kernel', 'magcon']

        self.kappa     = kappa
        self.gamma     = gamma
        self.fstar     = fstar
        self.nxpix     = nxpix
        self.nypix     = nypix
        self.xr        = xr
        self.yr        = yr
        self.nximg     = nximg
        self.nyimg     = nyimg
        self.bins      = bins
        self.nraymag1  = nraymag1
        self.stars     = stars
        self.cells     = cells
        self.defarr    = defarr
        self.magarr    = magarr
        self.rsrc      = 0.0
        self.stype     = ''
        self.kernel    = np.zeros([1,1])
        self.magcon    = magarr



    # -------- write magnification map and important parameters to fits
    def write_fits(self, filename, path=None):

        """ Write magnification map to fits file """

        # define output file
        path = "" if path==None else path if path[-1]=="/" else path + "/"
        fout = path + filename

        # set up HDU
        hdu = fits.PrimaryHDU(self.magarr)

        # add creation info to header
        now = datetime.date.today().strftime("%B %d, %Y")
        hdu.header['SIMPLE'] = (True, 'created by MULE via PyFITS on ' + now)

        # add important info to header
        keys = ['kappa', 'gamma', 'fstar', 'xmin', 'xmax', 'ymin',
                'ymax', 'nraymag1']
        vals = [self.kappa, self.gamma, self.fstar, self.xr[0],
                self.xr[1], self.yr[0], self.yr[1], self.nraymag1]
        coms = ['local convergence', 'local shear', 
                'fraction of surface density in stars', 
                'min source plane x-coord', 'max source plane x-coord', 
                'min source plane y-coord', 'max source plane y-coord', 
                '# rays per pixel for unity magnification']

        for i in range(len(keys)):
            hdu.header[keys[i]] = (vals[i], coms[i])

        # write to file
        hdu.writeto(fout, clobber=True)

        return



    # -------- convolve the magnification map with a source morphology
    def convmap(self, rsrc, stype=None):

        """ Convolve the magnification map with a source morphology """

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

          [distance] = map units, [angle] = degrees, [origin] = map units
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
            origin = [int(round((i-x0)*float(nxpix)/dx)) for i in origin]


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



    # -------- create a visualization of the map and light curve
    def view(self, fnum=None, **kwargs):

        """ Plot visualization of the magnification map and lightcurve. """





    # -------- return an item by its name
    def __getitem__(self, key):

        """ Return an item by its name."""

        data = { 
            'nxpix'     : self.nxpix,
            'nypix'     : self.nypix,
            'xr'        : self.xr,
            'yr'        : self.yr,
            'nximg'     : self.nximg,
            'nyimg'     : self.nyimg,
            'bins'      : self.bins,
            'nraymag1'  : self.nraymag1,
            'stars'     : self.stars,
            'cells'     : self.cells,
            'defarr'    : self.defarr,
            'magarr'    : self.magarr,
            'rsrc'      : self.rsrc,
            'stype'     : self.stype,
            'kernel'    : self.kernel,
            'magcon'    : self.magcon
        }

        return data[key]

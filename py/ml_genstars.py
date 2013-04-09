import numpy as np
from ml_cell_ring2nest import *

class starfield():

    """
      starfield class includes the following data.
        star parameters    : seed_pos, seed_rein, xstar, ystar, rein, icell
        field parameters   : nstars, xrange_st, yrange_st
        binning parameters : nside
    """

# -------- initialize the star field
    def __init__(self, nstars, xrange_st, yrange_st, nside, seed_pos, \
                     seed_rein, xstar, ystar, rein, icell):

        """ Initialize the star field parameters """

        # initialize data type and names
        self.dtype = 'starfield'
        self.names = ['seed_pos', 'seed_rein', 'xstar', 'ystar', 'rein', \
                          'icell', 'nstars', 'xrange_st', 'yrange_st', 'nside']


        # initialize star field parameters and seeds
        self.nstars    = nstars
        self.xrange_st = xrange_st
        self.yrange_st = yrange_st
        self.nside     = nside
        self.seed_pos  = seed_pos
        self.seed_rein = seed_rein


        # initialize star parameters
        self.xstar = xstar
        self.ystar = ystar
        self.rein  = rein
        self.icell = icell



# -------- return an item by its name
    def __getitem__(self, key):

        """ Return an item by its name. """

        data = { \
                'nstars'    : self.nstars, \
                'xrange_st' : self.xrange_st, \
                'yrange_st' : self.yrange_st, \
                'nside'     : self.nside, \
                'seed_pos'  : self.seed_pos, \
                'seed_rein' : self.seed_rein, \
                'xstar'     : self.xstar, \
                'ystar'     : self.ystar, \
                'rein'      : self.rein, \
                'icell'     : self.icell, \
                }

        return data[key]





def genstars(kappas, xrange_st, yrange_st, nside, seed_pos, seed_rein):

    """
    NAME:

    PURPOSE:

    CALLING SEQUENCE:

    INPUTS:

    OPTIONAL INPUTS:

    KEYWORDS:

    OUTPUTS:

    OPTIONAL OUTPUTS:

    EXAMPLES:

    COMMENTS:

    REVISION HISTORY:
      2013/03/26 - Written by Greg Dobler (KITP/UCSB)

    ------------------------------------------------------------
    """

# -------- make sure that inputs are of the correct type
    kappas    = float(kappas)
    xrange_st = np.array(xrange_st).astype(float)
    yrange_st = np.array(yrange_st).astype(float)
    nside     = int(nside)
    seed_pos  = int(seed_pos)
    seed_rein = int(seed_rein)



# -------- utilities
    delx   = xrange_st[1]-xrange_st[0]
    dely   = yrange_st[1]-yrange_st[0]
    nstars = int(round(kappas*delx*dely/np.pi))

    print "ML_GENSTARS: Creating star field with {0} stars".format(nstars)



# -------- draw positions and einstein radii
    np.random.seed(seed_pos)

    xstar = xrange_st[0] + delx*np.random.rand(nstars)
    ystar = yrange_st[0] + dely*np.random.rand(nstars)
    rein  = np.ones(nstars)



# -------- sort stars into ring-ordered cells and convert to nest
    dx         = delx/np.float(nside)
    dy         = dely/np.float(nside)
    xcell      = np.floor((xstar-xrange_st[0])/dx).astype(int)
    ycell      = np.floor((ystar-yrange_st[0])/dy).astype(int)
    icell_ring = xcell + nside*ycell
    icell      = ml_cell_ring2nest(icell_ring, nside=nside, quiet=True)



# -------- sort according to cell number
    ord = np.argsort(icell)

    xstar = xstar[ord]
    ystar = ystar[ord]
    rein  = rein[ord]
    icell = icell[ord]



# -------- put into star field class and return
    stars = starfield(nstars, xrange_st, yrange_st, nside, seed_pos, \
                          seed_rein, xstar, ystar, rein, icell)

    return stars

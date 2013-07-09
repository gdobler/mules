import numpy as np
from ml_cell_ring2nest import *

class starfield():

    """
      starfield class includes the following data.
        star parameters : seed_pos, seed_rein, xstar, ystar, rein, icell
        mass function parameters : beta, mbar, mrat
        field parameters : nstars, xrange_st, yrange_st
        binning parameters : nside
    """

    # -------- initialize the star field
    def __init__(self, nstars, xrange_st, yrange_st, nside, seed_pos,
                 seed_rein, xstar, ystar, rein, icell, beta, mbar,
                 mrat):

        """ Initialize the star field parameters """

        # initialize data type and names
        self.dtype = 'starfield'
        self.names = ['seed_pos', 'seed_rein', 'xstar', 'ystar',
                      'rein', 'icell', 'beta', 'mbar', 'mrat',
                      'nstars', 'xrange_st', 'yrange_st', 'nside']


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

        # initialize mass function parameters
        self.beta  = beta
        self.mbar  = mbar
        self.mrat  = mrat



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
                 'beta'      : self.beta, \
                 'mbar'      : self.mbar, \
                 'mrat'      : self.mrat \
                }

        return data[key]



def mf_draw(seed_rein, nstars, beta, mbar, mrat):

    """
      2013/06/05 - Written by Greg Dobler (KITP/UCSB)
      Draw Einstein radii (m^1/2) from a powerlaw mass function of
      stars with masses in units of solar masses.  E.g., m = 4 M_sun 
      implies Rein = 2.
    """

    # -------- check for flat mass function
    if (abs(beta) < 1.0e-5) or (abs(mrat)-1.0 < 1.0e-5):
        return mbar*np.ones(nstars)


    # -------- utilities
    bp1     = beta + 1.0
    bp2     = beta + 2.0
    rtbp1m1 = mrat**(beta+1) - 1.0
    rtbp2m1 = (rtbp1m1+1.0)*mrat - 1.0
    oobp1   = 1.0/bp1

    np.random.seed(seed_rein)
    ran = np.random.rand(nstars)


    # -------- draw from mass function
    mmin = mbar * bp2/bp1 * rtbp1m1/rtbp2m1
    m    = mmin * (ran*rtbp1m1 + 1.0)**oobp1

    return np.sqrt(m)



def genstars(kappas, xrange_st, yrange_st, nside, seed_pos, seed_rein,
             beta=None, mbar=None, mrat=None):

    """
    NAME:
      genstars

    PURPOSE:
      Generate a population of stars.

    CALLING SEQUENCE:
      stars = genstars(kappas, xrange_st, yrange_st, nside, seed_pos, 
                       seed_rein, [beta=, mbar=, mrat=])

    INPUTS:
      kappas    - surface density in stars
      xrange_st - the size of the stellar field in the x direction
      yrange_st - the size of the stellar field in the y direction
      nside     - number of cells per side
      seed_pos  - seed for the positions of the stars
      seed_rein - seed for the einstein radii of the stars

    OPTIONAL INPUTS:
      beta - power law in the stellar mass function (default 0.0)
      mbar - mean mass of the stars (default 1.0)
      mrat - ratio of lowest to highest mass in mass function (default 1.0)

    KEYWORDS:

    OUTPUTS:
      stars - star field of type "starfield" class"

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
    beta      = 0.0 if beta==None else float(beta)
    mbar      = 1.0 if mbar==None else float(mbar)
    mrat      = 1.0 if mrat==None else float(mrat)



    # -------- utilities
    delx   = xrange_st[1]-xrange_st[0]
    dely   = yrange_st[1]-yrange_st[0]
    nstars = int(round(kappas*delx*dely/np.pi))

    print("ML_GENSTARS: Creating star field with {0} stars".format(nstars))



    # -------- draw positions and einstein radii
    np.random.seed(seed_pos)

    xstar = xrange_st[0] + delx*np.random.rand(nstars)
    ystar = yrange_st[0] + dely*np.random.rand(nstars)
    rein  = mf_draw(seed_rein, nstars, beta, mbar, mrat)



    # -------- sort stars into ring-ordered cells and convert to nest
    dx         = delx/np.float(nside)
    dy         = dely/np.float(nside)
    xcell      = np.floor((xstar-xrange_st[0])/dx).astype(int)
    ycell      = np.floor((ystar-yrange_st[0])/dy).astype(int)
    icell_ring = xcell + nside*ycell
    icell      = ml_cell_ring2nest(icell_ring, nside=nside, quiet=True)



    # -------- sort according to cell number
    srt = np.argsort(icell)

    xstar = xstar[srt]
    ystar = ystar[srt]
    rein  = rein[srt]
    icell = icell[srt]



    # -------- put into star field class and return
    stars = starfield(nstars, xrange_st, yrange_st, nside, seed_pos,
                      seed_rein, xstar, ystar, rein, icell, beta,
                      mbar, mrat)

    return stars

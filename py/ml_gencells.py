import numpy as np
from ml_cellpos import *

class treecell():
    """
      treecell class includes the following data.
        ??? : ???
        ??? : ???
        ??? : ???
    """

    def __init__(self, icell, nstar, xcell, ycell, mcell, xside,
                 yside, cmom, smom, high, stind):

        self.dtype = 'treecell'
        self.names =['icell', 'nstar', 'xcell', 'ycell', 'mcell',
                     'xside', 'yside', 'cmom', 'smom', 'high',
                     'stind']

        self.icell = icell
        self.nstar = nstar
        self.xcell = xcell
        self.ycell = ycell
        self.mcell = mcell
        self.xside = xside
        self.yside = yside
        self.cmom  = cmom
        self.smom  = smom
        self.high  = high
        self.stind = stind


# -------- return an item by its name
    def __getitem__(self, key):

        """ Return an item by its name. """

        data = { \
                'icell' : self.icell, \
                'nstar' : self.nstar, \
                'xcell' : self.xcell, \
                'ycell' : self.ycell, \
                'mcell' : self.mcell, \
                'xside' : self.xside, \
                'yside' : self.yside, \
                'cmom'  : self.cmom, \
                'smom'  : self.smom, \
                'high'  : self.high, \
                'stind' : self.stind, \
                }

        return data[key]



def gencells(stars):

# -------- utilities
    xstar = stars.xstar
    ystar = stars.ystar
    rein  = stars.rein
    nside = stars.nside
    xrst  = stars.xrange_st
    yrst  = stars.yrange_st
    cpos  = ml_cellpos(xrst, yrst, nside)
    nlev  = np.log2(nside).astype(int)
    cmin  = np.sum(4**np.arange(nlev)).astype(int) - 1
    ncell = nside**2



# -------- find the maximum number of stars in a cell
    sicell = stars.icell
    maxst  = (np.histogram(sicell,bins=range(sicell.max())))[0].max()



# --------  loop through the cells and claulate the multipole moments
    for ilev in range(nlev):
        print("ML_GENCELLS: Building level {0} of {1}".format(ilev,nlev-1))

        tncell = 4**(ilev+1)
        tcmin  = np.sum(4**np.arange(ilev+1))-1
        tcmax  = tcmin + tncell
        dcell  = ncell/tncell

        txside = (xrst[1]-xrst[0])/2.0**(ilev+1)
        tyside = (xrst[1]-xrst[0])/2.0**(ilev+1)

        for icell in range(tncell):
            tcell  = tcmin + icell
            stind  = np.where((sicell >= cmin + icell*dcell) & \
                                  (sicell < cmin + (icell+1)*dcell))[0]
            tnstar = len(stind)

            if tnstar > 0:
                tposx = cpos[0,tcell] # mean x position of cell
                tposy = cpos[1,tcell] # mean y position of cell
                tmass = np.sum(rein[stind]*rein[stind]) # total cell mass

                delxi = tposx - xstar[stind]
                delyi = tposy - ystar[stind]
                ri    = np.sqrt(delxi*delxi + delyi*delyi)
                ti    = np.arctan2(delyi, delxi)
                mi    = stars.rein[stind]

                nmm = 20

                for m in range(1,nmm+1):
                    cc = np.sum(mi*mi * ri**float(m) * np.cos(m*ti))
                    ss = np.sum(mi*mi * ri**float(m) * np.sin(m*ti))

                    tcmom = cc if m==1 else np.append(tcmom,cc)
                    tsmom = ss if m==1 else np.append(tsmom,ss)

                thigh  = int(tcell >= cmin) # flag highest cell
                tstind = -np.ones(maxst,dtype=int)

                if thigh:
                    tstind[0:tnstar] = stind # inds of stars in cell
            else:
                tposx  = cpos[0,tcell] # mean x position of cell
                tposy  = cpos[1,tcell] # mean y position of cell
                tmass  = 0.0
                tcmom  = np.zeros(20)
                tsmom  = np.zeros(20)
                thigh  = int(tcell >= cmin) # flag highest cell
                tstind = -np.ones(maxst,dtype=int)


            # pack into a list of tree cell classes
            ttree = treecell(tcell, tnstar, tposx, tposy, tmass,
                             txside, tyside, tcmom, tsmom, thigh,
                             tstind)

            if (ilev==0) & (icell==0):
                cells = [ttree] 
            else:
                cells.append(ttree)



# --------  return cells list
    return cells

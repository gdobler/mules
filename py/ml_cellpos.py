import numpy as np
from .ml_cell_ring2nest import *

def ml_cellpos(xr, yr, nside):

    """
    NAME:
      ml_cellpos

    PURPOSE:
      Return the central positions of cartesian cells for a given xrange, 
      yrange, and nside.

    CALLING SEQUENCE:
      cpos = ml_cellpos(xr, yrnage, nside)

    INPUTS:
      xr     - the x range of the grid
      yr     - the y range of the grid
      nside  - the number of cells per side (total # of cells = nside^2)

    OPTIONAL INPUTS:

    KEYWORDS:

    OUTPUTS:
      cpos - central positions of cells on the grid.

    OPTIONAL OUTPUTS:

    EXAMPLES:

    COMMENTS:

    REVISION HISTORY:
      2013/03/26 - Written by Greg Dobler (KITP/UCSB)

    ------------------------------------------------------------
    """

    # -------- utilities
    nlev = np.log2(nside).astype(int)
    xmin, xmax = xr[0], xr[1]
    ymin, ymax = yr[0], yr[1]



    # -------- loop through the levels
    for ilev in range(nlev):
        tside = 2**(ilev+1)
        dx    = (xmax-xmin)/np.float(tside)
        dy    = (ymax-ymin)/np.float(tside)

        # get ring ordered cell numbers and x and y positions
        icell_ring = np.arange(tside**2)
        xpos       = xmin + dx*((icell_ring % tside) + 0.5)
        ypos       = ymin + dy*((icell_ring / tside) + 0.5)

        # convert ring cell numbers to nested cell numbers
        icell = ml_cell_ring2nest(icell_ring, quiet=True)

        # stack the cell numbers together
        sicell = icell if ilev==0 else np.concatenate([sicell,icell])
        sxpos  = xpos  if ilev==0 else np.concatenate([sxpos,xpos])
        sypos  = ypos  if ilev==0 else np.concatenate([sypos,ypos])



    # -------- order by cell number, concatenate, and return
    srt  = np.argsort(sicell)
    cpos = np.concatenate([[sxpos[srt]],[sypos[srt]]])

    return cpos

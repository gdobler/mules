import numpy as np
from ml_cell_ring2nest import *

def ml_cellpos(xrange, yrange, nside):

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

# -------- utilities
    nlev = np.log2(nside).astype(int)
    xmin, xmax = xrange[0], xrange[1]
    ymin, ymax = yrange[0], yrange[1]



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
    ord  = np.argsort(sicell)
    cpos = np.concatenate([[sxpos[ord]],[sypos[ord]]])

    return cpos

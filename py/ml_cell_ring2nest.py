import numpy as np

def ml_cell_ring2nest(cell_ring, nside=None, quiet=None):

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
    nside = nside if nside else np.sqrt(cell_ring.size).astype(int)
    xring = cell_ring % nside
    yring = cell_ring / nside

    if quiet==None:
        print("ML_CELL_RING2NEST: converting cells with nside = {0}"
              .format(nside))



# -------- convert ring cell numbers to nested cell numbers and return
    nnum = np.log2(nside).astype(int) - 1
    offr = 0 if nnum==0 else  np.cumsum(2.0**(2.0*np.arange(nnum)+1))
    offc = 0 if nnum==0 else -np.cumsum(2.0**(2.0*np.arange(nnum)+2))-2.0-nside

    row    = np.zeros(nside).astype(int)
    col    = np.zeros(nside).astype(int) - 2 - nside
    col[0] = 0 # special case

    for ioff in range(nnum):
        tind = 2**(ioff+2)*np.arange(nside/(2**(ioff+2))) + 2**(ioff+1)

        row[tind] = offr[ioff]
        col[tind] = offc[ioff]

    arow      = np.cumsum(row).astype(int)
    acol      = np.cumsum(col).astype(int)
    cell_nest = cell_ring + acol[yring] + arow[xring] + nside*nside - 2L

    return cell_nest

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
#  2012/01/19 - Written by Greg Dobler (KITP/UCSB)
#
#------------------------------------------------------------
def ml_cellpos(xrange, yrange, nside):

# -------- utilities
    from math import *
    from numpy import *

    nlev = round(log(float(nside))/alog(2.0))
    xmin = xrange[0]
    xmax = xrange[1]
    dx   = xmax - xmin
    ymin = yrange[0]
    ymax = yrange[1]
    dy   = ymax - ymin



# -------- loop through the levels
    for ilev in range(0L,nlev-1):
        tside = 2L^(ilev+1)

        # get ring ordered cell numbers and x and y positions
        icell_ring = arange(0L,tside**2,1)
        xring      = icell_ring % tside
        yring      = icell_ring / tside
        xpos       = xmin + dx*(xring+0.5)
        ypos       = ymin + dy*(yring+0.5)

        # convert ring cell numbers to nested cell numbers
        nnum = long(round(log(float(tside))/log(2.0))) - 1L
        offr = (cumsum(2.0**(2.0*arange(0,nnum,1,dtype=float)+1)) if nnum > 0 
                else 0)
        offc = (-cumsum(2.0**(2.0*(arange(0,nnum,1,dtype=float)+1)))-2.0-tisde 
                 if nnum > 0 else 0)
        row  = zeros(tside,dtype=int)
        col  = zeros(tside,dtype=int) - 2L - tside

        for ioff in range(0L,nnum-1):
            tind = 2L**(ioff+2)*arange(0,tside/(2L**(ioff+2))) + 2L**(ioff+1)
            print tind

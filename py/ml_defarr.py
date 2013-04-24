from ml_defang import *
import sys

def ml_defarr(defarr, xr, yr, nx, ny, cells, stars, scheme=None, bin=None):

# -------- defaults
    scheme = 1 if scheme==None else scheme



# -------- scheme 1: loop directly
    if scheme==1:
        print "ML_DEFARR: looping directly through 2D image plane..."

        for iximg in range(nx):
            print 'ML_DEFARR:   iximg = {0} out of {1}\r'.format(iximg+1,nx), 
            sys.stdout.flush()

            for iyimg in range(ny):
                ximg = xr[0] + (xr[1]-xr[0])*float(iximg)/float(nx)
                yimg = yr[0] + (yr[1]-yr[0])*float(iyimg)/float(ny)

                defarr[iximg,iyimg,:] = ml_defang(ximg,yimg,cells,stars,None)

        print
        return

# -------- linearly interpolate in fixed cells
    elif scheme==2:

        # -------- defaults and utilities
        if bin==None: bin = 64
        defur  = np.zeros(2)
        deful  = np.zeros(2)
        deflr  = np.zeros(2)
        defll  = np.zeros(2)
        bp1    = bin + 1                  #
        xind   = np.arange(bp1*bp1) % bp1 #
        yind   = np.arange(bp1*bp1) / bp1 # utilities to be used 
        xfrac  = xind/float(bin)          # for the interpolation
        yfrac  = yind/float(bin)          #
        ncx    = (nx-1)/bin
        ncy    = (ny-1)/bin
        x0, x1 = xr
        y0, y1 = yr
        ddx    = (x1-x0)/float(nx)
        ddy    = (y1-y0)/float(ny)

        print "ML_DEFARR: linearly interpolating 2D img plane with bin= ", bin

        for icx in range(ncx):
            print "ML_DEFARR:   icx = {0} out of {1}\r".format(icx+1,ncx), 
            sys.stdout.flush()

            # -------- get four corners in x
            iximgul = icx*bin
            iximgll = iximgul
            iximgur = iximgul+bin
            iximglr = iximgll+bin

            ximgul = x0 + ddx*float(iximgul)
            ximgur = x0 + ddx*float(iximgur)
            ximgll = x0 + ddx*float(iximgll)
            ximglr = x0 + ddx*float(iximglr)

            for icy in range(ncy):

                # -------- get four corners in y
                iyimgll = icy*bin
                iyimglr = icy*bin
                iyimgul = iyimgll+bin
                iyimgur = iyimglr+bin

                yimgul = y0 + ddy*float(iyimgul)
                yimgur = y0 + ddy*float(iyimgur)
                yimgll = y0 + ddy*float(iyimgll)
                yimglr = y0 + ddy*float(iyimglr)

                # -------- calculate the deflection angles, only need
                #          to do all the corners for special cases
                defur[:] = ml_defang(ximgur, yimgur, cells, stars, None)

                if icx==0:
                    deful[:] = ml_defang(ximgul, yimgul, cells, stars, None)
                else:
                    deful[:] = defarr[iximgul,iyimgul,:]

                if icy==0:
                    deflr[:] = ml_defang(ximglr, yimglr, cells, stars, None)
                else:
                    deflr[:] = defarr[iximglr,iyimglr,:]

                if (icx==0) & (icy==0):
                    defll[:] = ml_defang(ximgll, yimgll, cells, stars, None)
                else:
                    defll[:] = defarr[iximgll,iyimgll,:]


                # -------- bi-linearly interpolate from the four corners
                defarr[xind+iximgll,yind+iyimgll,0] = defll[0]*(1.0-xfrac)* \
                    (1.0-yfrac) + deflr[0]*xfrac*(1.0-yfrac)+ \
                    deful[0]*(1.0-xfrac)*yfrac + defur[0]*xfrac*yfrac
                defarr[xind+iximgll,yind+iyimgll,1] = defll[1]*(1.0-xfrac)* \
                    (1.0-yfrac) + deflr[1]*xfrac*(1.0-yfrac)+ \
                    deful[1]*(1.0-xfrac)*yfrac + defur[1]*xfrac*yfrac

        print
        return

    else:
        print "ML_DEFARR: ONLY SCHEMES #1 and #2 ARE IMPLEMENTED!!!"
        return

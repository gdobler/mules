from scipy import interpolate
from ml_defang_iter import *
from ml_defang import *
import sys
import multiprocessing

def ml_defarr(xr, yr, nx, ny, cells, stars, multi=None, recur=None,
              bins=None):

    # -------- defaults
    multi    = 4    if multi==None else multi
    recur    = True if recur==None else recur
    bins     = 0    if bins==None  else bins
    def_func = ml_defang if recur else ml_defang_iter


    # -------- utilities
    if multi:
        if ny % multi != 0:
            print("ML_DEFARR: ERROR: Number of image plane y-pixels must be")
            print("ML_DEFARR:  an integer number of processors!!!")
            print("ML_DEFARR:  (i.e., nximg mod multi = 0)")

            return -1

        nproc = multi
        dy    = (yr[1]-yr[0])/float(nproc)
        yrs   = [[i*dy+yr[0],(i+1)*dy+yr[0]] for i in range(nproc)]


    # -------- loop through the image plane
    def run_rect(conn,xr,yr,nx,ny,cells,stars,verbose=False):
        rect = np.zeros([nx,ny,2], dtype=np.float64)
        ximg = np.linspace(xr[0],xr[1],nx,endpoint=False)
        yimg = np.linspace(yr[0],yr[1],ny,endpoint=False)

        for i in range(nx):
            if verbose and ((i+1) % 5 == 0):
                print(('ML_DEFARR:     iximg = {0} out ' + 
                       'of {1}\r').format(i+1, nx)), 
                sys.stdout.flush()

            for j in range(ny):
                rect[i,j,:] = def_func(ximg[i],yimg[j],cells,stars,None)

        if verbose: print
        if multi:
            conn.send(rect)
            conn.close()
        else:
            return rect


    # -------- announce the method
    print("ML_DEFARR: looping through the 2D image plane with:")
    print("ML_DEFARR:   multiprocessing      : " + ("ON" if multi  else "OFF"))
    print("ML_DEFARR:   recursion            : " + ("ON" if recur  else "OFF"))
    print("ML_DEFARR:   subgrid inerpolation : " + ("ON" if bins>0 else "OFF"))


    # -------- calculate the deflection angle
    if multi:
        print("ML_DEFARR:   running {0} processes...".format(nproc))

        # -------- initialize the deflection array and processes
        defarr  = np.zeros([nx,ny,2], dtype=np.float64)
        parents, childs, ps = [], [], []

        # -------- initialize the pipes and processes, then start
        for ip in range(nproc):
            ptemp, ctemp = multiprocessing.Pipe()
            parents.append(ptemp)
            childs.append(ctemp)
            ps.append(multiprocessing.Process(target=run_rect,
                                              args=(childs[ip],xr,yrs[ip],
                                                    nx,ny/nproc,cells,stars),
                                              kwargs={'verbose':ip==0}))
            ps[ip].start()

        # -------- collect the results, put into defarr, and rejoin
        for ip in range(nproc):
            yind = [ny/nproc*ip,ny/nproc*(ip+1)]
            defarr[0:nx,yind[0]:yind[1],:] = parents[ip].recv()
            ps[ip].join()
    else:
        defarr = run_rect(-314,xr,yr,nx,ny,cells,stars,verbose=True)

    return defarr

"""
# -------- linearly interpolate in fixed cells
    elif scheme==2:

        # -------- defaults and utilities
        if bins==None: bins = 64
        defur  = np.zeros(2)
        deful  = np.zeros(2)
        deflr  = np.zeros(2)
        defll  = np.zeros(2)
        bp1    = bins + 1                  #
        xind   = np.arange(bp1*bp1) % bp1 #
        yind   = np.arange(bp1*bp1) / bp1 # utilities to be used 
        xfrac  = xind/float(bins)          # for the interpolation
        yfrac  = yind/float(bins)          #
        ncx    = (nx-1)/bins
        ncy    = (ny-1)/bins
        x0, x1 = xr
        y0, y1 = yr
        ddx    = (x1-x0)/float(nx)
        ddy    = (y1-y0)/float(ny)

        print "ML_DEFARR: linearly interpolating 2D img plane with bins= ", \
            bins

        for icx in range(ncx):
            print "ML_DEFARR:   icx = {0} out of {1}\r".format(icx+1,ncx), 
            sys.stdout.flush()

            # -------- get four corners in x
            iximgul = icx*bins
            iximgll = iximgul
            iximgur = iximgul+bins
            iximglr = iximgll+bins

            ximgul = x0 + ddx*float(iximgul)
            ximgur = x0 + ddx*float(iximgur)
            ximgll = x0 + ddx*float(iximgll)
            ximglr = x0 + ddx*float(iximglr)

            for icy in range(ncy):

                # -------- get four corners in y
                iyimgll = icy*bins
                iyimglr = icy*bins
                iyimgul = iyimgll+bins
                iyimgur = iyimglr+bins

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

# -------- evaluate grid and then spline
    elif scheme==3:

        print "ML_DEFARR: spline interpolating 2D grid..."

        # -------- defaults and utilities
        if bins==None: bins = 64

        cnx    = nx/bins
        cny    = ny/bins
        coarse = np.zeros([cnx,cny,2])
        xcoar  = np.linspace(xr[0],xr[1],num=cnx,endpoint=False)
        ycoar  = np.linspace(yr[0],yr[1],num=cny,endpoint=False)

        for iximg in range(cnx):
            print 'ML_DEFARR:   iximg = {0} out of {1}\r'.format(iximg+1,cnx), 
            sys.stdout.flush()

            for iyimg in range(cny):
                coarse[iximg,iyimg,:] = ml_defang(xcoar[iximg],ycoar[iyimg], 
                                                  cells,stars,None)
        print


        # -------- interpolate from the coarse grid onto the real grid
        ximg = np.linspace(xr[0],xr[1],num=nx,endpoint=False)
        yimg = np.linspace(yr[0],yr[1],num=ny,endpoint=False)

        defxint = interpolate.RectBivariateSpline(xcoar,ycoar,coarse[:,:,0])
        defyint = interpolate.RectBivariateSpline(xcoar,ycoar,coarse[:,:,1])

        defarr[:,:,0] = defxint(ximg,yimg)
        defarr[:,:,1] = defyint(ximg,yimg)

        return

    else:
        print "ML_DEFARR: ONLY SCHEMES #1, 2, and 3 ARE IMPLEMENTED!!!"
        return
"""

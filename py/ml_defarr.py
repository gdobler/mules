from scipy import interpolate
from ml_defang_iter import *
import sys
import multiprocessing

def ml_defarr(defarr, xr, yr, nx, ny, cells, stars, scheme=None, bins=None):

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

                defarr[iximg,iyimg,:] = ml_defang_iter(ximg,yimg,cells,stars,None)

        print
        return


    # -------- scheme 1: loop directly (threaded)
    elif scheme==-1:
        print "ML_DEFARR: looping directly through 2D image plane..."

        # -------- dummy function to thread
        def run_quadrant(conn,xr,yr,nx,ny,cells,stars,verbose=False):
            quad = np.zeros([nx,ny,2], dtype=np.float64)
            ximg = np.linspace(xr[0],xr[1],nx,endpoint=False)
            yimg = np.linspace(yr[0],yr[1],ny,endpoint=False)

            for i in range(nx):
                if verbose:
                    print(('ML_DEFARR:   iximg = {0} out ' + 
                           'of {1}\r').format(i+1, nx)), 
                    sys.stdout.flush()

                for j in range(ny):
                    quad[i,j,:] = ml_defang_iter(ximg[i],yimg[j],cells,stars,
                                                 None)

            if verbose: print
            conn.send(quad)
            conn.close()

        # -------- define the ranges and dummy maps
        xmin, xmax, xmid = xr[0], xr[1], 0.5*(xr[0]+xr[1])
        ymin, ymax, ymid = yr[0], yr[1], 0.5*(yr[0]+yr[1])

        xrul, yrul = [xmin,xmid], [ymid,ymax]
        xrur, yrur = [xmid,xmax], [ymid,ymax]
        xrll, yrll = [xmin,xmid], [ymin,ymid]
        xrlr, yrlr = [xmid,xmax], [ymin,ymid]

        # -------- initialize results queue
        parent_conn1, child_conn1 = multiprocessing.Pipe()
        parent_conn2, child_conn2 = multiprocessing.Pipe()
        parent_conn3, child_conn3 = multiprocessing.Pipe()
        parent_conn4, child_conn4 = multiprocessing.Pipe()

        # -------- initialize the threads (only print for first quadrant)
        threadul = multiprocessing.Process(target=run_quadrant, 
                                   args=(child_conn1,xrul,yrul,nx/2l,ny/2l, 
                                         cells,stars), 
                                   kwargs={'verbose':True})
        threadur = multiprocessing.Process(target=run_quadrant, 
                                   args=(child_conn2,xrur,yrur,nx/2l,ny/2l,
                                         cells,stars))
        threadll = multiprocessing.Process(target=run_quadrant, 
                                   args=(child_conn3,xrll,yrll,nx/2l,ny/2l,
                                         cells,stars))
        threadlr = multiprocessing.Process(target=run_quadrant, 
                                   args=(child_conn4,xrlr,yrlr,nx/2l,ny/2l,
                                         cells,stars))

        # -------- run the threads
        print "ML_DEFARR: running 4 threads..."
        threadul.start()
        threadur.start()
        threadll.start()
        threadlr.start()

        # -------- put the threads into the defarr
        defarr[0:nx/2,  ny/2:ny, :] = parent_conn1.recv()
        defarr[nx/2:nx, ny/2:ny, :] = parent_conn2.recv()
        defarr[0:nx/2,  0:ny/2,  :] = parent_conn3.recv()
        defarr[nx/2:nx, 0:ny/2,  :] = parent_conn4.recv()

        # -------- allow threads to complete
        threadul.join()
        print "ML_DEFARR: upper left  thread has joined..."
        threadur.join()
        print "ML_DEFARR: upper right thread has joined..."
        threadll.join()
        print "ML_DEFARR: lower left  thread has joined..."
        threadlr.join()
        print "ML_DEFARR: lower right thread has joined..."

        return

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

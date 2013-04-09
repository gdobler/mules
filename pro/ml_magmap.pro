;+
; NAME:
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;  2011/07/22 - Written by Greg Dobler (KITP/UCSB)
;
;------------------------------------------------------------
function ml_magmap, kappa, gamma, fstar, nxpix=nxpix, nypix=nypix, xrange=xr, $
                    yrange=yr, xseed=xseed, yseed=yseed, rseed=rseed, $
                    nside=nside, nximg=nximg, nyimg=nyimg, eps=eps, $
                    scheme=scheme, bin=bin, defarr=defarr, xsrc=xsrc, $
                    ysrc=ysrc, ixsrc=ixsrc, iysrc=iysrc, stars=stars, $
                    cells=cells

; -------- defaults
  if n_elements(nxpix) EQ 0 then nxpix = 300L
  if n_elements(nypix) EQ 0 then nypix = 300L
  if n_elements(xr) EQ 0    then xr = [-15,15] ; in stellar Rein
  if n_elements(yr) EQ 0    then yr = [-15,15] ; in stellar Rein
  if n_elements(nside) EQ 0 then nside = 16L ; nside of stellar population
  if n_elements(nximg) EQ 0 then nximg = 101L
  if n_elements(nyimg) EQ 0 then nyimg = 101L
  if n_elements(eps) EQ 0   then eps = 0.05 ; % of deflection in 6th moment



; -------- assure correct input type
  nximg = long(nximg)
  nyimg = long(nyimg)
  nxpix = long(nxpix)
  nypix = long(nypix)


; -------- utilities
  kappas   = kappa*fstar
  kappac   = kappa - kappas
  xmin     = xr[0]
  xmax     = xr[1]
  dx       = (xmax-xmin)/double(nxpix)
  ymin     = yr[0]
  ymax     = yr[1]
  dy       = (ymax-ymin)/double(nypix)
  xrimg    = 1.2*xr / abs((1.0 - kappa - gamma)) ; emperical factor
  yrimg    = 1.2*yr / abs((1.0 - kappa + gamma)) ; emperical factor
;  xrimg    = 1.2*xr / abs((1.0 - 0.3 - 0.3)) ; emperical factor
;  yrimg    = 1.2*yr / abs((1.0 - 0.3 + 0.3)) ; emperical factor
  xminimg  = xrimg[0]
  xmaximg  = xrimg[1]
  yminimg  = yrimg[0]
  ymaximg  = yrimg[1]
  dximg    = (xmaximg - xminimg)/double(nximg)
  dyimg    = (ymaximg - yminimg)/double(nyimg)
  xrst     = [xrimg[0] < yrimg[0], xrimg[1] > yrimg[1]]*2.2 ; empirical factor
  yrst     = xrst
  nraymag1 = dx*dy/(dximg*dyimg)

  splog, 'Shooting region size is ', dximg*nximg, ' by ', dyimg*nyimg, $
         form='(A,D6.2,A,D6.2)'

  splog, 'Number of rays per pixel for unity magnification is ', nraymag1, $
         form='(A,D6.2)'

  splog, 'size of stars region: ', xrst, yrst




; -------- initialize the stars and cells
  cells = ml_gencells(kappas, xrst, yrst, nside, xseed=xseed, yseed=yseed, $
                      rseed=rseed, stars=stars)

  splog, 'Mean number of stars per highest cell is ', $
         mean((cells.nstar)[where(cells.high)])



; -------- initialize the magnification and deflection angle maps
  defarr = dblarr(nximg, nyimg, 2)
  magarr = lonarr(nxpix, nypix)
  xpos   = xmin + dx*dindgen(nxpix)
  ypos   = ymin + dy*dindgen(nypix)



; -------- shoot the rays in the image plane (calculate deflection
;          angle at each point)
  common counters, cellcnt=cellcnt, starcnt=starcnt
  cellcnt = 0L
  starcnt = 0L


  ml_defarr, defarr, xrimg, yrimg, nximg, nyimg, cells, stars, scheme=scheme, $
             bin=bin

  splog, 'Total number of cells used:', cellcnt
  splog, 'Total number of stars used:', starcnt



; -------- gather rays into source plane pixels
  posgrid, ximg, yimg, xr=xrimg, yr=yrimg, nx=nximg, ny=nyimg


  xsrc  = ximg*(1.0-kappac-gamma) - defarr[*,*,0]
  ysrc  = yimg*(1.0-kappac+gamma) - defarr[*,*,1]
  ixsrc = reform(long(floor((xsrc-xmin)/dx)), nximg*nyimg)
  iysrc = reform(long(floor((ysrc-ymin)/dy)), nximg*nyimg)
  w     = where(ixsrc GE 0 and ixsrc LT nxpix and $
                iysrc GE 0 and iysrc LT nypix, nw)

  foo = where(ixsrc GE 0 and ixsrc LT nxpix, nfoo)
  bar = where(iysrc GE 0 and iysrc LT nypix, nbar)

  print, "numx land = ", nfoo
  print, "numy land = ", nbar
  print, "xmin, dx = ", xmin, dx
  print, "ymin, dy = ", ymin, dy
  print, "mean, std defx", mean(defarr[*,*,0]), stdev(defarr[*,*,0])
  print, "mean, std defy", mean(defarr[*,*,1]), stdev(defarr[*,*,1])

  if nw EQ 0 then message, 'NO RAYS DEFLECTED INTO SOURCE PLANE GRID!!!'

  isrc_ring = ixsrc[w] + iysrc[w]*nxpix
  magarr    = reform(magarr, long(nxpix)*long(nypix))

  magarr[isrc_ring]++

  magarr = reform(magarr, nxpix, nypix)

splog, 'Number of rays shot :   ', nximg*nyimg
splog, 'Number of rays landed : ', nw
splog, 'Double check : ', total(magarr)



; -------- convert number of rays into magnification
  magarr  /= nraymag1

  splog, 'Number of rays per unity magnification:', nraymag1, form='(A,D10.3)'

  return, magarr
end

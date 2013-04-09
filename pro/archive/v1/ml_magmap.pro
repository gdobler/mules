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
pro ml_magmap, kappa, gamma, fstar, nxpix=nxpix, nypix=nypix, xrange=xr, $
               yrange=yr, xseed=xseed, yseed=yseed, rseed=rseed, nside=nside, $
               nximg=nximg, nyimg=nyimg, eps=eps

; -------- defaults
  if n_elements(nxpix) EQ 0 then nxpix = 1000L
  if n_elements(nypix) EQ 0 then nypix = 1000L
  if n_elements(xr) EQ 0    then xr = [-100,100] ; in stellar Rein
  if n_elements(yr) EQ 0    then yr = [-100,100] ; in stellar Rein
  if n_elements(nside) EQ 0 then nside = 8L ; nside of stellar population
  if n_elements(nximg) EQ 0 then nximg = 100L
  if n_elements(nyimg) EQ 0 then nyimg = 100L
  if n_elements(eps) EQ 0   then eps = 0.05 ; % of deflection in 6th moment



; -------- utilities
  kappas = kappa*fstar
  kappac = kappa - kappas
  xmin   = xr[0]
  xmax   = xr[1]
  dx     = (xmax-xmin)/double(nxpix)
  ymin   = yr[0]
  ymax   = yr[1]
  dy     = (ymax-ymin)/double(nypix)
  xrimg  = xr * (1.0 - kappa - gamma)
  yrimg  = yr * (1.0 - kappa + gamma)
  xminimg = xrimg[0]
  xmaximg = xrimg[1]
  yminimg = xrimg[0]
  ymaximg = xrimg[1]
  dximg   = (xmaximg - xminimg)/double(nximg)
  dyimg   = (ymaximg - yminimg)/double(nyimg)



; -------- initialize the stars and cells
  cells = ml_gencells(kappas, fstar, xrange, yrange, nside, xseed=xseed, $
                      yseed=yseed, rseed=rseed, stars=stars)




; -------- initialize the magnification and deflection angle maps
  defarr = dblarr(nximg, nyimg)
  magarr = lonarr(nxpix, nypix)
  xpos   = xmin + dx*dindgen(nxpix)
  ypos   = ymin + dy*dindgen(nypix)



; -------- shoot the rays in the image plane (calculate deflection
;          angle at each point)
  ml_defarr, deffarr, xrimg, yrimg, nximg, nyimg



; -------- gather rays into source plane pixels

; -------- convert number of rays into magnification


  return
end

;+
; NAME:
;  ml_genstars
;
; PURPOSE:
;  Generate a stellar field given an input number of stars and a
;  window in the lens plane.  Output stellar field as the stars
;  positions, masses (from mass function [GGD: NOT WORKING YET!!!]),
;  and (highest) tree cell number.
;
; CALLING SEQUENCE:
;  stars = ml_genstars(nstars, xrange, yrange, nside [, xseed=,
;                      yseed=, rseed= ])
;
; INPUTS:
;  nstars - total number of stars
;  xrange - xrange of lens plane where stars are included
;  yrange - yrange of lens plane where stars are included
;  nside  - number of cells along one side of tree
;
; OPTIONAL INPUTS:
;  xseed  - seed for x positions (def: 314L)
;  yseed  - seed for y positions (def: 271L)
;  rseed  - seed for einstein radii GGD: NOT WORKING YET!!!
;
; KEYWORDS:
;
; OUTPUTS:
;  stars - structure containing stars position, rein, and highest cell number
;
; OPTIONAL OUTPUTS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;  2011/07/15 - Written by Greg Dobler (KITP/UCSB)
;
;------------------------------------------------------------
function ml_genstars, nstars, xrange, yrange, nside, xseed=xseed, $
                      yseed=yseed, rseed=rseed

; -------- defaults
  if n_elements(xseed) EQ 0 then xseed = 314L
  if n_elements(yseed) EQ 0 then yseed = 271L



; -------- generate the positions
  xstar = xrange[0] + (xrange[1]-xrange[0])*randomu(xseed, nstars)
  ystar = yrange[0] + (yrange[1]-yrange[0])*randomu(yseed, nstars)
  rein = dblarr(nstars) + 1.0 ; GGD: use mass function sooner or later



; -------- put stars into ring-ordered cells
  dx = (xrange[1]-xrange[0])/double(nside)
  dy = (yrange[1]-yrange[0])/double(nside)

  xcell = floor((xstar-xrange[0])/dx)
  ycell = floor((ystar-yrange[0])/dy)

  icell_ring = xcell + nside*ycell



; -------- now convert to nested cells (for the tree algorithm)
  nnum  = round(alog(double(nside))/alog(2.0)) - 1L
  offr  = nnum GT 0 ? total(2.0^(2.0*dindgen(nnum)+1), /cum) : 0
  offc  = nnum GT 0 ? -total(2.0^(2.0*(dindgen(nnum)+1)), /cum)-2.0-nside : 0

  row    = lonarr(nside)
  col    = lonarr(nside) - 2L - nside
  col[0] = 0 ; special case

  for ioff=0L, nnum-1 do begin
     tind = long(2L^(ioff+2)*lindgen(nside/(2^(ioff+2))) + 2L^(ioff+1))

     row[tind] = offr[ioff]
     col[tind] = offc[ioff]
  endfor

  arow = long(total(row, /cum))
  acol = long(total(col, /cum))

  icell = icell_ring + acol[ycell] + arow[xcell] + nside^2 - 2



; -------- put everything into a structure ordering by cell number
  ord = sort(icell)

  stars = {xstar : double(xstar[ord]), $
           ystar : double(ystar[ord]), $
           rein  : double(rein[ord]), $
           icell : long(icell[ord])}

  return, stars
end

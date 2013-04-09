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
;  2011/07/18 - Written by Greg Dobler (KITP/UCSB)
;
;------------------------------------------------------------
function ml_cellpos, xrange, yrange, nside

; -------- utilities
  nlev = round(alog(double(nside))/alog(2.0))



; -------- loop through the levels
  for ilev=0L, nlev-1 do begin
     tside = 2L^(ilev+1)
     xmin  = xrange[0]
     xmax  = xrange[1]
     dx    = (xmax-xmin)/double(tside)
     ymin  = yrange[0]
     ymax  = yrange[1]
     dy    = (ymax-ymin)/double(tside)

     ; get ring ordered cell numbers and x and y positions
     icell_ring = lindgen(tside^2)
     xring      = icell_ring mod tside
     yring      = icell_ring / tside
     xpos       = xmin + dx*(xring+0.5)
     ypos       = ymin + dy*(yring+0.5)

     ; convert ring cell numbers to nested cell numbers
     nnum  = round(alog(double(tside))/alog(2.0)) - 1L
     offr  = nnum GT 0 ? total(2.0^(2.0*dindgen(nnum)+1), /cum) : 0
     offc  = nnum GT 0 ? $
             -total(2.0^(2.0*(dindgen(nnum)+1)), /cum)-2.0-tside : 0
     row    = lonarr(tside)
     col    = lonarr(tside) - 2L - tside
     col[0] = 0 ; special case

     for ioff=0L, nnum-1 do begin
        tind = long(2L^(ioff+2)*lindgen(tside/(2^(ioff+2))) + 2L^(ioff+1))

        row[tind] = offr[ioff]
        col[tind] = offc[ioff]
     endfor

     arow = long(total(row, /cum))
     acol = long(total(col, /cum))

     icell = lonarr(tside^2)

     for i=0L, tside-1 do $
        icell[i*tside:(i+1)*tside-1] = icell_ring[i*tside:(i+1)*tside-1] + $
        acol[i] + arow + tside^2 - 2

     ; stack the cell numbers and positions together
     sicell = ilev EQ 0 ? icell : [sicell, icell]
     sxpos  = ilev EQ 0 ? xpos : [sxpos, xpos]
     sypos  = ilev EQ 0 ? ypos : [sypos, ypos]
  endfor



; -------- order by cell number, concatenate, and return
  ord  = sort(sicell)
  cpos = [[sxpos[ord]],[sypos[ord]]]

  return, cpos
end

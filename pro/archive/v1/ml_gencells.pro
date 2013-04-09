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
function ml_gencells, kappas, xrange, yrange, nside, xseed=xseed, $
                      yseed=yseed, rseed=rseed, stars=stars

; -------- utilities
  xmin  = xrange[0]
  xmax  = xrange[1]
  dx    = xmax-xmin
  ymin  = yrange[0]
  ymax  = yrange[1]
  dy    = ymax-ymin
  cpos  = ml_cellpos(xrange, yrange, nside)
  nstar = kappas ;GGD FIX!!!
  nlev  = long(alog(nside)/alog(2))



; -------- generate the stars and find the maximum number of stars in
;          a cell
  stars   = ml_genstars(nstar, xrange, yrange, nside, xseed=xseed, $
                        yseed=yseed, rseed=rseed)
  cmin    = min(stars.icell)
  cmax    = max(stars.icell)
  ncell   = nside^2
  histst  = histogram(stars.icell, bin=1)
  maxcell = min(stars.icell) + (where(histst EQ max(histst)))[0]
  dum     = where(stars.icell EQ maxcell, maxst) ; max # of stars in a cell



; -------- loop through the cells and calculate the multipole moments
  for ilev=0L, nlev-1 do begin
     splog, 'Building level ', ilev, ' of ', nlev-1, '...', $
            form='(A,I2,A,I2,A)'

     tncell = 4L^(ilev+1)
     tcmin  = long(total(4.0^lindgen(ilev+1))-1)
     tcmax  = tcmin + tncell
     dcell  = ncell/tncell

     for icell=0L, tncell-1 do begin
        tcell = tcmin + icell

        stind = where(stars.icell GE cmin + icell*dcell and $
                      stars.icell LT cmin + (icell+1)*dcell, tnstar)

        if tnstar GT 0 then begin
           tposx  = cpos[tcell,0] ; mean x position of cell
           tposy  = cpos[tcell,1] ; mean y position of cell
           tmass  = total(stars.rein[stind]^2) ; total mass (Rein^2) of cell

           delxi = stars.xstar[stind] - tposx
           delyi = stars.ystar[stind] - tposy
           ri    = sqrt(delxi^2+delyi^2)
           ti    = atan(delyi, delxi)
           mi    = stars.rein[stind]

           c1 = total(mi^2 * ri^1 * cos(1.0*ti))
           c2 = total(mi^2 * ri^2 * cos(2.0*ti))
           c3 = total(mi^2 * ri^3 * cos(3.0*ti))
           c4 = total(mi^2 * ri^4 * cos(4.0*ti))
           c5 = total(mi^2 * ri^5 * cos(5.0*ti))
           c6 = total(mi^2 * ri^6 * cos(6.0*ti))

           s1 = total(mi^2 * ri^1 * sin(1.0*ti))
           s2 = total(mi^2 * ri^2 * sin(2.0*ti))
           s3 = total(mi^2 * ri^3 * sin(3.0*ti))
           s4 = total(mi^2 * ri^4 * sin(4.0*ti))
           s5 = total(mi^2 * ri^5 * sin(5.0*ti))
           s6 = total(mi^2 * ri^6 * sin(6.0*ti))

           thigh  = tcell GE cmin ; flag highest cell (don't try to break up)
           tstind = lonarr(maxst) - 1
           if thigh then $
              tstind[0:tnstar-1] = stind ; indices of all stars in this cell
        endif else begin
           tposx  = cpos[tcell,0] ; mean x position of cell
           tposy  = cpos[tcell,1] ; mean y position of cell
           tmass  = 0.0d
           c1     = 0.0d
           c2     = 0.0d
           c3     = 0.0d
           c4     = 0.0d
           c5     = 0.0d
           c6     = 0.0d
           s1     = 0.0d
           s2     = 0.0d
           s3     = 0.0d
           s4     = 0.0d
           s5     = 0.0d
           s6     = 0.0d
           thigh  = tcell GE cmin ; flag highest cell (don't try to break up)
           tstind = lonarr(maxst) - 1
        endelse

        tstr = {icell : tcell, $
                nstar : tnstar, $
                xcell : tposx, $
                ycell : tposy, $
                mcell : tmass, $
                cmom  : [c1,c2,c3,c4,c5,c6], $
                smom  : [s1,s2,s3,s4,s5,s6], $
                high  : thigh, $
                stind : tstind }

        cells = n_elements(cells) EQ 0 ? tstr : [cells, tstr]
     endfor
  endfor

  return, cells
end

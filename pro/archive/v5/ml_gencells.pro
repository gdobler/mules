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
  nstar = kappas*dx*dy/!dpi ; Nstar = kappa*Lx*Ly/pi*bbar^2  GGD: fix bbar!!!
  cpos  = ml_cellpos(xrange, yrange, nside)
  nlev  = long(alog(nside)/alog(2))

  splog, 'Creating cell structure with ', nstar, ' stars...', form='(A,I10,A)'



; -------- generate the stars and find the maximum number of stars in
;          a cell
  stars   = ml_genstars(nstar, xrange, yrange, nside, xseed=xseed, $
                        yseed=yseed, rseed=rseed)
  cmin    = long(total(4.0^lindgen(nlev))-1)
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
           tposx = cpos[tcell,0] ; mean x position of cell
           tposy = cpos[tcell,1] ; mean y position of cell
           tmass = total(stars.rein[stind]^2) ; total mass (Rein^2) of cell

           delxi = tposx - stars.xstar[stind]
           delyi = tposy - stars.ystar[stind]
           ri    = sqrt(delxi^2+delyi^2)
           ti    = atan(delyi, delxi)
           mi    = stars.rein[stind]

           mm = dindgen(20) + 1.0d

           for imm=0L, n_elements(mm)-1 do begin
              tcmom = imm EQ 0 ? total(mi^2 * ri^mm[imm] * cos(mm[imm]*ti)) : $
                      [tcmom, total(mi^2 * ri^mm[imm] * cos(mm[imm]*ti))]

              tsmom = imm EQ 0 ? total(mi^2 * ri^mm[imm] * sin(mm[imm]*ti)) : $
                      [tsmom, total(mi^2 * ri^mm[imm] * sin(mm[imm]*ti))]
           endfor

           thigh  = tcell GE cmin ; flag highest cell (don't try to break up)
           tstind = lonarr(maxst) - 1
           if thigh then $
              tstind[0:tnstar-1] = stind ; indices of all stars in this cell
        endif else begin
           tposx  = cpos[tcell,0] ; mean x position of cell
           tposy  = cpos[tcell,1] ; mean y position of cell
           tmass  = 0.0d
           tcmom  = dblarr(20)
           tsmom  = dblarr(20)
           thigh  = tcell GE cmin ; flag highest cell (don't try to break up)
           tstind = lonarr(maxst) - 1
        endelse

        tstr = {icell : tcell, $
                nstar : tnstar, $
                xcell : tposx, $
                ycell : tposy, $
                mcell : tmass, $
                cmom  : tcmom, $
                smom  : tsmom, $
                high  : thigh, $
                stind : tstind }

        cells = n_elements(cells) EQ 0 ? tstr : [cells, tstr]
     endfor
  endfor

  return, cells
end

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
;  2011/07/25 - Written by Greg Dobler (KITP/UCSB)
;
;------------------------------------------------------------
function ml_defang, ximg, yimg, cells, stars, clist

  common counters


; -------- utilities
  ncellnum = n_elements(clist)
  defang   = dblarr(2)



; -------- Loop through input cell list
  for icellnum=0L, ncellnum-1 do begin
     tcellnum = clist[icellnum]

     if tcellnum LT 0 then begin
        if NOT cells[tcellnum].high then message, 'Cell list level ERROR!!!'

        stind = cells[tcellnum].stind
        stind = stind[where(stind GE 0)]

        starcnt += n_elements(stind)

        mi      = stars.rein[stind]
        delxi   = ximg - stars.xstar[stind]
        delyi   = yimg - stars.ystar[stind]
        ri      = sqrt(delxi^2+delyi^2)
        alphax  = total(mi*mi*delxi/ri^2)
        alphay  = total(mi*mi*delyi/ri^2)
        defang += [alphax, alphay]
     endif else begin
        cellcnt++

        delx   = ximg - cells[tcellnum].xcell
        dely   = yimg - cells[tcellnum].ycell
        rcell  = sqrt(delx^2 + dely^2)
        tcell  = atan(dely,delx)
        mcell  = cells[tcellnum].mcell
        cmom   = cells[tcellnum].cmom
        smom   = cells[tcellnum].smom
        nstar  = cells[tcellnum].nstar

        if nstar EQ 0 then message, 'Cell list using empty cell!!!'

        nmm      = n_elements(cmom)
        mm       = dindgen(nmm)+1.0
        oormp1   = 1.0d/rcell^(mm+1.0d)
        alphar0  = mcell/rcell
        alphar   = (cmom*cos(mm*tcell)+smom*sin(mm*tcell))*oormp1
        alphat   = (cmom*sin(mm*tcell)-smom*cos(mm*tcell))*oormp1
        talphar  = alphar0 + total(alphar)
        talphat  = total(alphat)
        alphax   = cos(tcell)*talphar - sin(tcell)*talphat
        alphay   = sin(tcell)*talphar + cos(tcell)*talphat
        defang  += [alphax, alphay]
     endelse
  endfor

  return, defang
end

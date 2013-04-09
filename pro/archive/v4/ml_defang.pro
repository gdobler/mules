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
function ml_defang, ximg, yimg, cells, stars, cellnum, clist=clist

  common counters

; -------- initialize the deflection angle and run lowest level cell
  if n_elements(cellnum) EQ 0 then begin
     defang0 = ml_defang(ximg, yimg, cells, stars, 0, clist=clist0)
     defang1 = ml_defang(ximg, yimg, cells, stars, 1, clist=clist1)
     defang2 = ml_defang(ximg, yimg, cells, stars, 2, clist=clist2)
     defang3 = ml_defang(ximg, yimg, cells, stars, 3, clist=clist3)

     clist = [clist0, clist1, clist2, clist3]
     clist = clist[sort(abs(clist))]

     return, defang0 + defang1 + defang2 + defang3
  endif



; -------- utilities
  delx   = ximg - cells[cellnum].xcell
  dely   = yimg - cells[cellnum].ycell
  rcell  = sqrt(delx^2 + dely^2)
  tcell  = atan(dely,delx)
  mcell  = cells[cellnum].mcell
  cmom   = cells[cellnum].cmom
  smom   = cells[cellnum].smom
  nstar  = cells[cellnum].nstar
  defang = dblarr(2)

  if nstar EQ 0 then return, defang



; -------- calculate the deflection angle up to order 20
  nmm     = n_elements(cmom)
  mm      = dindgen(nmm)+1.0
  oormp1  = 1.0d/rcell^(mm+1.0d)
  alphar0 = mcell/rcell
  alphar  = (cmom*cos(mm*tcell)+smom*sin(mm*tcell))*oormp1
  alphat  = (cmom*sin(mm*tcell)-smom*cos(mm*tcell))*oormp1

  if total((sqrt(alphar^2+alphat^2))[nmm-11:nmm-1]) GE 0.1 then begin
     if cells[cellnum].high then begin
        stind = cells[cellnum].stind
        stind = stind[where(stind GE 0)]

        starcnt += n_elements(stind)

        mi     = stars.rein[stind]
        delxi  = ximg - stars.xstar[stind]
        delyi  = yimg - stars.ystar[stind]
        ri     = sqrt(delxi^2+delyi^2)
        alphax = total(mi*mi*delxi/ri^2)
        alphay = total(mi*mi*delyi/ri^2)
        defang = [alphax, alphay]

        clist = n_elements(clist) EQ 0 ? -cellnum : [clist,-cellnum]
     endif else for icell=0L, 3 do $
        defang += ml_defang(ximg, yimg, cells, stars, 4L*(cellnum+1L)+icell)
  endif else begin
     cellcnt++

     talphar = alphar0 + total(alphar)
     talphat = total(alphat)
     alphax  = cos(tcell)*talphar - sin(tcell)*talphat
     alphay  = sin(tcell)*talphar + cos(tcell)*talphat
     defang  = [alphax, alphay]

     clist = n_elements(clist) EQ 0 ? cellnum : [clist,cellnum]
  endelse

  return, defang
end

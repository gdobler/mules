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
function ml_defang, ximg, yimg, cells, stars, cellnum

; -------- initialize the deflection angle and run lowest level cell
  if n_elements(cellnum) EQ 0 then begin
     defang0 = ml_defang(ximg, yimg, cells, stars, 0)
     defang1 = ml_defang(ximg, yimg, cells, stars, 1)
     defang2 = ml_defang(ximg, yimg, cells, stars, 2)
     defang3 = ml_defang(ximg, yimg, cells, stars, 3)

     return, defang0 + defang1 + defang2 + defang3
  endif



; -------- utilities
  delx  = cells[cellnum].xcell - ximg
  dely  = cells[cellnum].ycell - yimg
  rcell = sqrt(delx^2 + dely^2) > 1d-6 ; GGD: careful here, this is a scale!!!
  tcell = atan(dely,delx)
  mcell = cells[cellnum].mcell
  cmom  = cells[cellnum].cmom
  smom  = cells[cellnum].smom



; -------- calculate the deflection angle up to order 6
  mm      = dindgen(6)+1.0
  oormp1  = 1.0d/rcell^(mm+1.0d)
  alphar0 = mcell*alog(rcell)
  alphar  = (cmom*cos(mm*tcell)+smom*sin(mm*tcell))*oormp1
  alphat  = (cmom*sin(mm*tcell)-smom*cos(mm*tcell))*oormp1

  print, cellnum, form='(I5)'
;  print, total(alphar), 100.*alphar/total(alphar), form='(7D8.2)'
;  print, total(alphat), 100.*alphat/total(alphat), form='(7D8.2)'
  print, alphar, form='(6D10.2)'
  print, alphat, form='(6D10.2)'


; Cell number i contains cells numbered 4*(cell+1) to 4*(cell+1)+3
; For this cell evaluate the deflection angle including up to order 6
; If the 6th moment represents 10% of the total deflection angle, then
; do 1 of 2 things:
;   if highest keyword is set, use all the stars else recall the
;   function with cells numbered 4*(icell+1) to 4*(icell+1)+3.

;  print, cellnum, delx, dely, cells[cellnum].xcell, cells[cellnum].ycell, $
;         form='(I4,4D6.2)'

  if cells[cellnum].high then return, 1.0 else defang = 0.0

  for icell=0L, 3 do $
     defang += ml_defang(ximg, yimg, cells, stars, 4L*(cellnum+1L)+icell)

  return, defang
end

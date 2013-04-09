; This will be the deflection array and needs to be recursive.  First
; we get the deflection angle for the ul, ur, ll, lr.  Each one of
; these evaluations should also include a list of cells used in the
; calculation of the deflection angle.  Recursively go through the
; entire thing.  If four corners have the exact same cell list, pass
; that list to all of the following subgrids.

pro ml_defarr_recur, defarr, i1, i2, j1, j2, mode, xr, yr, nx, ny

; -------- evaluate the necessary corners of the grid
  if mode EQ 0 then begin
     v1 = ml_defang(i1,j1)
     v2 = ml_defang(i2,j1)
     v3 = ml_defang(i1,j2)
     v4 = ml_defang(i2,j2)
  endif else if mode EQ 1 then begin
     v2 = ml_defang(i2,j1)
     v3 = ml_defang(i1,j2)
     v4 = ml_defang(i2,j2)
  endif else if mode EQ 2 or mode EQ 3 then begin
     v4 = ml_defang(i2,j2)
  endif
  
  return
end


pro ml_deffarr, defarr, xr, yr, nx, ny

; -------- GGD: ADD RECURSIVE LOOP FOR SPEED!!!

; -------- loop directly
  for iximg=0L, nx-1 do begin
     for iyimg=0L, ny-1 do begin

        ximg = xr[0] + (xr[1]-xr[0])*double(iximg)/double(nx)
        yimg = yr[0] + (yr[1]-yr[0])*double(iyimg)/double(ny)

        defarr[iximg,iyimg] = ml_defang(ximg, yimg, cells, stars)

     endfor
  endfor






  return
end

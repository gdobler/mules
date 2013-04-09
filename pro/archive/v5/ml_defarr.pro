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


pro ml_defarr, defarr, xr, yr, nx, ny, cells, stars, scheme=scheme, bin=bin

; -------- defaults
  if n_elements(scheme) EQ 0 then scheme=1


; -------- loop directly
  if scheme EQ 1 then begin
     for iximg=0L, nx-1 do begin
        splog, 'iximg = ', iximg, form='(A,I4)'

        for iyimg=0L, ny-1 do begin
           ximg = xr[0] + (xr[1]-xr[0])*double(iximg)/double(nx)
           yimg = yr[0] + (yr[1]-yr[0])*double(iyimg)/double(ny)

           defarr[iximg,iyimg,*] = ml_defang(ximg, yimg, cells, stars)
        endfor
     endfor

     return
  endif



; -------- linearly interpolate in fixed cells
  if scheme EQ 2 then begin
     if n_elements(bin) EQ 0 then bin = 64L

     ncx = (nx-1)/bin
     ncy = (ny-1)/bin

     for icx=0L, ncx-1 do begin
        splog, 'icx = ', icx, ' out of ', ncx, form='(A,I5,A,I5)'

        for icy=0L, ncy-1 do begin

           ; get four corners

           iximgul = icx*bin
           iyimgul = (icy+1)*bin
           iximgur = (icx+1)*bin
           iyimgur = (icy+1)*bin
           iximgll = icx*bin
           iyimgll = icy*bin
           iximglr = (icx+1)*bin
           iyimglr = icy*bin

           ximgul = xr[0] + (xr[1]-xr[0])*double(iximgul)/double(nx)
           yimgul = yr[0] + (yr[1]-yr[0])*double(iyimgul)/double(ny)
           ximgur = xr[0] + (xr[1]-xr[0])*double(iximgur)/double(nx)
           yimgur = yr[0] + (yr[1]-yr[0])*double(iyimgur)/double(ny)
           ximgll = xr[0] + (xr[1]-xr[0])*double(iximgll)/double(nx)
           yimgll = yr[0] + (yr[1]-yr[0])*double(iyimgll)/double(ny)
           ximglr = xr[0] + (xr[1]-xr[0])*double(iximglr)/double(nx)
           yimglr = yr[0] + (yr[1]-yr[0])*double(iyimglr)/double(ny)

           ; calculate deflection angles, only need to do all the 
           ; corners for special cases
           defur = ml_defang(ximgur, yimgur, cells, stars)

           if icx EQ 0 then $
              deful = ml_defang(ximgul, yimgul, cells, stars) $
           else $
              deful = defarr[iximgul,iyimgul,*]

           if icy EQ 0 then $
              deflr = ml_defang(ximglr, yimglr, cells, stars) $
           else $
              deflr = defarr[iximglr,iyimglr,*]

           if icx EQ 0 and icy EQ 0 then $
              defll = ml_defang(ximgll, yimgll, cells, stars) $
           else $
              defll = defarr[iximgll,iyimgll,*]


           ; funny stacking... flipped from what I'd expect
           defx  = [[defll[0],deflr[0]], [deful[0],defur[0]]]
           defy  = [[defll[1],deflr[1]], [deful[1],defur[1]]]

           defxint = interpolate(defx, dindgen(bin+1)/double(bin), $
                                 dindgen(bin+1)/double(bin), /grid)
           defyint = interpolate(defy, dindgen(bin+1)/double(bin), $
                                 dindgen(bin+1)/double(bin), /grid)

           defarr[iximgll:iximgur,iyimgll:iyimgur,0] = defxint
           defarr[iximgll:iximgur,iyimgll:iyimgur,1] = defyint
        endfor
     endfor

     return
  endif




  return
end

pro runit, defarr

  cells = ml_gencells(0.45, [-1000.,1000.], [-1000.,1000.], 32L, stars=stars)

  defarr = dblarr(100,100,2)

  xr = [-500.,500]
  yr = [-500.,500]
  nx = 100L
  ny = 100L

  ml_defarr, defarr, xr, yr, nx, ny, cells, stars

  return
end

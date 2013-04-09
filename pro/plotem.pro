pro plotem, cells

  i = 0
  plot, [cells[i].xcell], [cells[i].ycell], xr=[-15,15], yr=[-15,15], $
        symsize=2, psym=1

  wait, 0.5

  for i=1, n_elements(cells)-1 do begin
     oplot, [cells[i].xcell], [cells[i].ycell], symsize=2, psym=1
     wait, 0.5
  endfor


  return
end

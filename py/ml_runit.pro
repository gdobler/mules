;+
;
;  2011/10/07 - Test case for PMN 1632 (see Dobler, Keeton, &
;               Wambsganss, 2007)
;
pro ml_runit

  kapA = 0.370
  kapB = 3.309
  kapC = 16.517

  gamA = 0.422
  gamB = 2.859
  gamC = 13.796

  mapA100 = ml_magmap(kapA, gamA, 1.0, xr=[-50,50], yr=[-50,50], $
                      scheme=2, nximg=8001L, nyimg=8001L, nside=32l)

  mapB100 = ml_magmap(kapB, gamB, 1.0, xr=[-50,50], yr=[-50,50], $
                      scheme=2, nximg=8001L, nyimg=8001L, nside=32l)

  mapC100 = ml_magmap(kapC, gamC, 1.0, xr=[-50,50], yr=[-50,50], $
                      scheme=2, nximg=8001L, nyimg=8001L, nside=32l)



  return
end

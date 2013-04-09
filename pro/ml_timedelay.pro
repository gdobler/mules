pro ml_timedelay

  
; -------- make a light curve and shift it by 118.5 days then resample
;          at 1 day cadence

  mag0 = 20.d
  td   = 118.5d

  lc = lc_car_gen(314l, meanmag=mag0, year=100L, time=time, tau=10.d^3.5)
  lcA = lc
  lcB = shift(lc, long(td*10.))

  kapA = 0.370
  gamA = 0.422
  kapB = 3.309
  gamB = 2.859

  muA = 1.0/((1.0-kapA)^2-gamA^2)
  muB = 1.0/((1.0-kapB)^2-gamB^2)

  magA = 18.93 ; from castles
  magB = 20.94 ; from castles

  lcA  += magA-mag0
  lcB  += MagB-mag0

  ; rebin and excise the middle 2 years
  index = 10L*lindgen(round(n_elements(lcA)/10.))

  lcA_bin = lcA[index]
  lcB_bin = lcB[index]
  lcA_cut = lcA_bin[500:500+730]
  lcB_cut = lcB_bin[500:500+730]

  days = indgen(n_elements(lcA_cut))


; -------- read in the magmaps and apply fluctuations
  magmapA  = mrdfits('~/microlens/src/1632A.fits', 0, hA)
  rayamp1A = sxpar(hA,'RAYAMP1')
  magmapA /= rayamp1A*muA

  magmapB  = mrdfits('~/microlens/src/1632B.fits', 0, hB)
  rayamp1B = sxpar(hB,'RAYAMP1')
  magmapB /= rayamp1B*muB

  magmapB  = -2.5*alog10(abs(magmapB))

  ; start in the middle (forget A since there are barely any fluctuations)
  midiB = round(n_elements(magmapB[*,0])/2.0)


  ; each pixel is 0.03 Rein so half an Rein is 17 pixels

  mlmagB = congrid(magmapB[midi:midi+17], n_elements(days), cubic=-0.5)

  lcB_cut_ml = lcB_cut + mlmagB

  





  return
end

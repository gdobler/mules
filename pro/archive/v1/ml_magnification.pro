;+
; NAME:
;  ml_magnification
;
; PURPOSE:
;  Calculate the magnification for a point in the source plane for a
;  given realization of stars.
;
; CALLING SEQUENCE:
;  mag = ml_magnification(x, y, kappa, gamma, stars, cells)
;
; INPUTS:
;  x     - source plane x position in units of stellar einstein radii
;  y     - source plane y position in units of stellar einstein radii
;  kappa - local smooth convergence
;  gamma - local shear (aligned with x-axis)
;  stars - structure containing stars (see ml_genstars)
;  cells - structure containing cells (see ml_gencells)
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
;  2011/07/20 - Written by Greg Dobler (KITP/UCSB)
;
;------------------------------------------------------------
function ml_magnification, x, y, kappa, gamma, stars, cells



; -------- calculate the magnification tensor and its inverse det
  muinv      = dblarr(2,2)
  muinv[0,0] = 1.0 - phixx - kappa - gamma
  muinv[1,0] = -phixy
  muinv[0,1] = -phixy
  muinv[1,1] = 1.0 - phiyy - kappa + gamma

  mag = 

  return, mag
end

function jrandomn, seed, d1, d2, d3, d4, d5, d6, d7, d8, use_seed = use, _extra = extra
  
;Fixes a bug in randomn:  In that routine, if seed is not specified, 
;a new one is generated.  Unfortunately, this new one is not generated randomly, but based on an internal sequence.  This means that if you rely on this new seed to be random (as most users do), that sequential calls to your routine that calls randomn will have very similar outputs:
;For instance, try the following code:
;
; delvar, seed, seed2
; print, randomn(seed,2)
; print, randomn(seed2,2)
; print, randomn(seed,2)
;
; the undefined variable "seed2" has not been given a random seed -- it has the next value of "seed"

;
; jrandom[u|n] addresses this by ensuring that:
;   1) sequential calls with undefined seeds to jrandom ALWAYS use a new seed by keeping track of the current seed with a global variable "!random_seed"
;   2) ignoring input values of the seed (which may be "stale" values already used by some other call to random) unless explicitly told otherwise.
;
; Usage:  IDENTICAL to randomn, except:
  
;IF YOU WISH TO SPECIFY THE SEED EXPLICITLY, YOU MUST MUST MUST SET USE_SEED=1

;Jason Wright 
;19 Jan 2009 (less than 24 hours now...)
  
  defsysv, '!random_seed', exists = x   ;does the global exist?

  if keyword_set(use) and n_elements(seed) ne 0 then sx = 1 else sx = 0;default to get new seeds

  if ~sx then begin  ;no existing seed -- get new one
    if x then seed = !random_seed else begin ;if present, use the global value
      r = randomn(blah, 1)               ;otherwise get a new seed
      seed = blah
    endelse
  endif 

  if ~x then begin                      ;if not
    defsysv, '!random_seed', seed    ;make a global
  endif

  case n_params()-1 of                  ;call randomn.pro
    0:ans = randomn(seed, _extra = extra)
    1:ans = randomn(seed, d1, _extra = extra)
    2:ans = randomn(seed, d1, d2, _extra = extra)
    3:ans = randomn(seed, d1, d2, d3, _extra = extra)
    4:ans = randomn(seed, d1, d2, d3, d4, _extra = extra)
    5:ans = randomn(seed, d1, d2, d3, d4, d5, _extra = extra)
    6:ans = randomn(seed, d1, d2, d3, d4, d5, d6, _extra = extra)
    7:ans = randomn(seed, d1, d2, d3, d4, d5, d6, d7, _extra = extra)
    else: ans = randomn(seed, d1, d2, d3, d4, d5, d6, d7, d8, _extra = extra)
  endcase
  !random_seed = seed                   ;update the global
  return, ans
end
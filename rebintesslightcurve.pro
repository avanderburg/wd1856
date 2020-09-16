function rebintesslightcurve, ts, fs, qs, es = es, num = num, cadence = cadence, ngood = ngood, method1 = method1
  if 1 - keyword_set(cadence) then cadence = 2d ; usually working with 2 minute cadence data.
  if 1 - keyword_set(num) then num = 15d; change TESS 2 minute cadence to 30 minute cadence
  if n_elements(qs) eq 0 then qs = dblarr(n_elements(ts))
  ; This is not going to work quite right where there is missing data, but it's good enough to do precision tests with
  method2 = 1-keyword_set(method1)
  if 1 - keyword_set(ngood) then ngood = floor(num*2d/3d) ; usually want 10 observations out of 15 in order to bin 2 minute tess data into a 30 minute cadence light curve

  if 1 - method2 then begin
    t = dblarr(floor(n_elements(ts)/num))
    f = t
    q = t

    for  j = 0, n_elements(t) -1 do begin
      ;thisbin = [i*15: (i+1)*15-1]
      if  (j+1)*num-1 lt n_elements(ts) then begin
        t[j] = mean(ts[j*num: (j+1)*num-1])
        f[j] = total(fs[j*num: (j+1)*num-1])/num
        q[j] = total(qs[j*num: (j+1)*num-1])/num
      end
;      if  (j+1)*15-1 ge n_elements(ts) then begin
;        thisnum = n_elements(ts) - j*15
;        t[j] = mean(ts[j*num:*])
;        f[j] = total(fs[j*num:*])/thisnum
;        q[j] = total(qs[j*num:*])/thisnum
;      end
    endfor

    out = {t:t, f:f, q:q}
  end


  bint = cadence * num
  if method2 then begin
    ut = fillarr(bint/24d/60d, min(ts), max(ts))
    nut = intarr(n_elements(ut))
    fut = dblarr(n_elements(ut))
    qut = dblarr(n_elements(ut))
    if keyword_set(es) then  eut = dblarr(n_elements(ut))
    


    for i = 0, n_elements(ut)-1 do begin
      thisone = where(abs(ut[i] - ts) le bint/24d/60d/2)
      if thisone[0] gt -1 then begin
        ;fut[i] = mean(rf[thisone]/rs[thisone] );- model[thisone] + 1
        fut[i] = mean(fs[thisone])
        nut[i] = n_elements(thisone)
        qut[i] = total(qs(thisone))
        if keyword_set(es) then eut[i] = sqrt(1 / total(1 / es[thisone]^2))
      endif
    endfor

    good = where(nut gt ngood)
    ut = ut(good)
    nut = nut(good)
    fut = fut(good)
    qut = qut(good)
    if keyword_set(es) then eut = eut(good)

    out = {t:ut, f:fut, q:qut}
    if keyword_set(es) then out = {t:ut, f:fut, q:qut, e:eut}


  endif



  return, out

end
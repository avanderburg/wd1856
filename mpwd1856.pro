function mpwd1856, p, tspit=tspit,fspit=fspit, espit=espit,tgtc=tgtc,fgtc=fgtc, egtc=egtc, modelspit = modelspit, modelgtc = modelgtc, $
  mp = mp, b1 = b1, cosi = cosi, lowfrgtc=lowfrgtc, lowfrspit = lowfrspit, dolowfr = dolowfr, esincos = esincos, verbose = verbose;, arst1 = arst1, arst2 = arst2, arst3 = arst3

  pm=p
  ;Roche limit stuff: mwd / (4d/3d * !dpi * rwd^3). For a 0.516 Msun, 0.0131 Rsun WD, rho =  323600 grams per cubic centimeter
  ; https://www.wolframalpha.com/input/?i=0.516+*+mass+of+sun+%2F+%284%2F3+*+pi+*+%280.0131+*+radius+of+sun%29%5E3%29
  ; For a 15 Mj planet, 1Rj, rho = 19.86 grams per cubic centimeter.  https://www.wolframalpha.com/input/?i=15+*+mass+of+jupiter+%2F+%284%2F3+*+pi+*+%28radius+of+jupiter%29%5E3%29
  ;
  ;

  ;  p0[0] = coeffsgtc[0]; u1gtc
  ;  p0[1] = coeffsgtc[1]; u2gtc
  ;  p0[2] = pnew[0]; period1
  ;  p0[3] = pnew[1]; tp1
  ;  p0[4] = 0; ecc1
  ;  p0[5] = !dpi/2d; omega1
  ;  p0[6] = cos(1.535); inclination1
  ;  p0[7] = pnew[6]; rprst1
  ;  p0[8] = 25; arst1
  ;  p0[9] = 1d-2; gtcscatter
  ;  p0[10] = coeffsspitzer[0]
  ;  p0[11] = coeffsspitzer[1]
  ;  p0[12] = 0d; spitzer dilution
  ;  p0[13] = 1d-1; spitzer scatter


  jitterpargtc = p[9]
  jitterparspit = p[13]
  dilutspit = p[12]

  ;G = 6.67259d-8
  lc = 0
  osamp = 6
  qs=0

  ; Stellar Parameters
  u1gtc = p[0]
  u2gtc =p[1]

  u1spitzer = p[10]
  u2spitzer = p[11]



  q1gtc = (u1gtc + u2gtc)^2
  q2gtc = 0.5* u1gtc /(u1gtc + u2gtc)
  q1spitzer = (u1spitzer + u2spitzer)^2
  q2spitzer = 0.5* u1spitzer /(u1spitzer + u2spitzer)
  
  if q1spitzer lt 0 or q2spitzer lt 0 or q1spitzer gt 1 or q2spitzer gt 1 then return, !values.d_infinity
  if q1gtc lt 0 or q2gtc lt 0 or q1gtc gt 1 or q2gtc gt 1 then return, !values.d_infinity
  if u1gtc lt 0 or u1spitzer lt 0 then return, !values.d_infinity
  if u1spitzer gt 0.2 or u2spitzer gt 0.3 then return, !values.d_infinity; a test for more physically reasonable values of spitzer LD coefficients. 
  ;if u1gtc lt 0 or u1gtc gt 1 or u2gtc lt -1 or u2gtc gt 1 then return, !values.d_infinity
  ;if u1spitzer lt 0 or u1spitzer gt 1 or u2spitzer lt -1 or u2spitzer gt 1 then return, !values.d_infinity


  ;dilutspitzer = p[12]
  



  ; Planet 1 transit parameters
  period1 = p[2]
  tp1 = p[3]
  ecc1 = p[4]
  omega1 = p[5]
  inclination1 = p[6]
  rprs1 = p[7]
  arst1 = p[8]

  if rprs1 lt 0 then begin
    if keyword_set(verbose) then print, 'rprs lt 0'
    return, !values.d_infinity

  endif

  if 1 - keyword_set(esincos) then begin ; this only works if we are fitting in inclination, not impact parameter
    if keyword_set(cosi) then begin
      inclination1 = acos(inclination1)
    endif
  end

  if keyword_set(esincos) then begin
    rtesinw1 = ecc1
    rtecosw1 = omega1
    ecc1 = rtesinw1^2 + rtecosw1^2
    omega1 = atan(rtesinw1, rtecosw1)



    if rtesinw1 gt 1  or rtecosw1 gt 1 then begin
      if keyword_set(verbose) then print, 'rtesinw1 gt 1  or rtecosw1 gt 1'
      return, !values.d_infinity
    endif
    if rtesinw1 lt -1  or rtecosw1 lt -1  then begin
      if keyword_set(verbose) then print, 'rtesinw1 lt -1  or rtecosw1 lt -1'
      return, !values.d_infinity
    endif

    tp1 = getperiastrontime(period1, omega1, ecc1, tp1)

    b1 = rprs1 - inclination1
    
    inclination1 = acos(b1/arst1 * (1 + ecc1 * sin(omega1)) / (1-ecc1^2))

    
  endif

  ;For the transit model,
  ;period = pmx[0]
  ;arst = pmx[1]
  ;tp = pmx[2]
  ;ecc = pmx[3]
  ;omega = pmx[4]
  ;inclination = pmx[5]
  ;u1 = pmx[6]
  ;u2=pmx[7]
  ;rprs = pmx[8]

  pm1 = dblarr(9)
  pm2 = dblarr(9)




  ;Planet 1 transit parameters
  pm1[0]=period1
  pm1[1]=arst1
  pm1[2]=tp1
  pm1[3]=ecc1
  pm1[4]=omega1
  pm1[5]=inclination1
  pm1[6]=u1gtc
  pm1[7]=u2gtc
  pm1[8]=rprs1



  pmgtc = pm1
  pmspit = pm1
  pmspit[[6,7]] = [u1spitzer, u2spitzer]




  b1  =  pm1(1) * cos(pm1(5)) * (1-pm1[3]^2) / (1 + pm1[3] * sin(pm1[4])) 
  
  rstar =  911367006.84204698; cm
  sma = pm1[1] * rstar
  mpl = 15 * 0.0009543; msun
  mst = 0.518 ; msun
  rpl = pm1[8] * rstar
  RL = (1-ecc1)*sma*0.46*(Mpl/Mst)^(1d/3d)
  if rpl gt rl then begin
    if keyword_set(verbose) then print, '  Roche lobe overflow'
    return, !values.d_infinity
  endif
  


  if jitterpargtc lt 0 then begin

    if keyword_set(verbose) then print, '  jitterpargtc lt 0'
    return, !values.d_infinity
  endif
  if jitterparspit lt 0 then begin
    if keyword_set(verbose) then print, ' jitterparspit lt 0'
    return, !values.d_infinity
  endif



  if inclination1 gt !dpi/2 then begin

    if keyword_set(verbose) then print, ' inclination1 gt !dpi/2'
    return, !values.d_infinity
  endif
  if b1 gt 1d + pm1[8] then begin
    
    if keyword_set(verbose) then print, ' b1 gt 1d + pm1[8]'
    return, !values.d_infinity; if the planet does not transit the star at all...
  endif



  if  pm1(1) lt 0 then begin
    if keyword_set(verbose) then print, 'pm1(1) lt 0 (negative a/r*)'
    return, !values.d_infinity ; exclude negative a/r* solutions
  endif

  ;if b gt 1d + p[8] then return, !values.d_infinity; if the planet does not transit the star at all...

  if pm1[3] lt 0 then begin
    if keyword_set(verbose) then print,  'pm1[3] lt 0 (ecc)'
    return, !values.D_INFINITY; Bound eccentricity
  endif
  if pm1[3] ge 1 then begin
    if keyword_set(verbose) then print, ' pm1[3] ge 1  (ecc)'
    return, !values.D_INFINITY; Bound eccentricity
  endif


  modelspit = transitmodel(pmspit,t=tspit,bint = 10d/24d/3600d, nos = 6, qs = qs)


  modelspit = (modelspit + dilutspit) / (1 + dilutspit)

  modelgtc = transitmodel(pmgtc,t=tgtc,bint = 10d/24d/3600d, nos = 6, qs = qs)
  
  

  ;  !p.multi = [0,2,1]
  ;  plot, tspit, fspit, ps=3,/yno
  ;  ;oplot, t * tbin, fbin, ps = 8, syms = 0.5, color = cgcolor('magenta')
  ;  oplot, tspit, modelspit, color = cgcolor('red')
  ;
  ;  plot, tgtc, fgtc, ps = 8,/yno
  ;  oplot, tgtc, modelgtc, color = cgcolor('blue')
  ;  !p.multi = 0
  ;  stop
  ;stop

  ;  if 1 - keyword_set(u1prior) then u1prior = 0.46290202118735219
  ;  if 1 - keyword_set(u2prior) then u2prior =    0.19287981261228282
  ;  pu1 = ((u1 - u1prior)^2 / (2 * 0.07^2));
  ;  pu2 = ((u2 -   u2prior)^2 / (2 * 0.07^2))

  G = 6.67259d-8
  rhostar = 3 * !dpi * arst1^3 / (G * (period1 * 24d * 3600d)^2)

  prho1 = ((rhostar - 323990.59)^2 / (2 * 54053.884^2));



  if keyword_set(dolowfr) then begin
    lowfrgtc = scaledpolyfit(tgtc, fgtc-modelgtc, 2)
    lowfrspit = scaledpolyfit(tspit, fspit-modelspit, 2)
  endif
  if 1 - keyword_set(dolowfr) then begin
    lowfrgtc = 0d
    lowfrspit = 0d
  endif


  newsigspit = sqrt(espit^2 + jitterparspit^2)
  newsiggtc = sqrt(egtc^2 + jitterpargtc^2)



  return, total(0.5 * (fspit - modelspit - lowfrspit)^2 / newsigspit^2 + alog(newsigspit)) + total(0.5 * (fgtc -modelgtc - lowfrgtc)^2 / newsiggtc^2 + alog(newsiggtc))  + prho1
  ;if keyword_set(mp) then return, (model - f) / e

end
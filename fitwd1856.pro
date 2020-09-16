pro fitwd1856

  rerun = 1
  dolowfr = 1
  ldfree = 0
  ldfreespit = 0
  testremovepoint = 0; check whether removing a single point ~4-5 sigma outlier makes a signifiant difference - changes dilution to -.004
  ld0 = 0
  ecc = 0
  syst = systime()
  closeall
  !p.multi = 0

  mwd =  (0.516  + jrandomn(seed, 10000) * 0.018) * 1.989d33 ; grams
  rwd = (0.0131+ jrandomn(seed, 10000) * 0.0003) * !rsunmeters * 100; cm

  rhowd = mwd / (4d/3d * !dpi * rwd^3)

  readcol, '~/wd1856_Spitzer-IRAC2_photometry.csv', tspitzer, fspitzer,goodpoints, xspitzer, yspitzer, format = 'D,D,D,D,D'
  restore,  '~/correctspitzerwd1856times.idl'
  oldbadtspitzer = tspitzer
  deltatspitzer = median(spbjd - tspitzer)
  tspitzer = spbjd
  ; print, 'NOT CORRECTING SPITZER TIMES'

  ;kspit = rebintesslightcurve(tspitzer, fspitzer, num = 5, cadence = 2d/60d)
  ;kspit2 = rebintesslightcurve(tspitzer, fspitzer, num = 5, cadence = 2d/60d, method1 = 1)
  ;  tspitzer = kspit.t
  ;  fspitzer = kspit.f
  ;  stop

  tspitzer = tspitzer - 2457000d
  ;oldbadtspitzer = oldbadtspitzer - 2457000d
  espitzer = sqrt(fspitzer)
  level = median(fspitzer[1200:*])
  fspitzer = fspitzer / level
  espitzer = espitzer/level


  x = where(tspitzer lt 1834.3+deltatspitzer and tspitzer gt 1834.27 +deltatspitzer)
  ;x2 = where(oldbadtspitzer  lt 1834.3 and oldbadtspitzer gt 1834.27)
  ;stop
  tspitzer = tspitzer[x]
  fspitzer = fspitzer[x]
  espitzer = espitzer[x]

  kspit = rebintesslightcurve(tspitzer, fspitzer, es = espitzer, num = 5, cadence = 2d/60d)

  tspitzer = kspit.t
  fspitzer = kspit.f
  espitzer = kspit.e
  

  
  

  espitzer = sqrt(fspitzer) * p2pscatter(fspitzer[170:250])
stop
  if testremovepoint then begin
    good  = where(tspitzer ne tspitzer[99])
    tspitzer = tspitzer[good]
    fspitzer = fspitzer[good]
    espitzer = espitzer[good]
    
    
  endif

  t0 = 2458708.978112d -2457000d
  p = 1.4079342d
  ;stop


  ;readcol, '~/GTC_WD1856_IRAF_Photometry.txt', tgtc, fgtc, egtc
  restore, '~/GTC_WD1856_IRAF_Photometry_BJDCORRECTED.idl'
  ;tgtc = utc2bjd(tgtc, 284.41393, 53.509250)

  tgtc = tgtc - 2457000d



  x = where(tgtc gt 1779.378)
  egtc = sqrt(fgtc) * p2pscatter(fgtc[x[0:90]])



  x = where(tgtc lt 1779.382 and tgtc gt  1779.369)
  tgtc = tgtc[x]
  fgtc = fgtc[x]
  egtc = egtc[x]


  pfit = 1.40797303919d
  t0fit = 1708.97643079523



  if 1 - dolowfr then begin
    tmgtc = maketmodn(tgtc, pfit, t0fit)

    lowfrgtc2 = scaledpolyfit(tgtc, fgtc, 1, include = where(abs(tmgtc) gt 5d/24d/60d))
    fgtc = fgtc / lowfrgtc2

    tmspit = maketmodn(tspitzer, pfit, t0fit)
    lowfrspit2 = scaledpolyfit(tspitzer, fspitzer, 1, include = where(abs(tmspit) gt 5d/24d/60d))
    fspitzer = fspitzer / lowfrspit2

  endif


  pnew = [ p    ,   t0   ,  0.078789159  ,    0.97942000   ,    0.0000000   ,  -0.50000000  ,    0.16339442]



  ;coeffgtc = 0.4475; linear for lsst g band; this is pre-referee
  ;coeffsgtc = [0.0585, 0.516]; quadratic for LSST g band; this is pre-referee
  coeffsgtc = [0.0477, 0.5159]

  if ld0 then begin
    coeffsgtc = [0,0]
    coeffsspitzer = [0,0]
  endif

  ;coeffsspitzer = quadld(7.75, 4300, 0, 'Spit45'); this doesn't actually extrapolate to log(g) of 7.75, but it's pretty minor changes all the way out
  ;coeffsspitzer = quadld(7.75, 4300, 0, 'Spit45'); pre-referee. this gives [0.05874], [0.17826]
  coeffsspitzer = quadld(5, 4710, -4.5, 'Spit45')
  coeffsspitzer[0] = 0; to keep it physical
  
  ;  q1 = (coeffs[0] + coeffs[1])^2
  ;  q2 = 0.5* coeffs[0] /(coeffs[0] + coeffs[1])

  p0 = dblarr(14)

  p0[0] = coeffsgtc[0]; u1gtc
  p0[1] = coeffsgtc[1]; u2gtc
  p0[2] = pnew[0]; period1
  p0[3] = pnew[1]; tp1
  p0[4] = 0; ecc1
  p0[5] = !dpi/2d; omega1
  p0[6] = cos(!dpi/2 - 2.3d-2); inclination1
  p0[7] = 8; rprst1
  p0[8] = 350; arst1
  p0[9] = 1d-2; gtcscatter
  p0[10] = coeffsspitzer[0]
  p0[11] = coeffsspitzer[1]
  p0[12] = 0d; spitzer dilution
  p0[13] = 1d-1; spitzer scatter

  restore, '~/wd1856mcmcstartingpars.idl'
  p0 = transpose(pbest)

  p0[2] =  1.40793806257d
  p0[3] = 1746.992434d; uncomment to fix spitzer times


  if ld0 then begin
    p0[[0,1,10,11]] = 0
  endif


  widthp0 = dblarr(n_elements(p0))
  widthp0[0] = 0; u1gtc
  widthp0[1] = 0; u2gtc
  widthp0[2] =0.05 * 1d/24d/3600d; period1
  widthp0[3] = 20d/24d/3600; tp1
  widthp0[4] = 0; ecc1
  widthp0[5] = 0; omega1
  widthp0[6] = .00003; inclination1
  widthp0[7] = p0[7]/30; rprst1
  widthp0[8] = 3; arst1
  widthp0[9] = 1d-2; gtcscatter
  widthp0[10] = 0d; u1spitzer
  widthp0[11] = 0d; u2spitzer
  widthp0[12] = 1d-2; spitzer dilution
  widthp0[13] = 1d-2; spitzer scatter

  if ldfree then widthp0[[0,1,10,11]] = 0.2
  if ldfreespit then widthp0[[10,11]] = 0.2

  if ecc then begin
    p0[[4,5]] = [0,0]
    widthp0[[4,5]] = [0.5,0.5]

    impact0 = p0[8] * p0[6]; change inclination into impact parameter
    p0[6] = p0[7] - impact0; fit in the difference of Rp/R* and impact parameter
    widthp0[6] = 0.03
    widthp0[7] = 3

  endif

  width = widthp0


  ;tmodn1 = maketmodn(t, p0[2], p0[3])
  ;smalltmodn = maketmodn(t, p0[2], p0[3])

  extra = {tspit:tspitzer, fspit:fspitzer, espit:espitzer,tgtc:tgtc, fgtc:fgtc, egtc:egtc, cosi:1, dolowfr:dolowfr, esincos:ecc}

  junk = mpwd1856(p0, tspit = extra.tspit, fspit=extra.fspit, espit=extra.espit, tgtc = extra.tgtc, fgtc=extra.fgtc, egtc=extra.egtc,$
    modelspit = modelspit, modelgtc = modelgtc, cosi = extra.cosi, b1 = b1, dolowfr=extra.dolowfr, esincos = extra.esincos)

  print, b1

  plot, tspitzer, fspitzer, ps=3,/yno
  ;oplot, t * tbin, fbin, ps = 8, syms = 0.5, color = cgcolor('magenta')
  oplot, tspitzer, modelspit, color = cgcolor('red')
  figure
  plot, tgtc, fgtc, ps = 8,/yno
  oplot, tgtc, modelgtc, color = cgcolor('blue')



  if rerun then begin
    nlink =  200000
    nburnin = 50000
    nwalk = 50

    filename = '/data/k2/misc/wd1856transitpars-'+syst+'-nlink'+trim(nlink)+'-nwalk'+trim(nwalk)+'-nburn'+trim(nburnin)
    if ld0 then filename = filename + '-ld0'
    if 1- dolowfr then filename = filename + '-nolowfr'
    if ecc then filename = filename + '-ecc'
    if ldfree then filename = filename + '-ldfree'
    if ldfreespit then filename = filename + '-ldfreespit'
    if testremovepoint then filename = filename + '-outliercliped'

    filename = filename + '.idl'
    filename = repstr(filename, ' ', '_')
    print, 'File will be written to: ' + filename

    thispnew = affinemcmcfit('mpwd1856',p0, functargs = extra, $
      perror = perror, nwalkers = nwalk, nlink = nlink,nburnin = nburnin, chains = chains, width = width,$
      allneglogl = neglogl, whichwalker = whichwalker, whichlink = whichlink, walkers = walkers)

    walkersout = walkers

    save, thispnew, perror, chains, widthp0, neglogl, whichwalker, whichlink, walkersout,nwalk, nlink, nburnin, walkers, $
      filename = filename
    stop
  end
  if 1 - rerun then begin
    ;restore, '/Users/andrew/Documents/UT/TESS/wd1856/transitpars-Mon Jan 20 16:18:10 2020-nlink10000-nwalk50-nburn1000.idl'
    if 1 - ld0 and 1-ecc then begin

      restore, '/Users/andrew/Documents/UT/TESS/wd1856/transitpars-Mon Feb 10 09:57:53 2020-nlink10000-nwalk50-nburn1000.idl'
      ;restore, '/Users/andrew/Documents/UT/TESS/wd1856/transitpars-Sun Feb  2 16:44:16 2020-nlink10000-nwalk50-nburn1000.idl'
      ;restore, '/Users/andrew/Documents/UT/TESS/wd1856/transitpars-Sun Feb  2 16:13:49 2020-nlink10000-nwalk50-nburn1000.idl'
      ;restore, '/Users/andrew/Documents/UT/TESS/wd1856/transitpars-Mon Jan 20 17:36:21 2020-nlink100000-nwalk50-nburn10000.idl'
    end
    if ld0 then restore, '~/Documents/UT/TESS/wd1856/transitpars-Tue Jan 21 15:10:15 2020-nlink10000-nwalk50-nburn1000-ld0.idl'
    if ecc then restore, '/data/k2/misc/wd1856transitpars-Mon Feb 10 19:08:24 2020-nlink200000-nwalk50-nburn10000-ecc.idl'
  endif
  excludewalkers, chains, neglogl, whichwalker, pbest = pbest, pnew = pnew, perror = perror
  ;save, pbest, filename = '~/Documents/UT/TESS/wd1856/mcmcstartingpars.idl'


  junk = mpwd1856(pbest, tspit = extra.tspit, fspit=extra.fspit, espit=extra.espit, tgtc = extra.tgtc, fgtc=extra.fgtc, egtc=extra.egtc,$
    modelspit = modelspit, modelgtc = modelgtc, cosi = extra.cosi, b1 = b1, lowfrgtc=lowfrgtc, lowfrspit = lowfrspit, dolowfr = extra.dolowfr, esincos = extra.esincos)


  print, b1
  closeall
  figure, xsize = 1400
  !p.multi = [0,2,1]
  plot, tspitzer, fspitzer - lowfrspit, ps=8,/yno
  ;oplot, t * tbin, fbin, ps = 8, syms = 0.5, color = cgcolor('magenta')
  oplot, tspitzer, modelspit, color = cgcolor('red')

  plot, tgtc, fgtc - lowfrgtc, ps = 8,/yno
  oplot, tgtc, modelgtc, color = cgcolor('blue')
  !p.multi = 0

  tmgtc = maketmodn(tgtc, pbest[2], pbest[3])
  tmspitzer = maketmodn(tspitzer, pbest[2], pbest[3])

  ;figure
  ;  psopen, '~/Documents/UT/TESS/wd1856/gtcspitzer.eps', xsize = 6, ysize = 6, /inches, /color,/encapsulated
  ;
  ;  plot, tmgtc * 24d * 60d, fgtc - lowfrgtc, ps = 8, syms = 1.5,/nodata, xtitle = 'Minutes from Midtransit', ytitle = 'Relative Brightness', xr = [-7,7],/xs
  ;  oploterror, tmspitzer * 24d * 60d, fspitzer - lowfrspit,sqrt(espitzer^2 + pbest[13]^2), ps = 8,syms = 1.5, color = cgcolor('black'),/nohat
  ;  oplot, tmspitzer * 24d * 60d, fspitzer - lowfrspit, ps = 8,syms = 1.0, color = cgcolor('maroon')
  ;
  ;  oplot, tmspitzer* 24d * 60d, modelspit, color = cgcolor('pur8'), thick = 8
  ;
  ;  oploterror, tmgtc * 24d * 60d, fgtc-lowfrgtc, egtc, ps = 8,syms = 1.5
  ;  oplot, tmgtc * 24d * 60d, fgtc-lowfrgtc, ps = 8, color = cgcolor('dodger blue')
  ;  oplot, tmgtc* 24d * 60d, modelgtc, color = cgcolor('red'), thick = 8
  ;  oplot, tmspitzer* 24d * 60d, modelspit, color = cgcolor('pur8'), thick = 8
  ;
  ;
  ;
  ;  psclose
  ;
  ;  psopen, '~/Documents/UT/TESS/wd1856/gtcspitzersidebyside.eps', xsize = 9, ysize = 4, /inches, /color,/encapsulated
  ;  !p.multi = [0,2,1]
  ;  plot, tmgtc * 24d * 60d, fgtc - lowfrgtc, ps = 8, syms = 1.5,/nodata, xtitle = 'Minutes from Midtransit', ytitle = 'Relative Brightness', xr = [-7,7],/xs,/ys, yr = [0.2, 1.2]
  ;  oploterror, tmgtc * 24d * 60d, fgtc-lowfrgtc, egtc, ps = 8,syms = 1.5
  ;  oplot, tmgtc * 24d * 60d, fgtc-lowfrgtc, ps = 8, color = cgcolor('grey')
  ;  oplot, tmgtc* 24d * 60d, modelgtc, color = cgcolor('pur8'), thick = 8
  ;
  ;  plot, tmspitzer * 24d * 60d, fspitzer - lowfrspit,ps = 8, syms = 1.5,/nodata, xtitle = 'Minutes from Midtransit', ytitle = 'Relative Brightness', xr = [-7,7],/xs,/ys, yr = [0.2, 1.2]
  ;  oploterror, tmspitzer * 24d * 60d, fspitzer - lowfrspit,sqrt(espitzer^2 + pbest[13]^2), ps = 8,syms = 1.5, color = cgcolor('black'),/nohat
  ;  oplot, tmspitzer * 24d * 60d, fspitzer - lowfrspit, ps = 8,syms = 1.0, color = cgcolor('grey')
  ;
  ;  oplot, tmspitzer* 24d * 60d, modelspit, color = cgcolor('pur8'), thick = 8
  ;  !p.multi = 0
  ;  psclose
  rwd = (0.0131+ jrandomn(seed, n_elements(chains(*,7))) * 0.0003) * !rsunmeters * 100; cm

  radius = rwd * chains(*,7) / 100 / !rearthmeters

  ecc = chains(*,4)^2 + chains(*,5)^2
  omega = atan(chains(*,4), chains(*,5))

  habitablezone, 4785, 8.0716433e-05, consinner=consinner, consouter = consouter, optinner = optinner, optouter = optouter

  auincm = 1.496d13


  stop

end
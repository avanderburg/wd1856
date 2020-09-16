pro wd1856fluxes
  closeall
  usesonora = 0
  usespitzerflux = 1

  ;    IDL> print, getpercentile(chains(*,12), 95)
  ;    % Compiled module: GETPERCENTILE.
  ;    0.065817240
  ;    IDL> print, getpercentile(chains(*,12), 99)
  ;    0.089895203
  ;    IDL> print, getpercentile(chains(*,12), 99.7)
  ;    0.10566167

  ;restore, '/Users/andrew/Documents/UT/TESS/wd1856/transitpars-Sun Feb  2 16:44:16 2020-nlink10000-nwalk50-nburn1000.idl'; with old spitzer errors
  ;restore, '/Users/andrew/Documents/UT/TESS/wd1856/transitpars-Mon Feb 10 09:57:53 2020-nlink10000-nwalk50-nburn1000.idl'; with new (correct) spitzer errors
  ;restore, '/Users/andrew/Documents/UT/TESS/wd1856/transitpars-Mon Feb 10 12:52:19 2020-nlink10000-nwalk50-nburn1000.idl'; new spitzer errors, no lowfr
  ;restore, '/Users/andrew/Documents/UT/TESS/wd1856/wd1856transitpars-thinned-Feb11-11:06:15.idl' - used for the initial submission
  ;restore, '/Users/andrew/Documents/UT/TESS/wd1856/wd1856transitpars-ThuMar5_13:04:41_2020-nlink500000-nwalk50-nburn100000-thinned.idl'; had shifted BJDs
  restore, '/Users/andrew/Documents/UT/TESS/wd1856/wd1856transitpars-Mon-May-4-19:27:37-2020-nlink200000-nwalk50-nburn100000-THIN.idl'
  dg0 = 1; use only samples with positive dilution parameter
  if dg0 then begin
    pos = where(chains(*,12) gt 0)
    chains = chains[pos, *]
  endif
  
  
  confidencelevel = 95; %
 ; confidencelevel = 99.7
 ; confidencelevel = 68

  ;limit = 0.10566167; 99.7%
  ;limit = 0.089895203; 99%
  ;limit = 0.065817240; 95%
  ;limit = 0.0544; %90%
  ;limit = 0.027015570; 68%


  if 1-usesonora then begin
    restore, '/Users/andrew/Documents/UT/rvexomoons/BD_MESA_complete_set.data.idl'

    g = 6.67259d-8
    rightage = where(age gt 5.6d9 and age lt 13.8d9 and mass lt 0.094)
    ;rightage = where(age gt 5.6d9 and age lt 10d9 and mass lt 0.094)
    ;rightage = where(age gt 2d8)
    ;rightage = where(age lt 5.6d9 and age lt 13.82d9 and mass lt 0.03)

    ;rightage = where(age gt 7d9 and age lt 8d9)

    mesarad = rad[rightage]
    mesamass = mass[rightage]
    mesaage = age[rightage]
    mesateff = teff[rightage]

    msun = 1.989d33
    rsun = 69.634d9
    mesalogg = alog10(g * mesamass * msun / (rsun * mesarad)^2)
  end
  if usesonora then begin
    ;readcol, '/Users/andrew/Documents/UT/TESS/wd1856/sonoratables_2018/nc+0.0_co1.0_age', age ,  mass , lum  ,teff  ,logg , rad, format ='D,D,D,D,D,D', delimiter = ' '
    readcol, '/Users/andrew/Documents/UT/TESS/wd1856/sonoratables_2018/nc+0.0_co1.0_lbol', loglum,  mass,   age,  Teff,  logg , rad, format = 'D,D,D,D,D,D', delimiter = ' '

    g = 6.67259d-8
    age = age * 1d9
    rightage = where(age gt 5.6d9 and age lt 13.82d9 and mass lt 0.094)
    rightage = where(age gt 5.6d9 and age lt 13.82d9 and mass lt 0.094)
    ;rightage = where(age lt 5d20)

    ;rightage = where(age lt 5.6d9 and age lt 13.82d9 and mass lt 0.03)

    ;rightage = where(age gt 7d9 and age lt 8d9)

    mesarad = rad[rightage]
    mesamass = mass[rightage]
    mesaage = age[rightage]
    mesateff = teff[rightage]

    msun = 1.989d33
    rsun = 69.634d9
    mesalogg = alog10(g * mesamass * msun / (rsun * mesarad)^2)



  endif


  closeall


  rsun = 6.995d10 ; cm
  parsec = 3.08567758128d18
  distance = 1000/(40.3983d + randomn(seed, 1000) * 0.0704537) * parsec


  readcol, '/Users/andrew/Documents/UT/TESS/wd1856/WD1856SED.csv', ra, dec, cat, freq, newflux, newerrflux, band, format = 'D,D,A,D,D,D,A',comment = '#'
  ; fluxes are in janskys

  freq = freq * 1d9 ; into hertz
  ; lambda = c/freq
  lambda = 3d8 / freq; into meters

  readcol, '/Users/andrew/Documents/UT/TESS/wd1856/WD1856_synth_spec.txt', wlsimon, hnusimon
  fnusimon = hnusimon * 4 * !dpi * ((0.0131 * rsun)/(median(distance)))^2

  ;flambda = fnu * dnu/dlambda
  flambdasimon = fnusimon *  (3d18/wlsimon^2)

  plot, lambda * 1d10, newflux, ps = 8,xtitle = 'wavelength (angstrom)', ytitle = 'Fnu (Jansky)'
  oploterror, lambda * 1d10, newflux, newerrflux, ps = 8, color = cgcolor('red')
  oplot, wlsimon, fnusimon * 1d23

  readcol,'/Users/andrew/Documents/UT/TESS/wd1856/sonora_flux_table.txt', Teff,  logg ,  mass , RRsun,   Y,   logKzz, Ymag ,Z, J, H, K,Lprime , Mprime, twomassJ, twomassH, twomassKs, keckKs, keckLprime ,keckMs ,gprime,rprime, iprime, zprime , irac36, irac45, irac57, irac79, W1, W2, W3, W4, delimiter = ' ', comment = '#'






  ;  subset = where(logg gt 4.4 and logg lt 4.6)
  ;  teff = teff[subset]
  ;  logg = logg[subset]
  ;  mass = mass[subset]
  ;  rrsun = rrsun[subset]
  ;  irac45 = irac45[subset]

  mass = mass * 0.0009543
  distancepc = 1000/(40.3983d)

  sedwdflux45 = interpol(fnusimon * 1d23, wlsimon, 4.5d4)
  if 1 - usespitzerflux then begin
    wdflux45 = interpol(fnusimon * 1d23, wlsimon, 4.5d4)
    ;bdflux = wdflux45 * limit
    bdflux = getpercentile(wdflux45 * chains(*,12), confidencelevel)
  end
  ;wdflux45 = mean(fnusimon[where(wlsimon lt 5d4 and wlsimon gt 4d4)]) * 1d23

  if usespitzerflux then begin
    wdflux45 = 173d-6 + jrandomn(seed, n_elements(chains(*,12))) * 10d-6
    ;bdflux = wdflux45 / (1 + 1/limit)
    bdflux = getpercentile(wdflux45 / (1 + 1/chains(*,12)), confidencelevel)

  end



  ;oplot, [4.5d4],[sedwdflux45],ps = 8, color = cgcolor('blue')
  oploterror, [4.5d4],[173d-6], [10d-6], color = cgcolor('magenta'), ps = 8






  figure
  irac45fluxes = 10^(irac45) * 1d-3 *  (10d/distancepc)^2; should be in jansky
  ;  plot, teff, irac45fluxes,/yl, ps = 8, xtitle = 'teff', ytitle = 'IRAC 4.5 micron fluxes'
  ;  oplothoriz, wdflux45 *limit

  ;  x = where(irac45fluxes lt wdflux45 * limit)
  ;  mesa45fluxesnearest = dblarr(n_elements(mesaage))
  ;  for i = 0, n_elements(mesaage)-1 do begin
  ;    dis = sqrt((mesalogg[i] - logg)^2 + (mesateff[i] - teff)^2 + (mesamass[i] - mass)^2)
  ;    junk = min(dis, minind)
  ;    mesa45fluxesnearest[i] = irac45fluxes[minind]
  ;    ;print, logg[minind]
  ;  endfor
  ;
  ;
  ;  plot, mesamass /0.0009543, mesa45fluxesnearest, ps = 3,xr = [10,20]
  ;  hline, bdflux



  triangulate, logg, teff, triangles
  gridfluxes = trigrid(logg, teff, irac45fluxes, triangles, xgrid = xg, ygrid = yg, xout = xout, yout = yout)

  gridlogg = interpol(findgen(51), xg, mesalogg)
  gridteff = interpol(findgen(51), yg, mesateff)

  mesaflux45 = interpolate(gridfluxes, gridlogg, gridteff)

  ;mesaflux45test = griddata(logg, teff, irac45fluxes, xout = mesalogg, yout = mesateff,method = 'NearestNeighbor', triangles = triangles)

  plot, mesamass /0.0009543, mesaflux45, ps = 3,xr = [10,20]
  hline, bdflux

  x = where(mesaflux45 lt bdflux and mesamass lt 0.05)
  print, trim(max(mesamass[x]) / 0.0009543) + ' Mjup'
  print, 'At '+ trim(confidencelevel) + '% confidence, the mass must be less than ' + trim(max(mesamass[x])/0.0009543) + ' M_j.'

  print,  trim(max(mesateff[x])) + ' K'

  figure
  plot, mesaage/1d9, mesamass, ps = 3, xtitle = 'Age (gyr)', ytitle = 'Mass'
  oplot, mesaage[x]/1d9, mesamass[x], ps = 3, color = cgcolor('red')
  y = where(mesamass lt 0.013 and mesaflux45 lt bdflux)
  oplot, mesaage[y]/1d9, mesamass[y], ps = 3, color = cgcolor('blue')


  ;
  ;    triangulate, logg, rrsun, triangles
  ;    gridfluxes = trigrid(mass, rrsun, irac45fluxes, triangles, xgrid = xg, ygrid = yg, xout = xout, yout = yout)
  ;
  ;    gridmass = interpol(findgen(51), xg, mesamass)
  ;    gridradius = interpol(findgen(51), yg, mesarad)
  ;
  ;    mesaflux45 = interpolate(gridfluxes, gridmass, gridradius)
  ;
  ;    mesaflux45nearest = dblarr(n_elements(mesaage))
  ;
  ;    for i = 0, n_elements(mesaage)-1 do begin
  ;      dis = sqrt()
  ;    endfor
  ;
  ;
  ;


  ;    test = wraptriinterp2d( logg, teff, irac45fluxes, logg+0.0001, teff + 0.0001)
  ;    test2 = griddata(logg, teff, irac45fluxes, xout = logg+0.0001, yout = teff + 0.0001)
  ;    figure, xsize = 1200
  ;    !p.multi = [0,2,1]
  ;    heatmap,  mesaage, mesamass, alog(mesaflux45), title = 'original triangulation'
  ;    heatmap,  mesaage, mesamass, alog(mesaflux45test),  title = 'test';, min = min(alog(mesaflux45)), max = max(alog(mesaflux45))
  ;    !p.multi = 0
  stop


end
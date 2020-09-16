

function transitmodel, p, t = t, bint = bint, nos = nos, lc=lc, sc=sc, qs = qs, impact = b, zplanet = zplanet

  if keyword_set(lc) then begin
    bint = 1764.944d/(24d * 3600d);used to be 29.4244d/(24d * 60d)
    if 1 - keyword_set(nos) then nos = 4
  end
  if keyword_set(sc) then begin
    bint = 58.34d/(24d * 3600d); used to be 58.85 
    if 1 - keyword_set(nos) then nos = 4
  end
  
  
  
  
  period = p[0]
  arst = p[1]
  tp = p[2]
  ecc = p[3]
  omega = p[4]
  inclination = p[5]
  u1 = p[6]
  u2=p[7]
  rprs = p[8]
  
  if keyword_set(qs) then begin
    q1 = u1 ; sample the mcmc chain in q-space instead of quadratic law coefficients from kipping's limb darkening sampling paper
    q2 = u2
    
    u1 = 2 * sqrt(q1) * q2
    u2 = sqrt(q1) * (1 - 2 * q2)
  end
  
  flux = dblarr(n_elements(t))
  
  if 1 - keyword_set(bint) then begin
  
    b = exofast_getb(t, i=inclination,a=arst, tperiastron=tp,period=period,e=ecc,omega=omega, x=xplanet,y=yplanet,z=zplanet)
    b = b * (zplanet ge 0) + 10 * (zplanet lt 0)
    exofast_occultquad, b, u1, u2, rprs, fluxi
    return, fluxi
    
  end
  
  
  if keyword_set(bint) then begin
    if 1 - keyword_set(nos) then nos = 4
    
    dt = bint / nos
    
    fluxes = dblarr(n_elements(t), nos+1)
    ;tzeros = dblarr((nos) + 1)
    for i = 0, nos do begin
      ti = t  - bint/2. + (i) * dt
      ;tzeros(i) = ti(0)
      b = exofast_getb(ti, i=inclination,a=arst, tperiastron=tp,period=period,e=ecc,omega=omega, x=xplanet,y=yplanet,z=zplanet)
      b = b * (zplanet ge 0) + 10 * (zplanet lt 0)
      
      exofast_occultquad, b, u1, u2, rprs, fluxi      
      fluxes(*,i) = fluxi
    ;flux = flux + fluxi
    end
    ;stop
    ;flux = flux / (nos + 1d)
    flux = total(fluxes, 2) - fluxes(*,0) / 2d - fluxes(*,nos) / 2d ; trapezoidal integration
    flux = flux / nos
    return, flux
  end
end




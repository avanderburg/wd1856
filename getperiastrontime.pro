function getperiastrontime, period, wp, e, tt

  ;calculates the time of transit given the period, angle of periastron (wp), eccentricity (e), and time of transit (tt).
  
  f = !dpi/2 - wp
  ea = 2 * atan(sqrt((1-e)/(1+e)) * tan(f/2))
  te = period / (2 * !dpi) * (ea - e * sin(ea))
  
  ;tt = tp + te
  tp = tt - te
  
  return, tp
  
  
end
pro excludewalkers, chains, neglogl, whichwalker, oldchains = oldchains, oldneglogl = oldneglogl, oldwhichwalker = oldwhichwalker, pbest = pbest, perror = perror, pnew = pnew

  minnll = min(neglogl, bestind)
  nwalkers = max(whichwalker)
  excluded = intarr(nwalkers+1)
  ksexcluded = excluded
  nlinksdone = n_elements(neglogl)
  goodchains = [-1]
  ksgoodchains = [-1]
  oldchains = -1
  oldneglogl = -1
  oldwhichwalker = -1
  bestwalker = whichwalker(bestind)
  bestsample = where(whichwalker eq bestwalker)
  for ni = 0, nwalkers do begin
    thisone = where(whichwalker eq ni)
    prob = exp(double(minnll) - min(neglogl(thisone)))
    
    kstwo, neglogl(bestsample), neglogl(thisone), deviation, ksprob
    if deviation gt 0.5 then ksexcluded(ni) = 1 ; if the KS test gives significant deviations (cumulative probability differs by greater than 0.3),reject walker
    if deviation le 0.5 then ksgoodchains = [ksgoodchains, thisone]
    
    if prob lt 1d/(100d * double(nlinksdone)) then excluded(ni) = 1
    if prob ge 1d/(100d * double(nlinksdone)) then goodchains = [goodchains, thisone]
  end
  if total(excluded) gt 0 and total(excluded) ne nwalkers+1 then begin
  
    nexcluded = total(excluded)
    
    goodchains = goodchains[1:*]
    ksgoodchains = ksgoodchains[1:*]
    
    if total(ksexcluded) lt floor(nwalkers/2) then begin
      goodchains = cgsetintersection(goodchains, ksgoodchains)
      nexcluded = total(excluded or ksexcluded)
    end
    
    print, ''
    print, strtrim(fix(nexcluded),2) + ' walkers rejected due to time spent in a bad miminum.'
    print, ''
    
    oldchains = chains
    oldneglogl = neglogl
    oldwhichwalker = whichwalker
    chains = chains(goodchains, *)
    whichwalker = whichwalker(goodchains)
    neglogl = neglogl(goodchains)
    pnew = dindgen(nrows(chains))
    for ni = 0, n_elements(pnew)-1 do begin
      pnew[ni] = median(chains[*,ni])
      perror[ni] = stdev(chains[*,ni])
    end
  end

  minnll = min(neglogl, bestind)

  pbest = chains(bestind, *)









end
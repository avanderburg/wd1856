function scaledpolyfit, t,y,deg,yfit,yfit1,covariance=cov,weight=w, includeindsin = includeindsin

  t2 = (t - min(t)) / (max(t) - min(t))
  if 1 - keyword_set(includeindsin) then begin
    includeinds = indgen(n_elements(t))
  endif
  if keyword_set(includeindsin) then begin
    includeinds = includeindsin
  endif
  c = polyfit(t2[includeinds],y[includeinds],deg,yfit,yfit1,covariance=cov,weight=w)
  
  
  
  return, polyexpand(t2, c)



end
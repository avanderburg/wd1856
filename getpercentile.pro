function getpercentile, data, percentile, onesigma = onesigma
  a = sort(data)
  sd = data(a)
  if 1 - keyword_set(onesigma) then return, sd(floor(double(percentile) / 100d * (n_elements(sd) - 1)))
  if keyword_set(onesigma) then return, [sd(floor(15.865d / 100d * (n_elements(sd) - 1))), sd(floor(84.135d / 100d * (n_elements(sd) - 1)))]
end
function p2pscatter, f ; calculates the point to point scatter and scales it to estimate the standard deviation of the time series. This is robust to slow trends and outliers. 

  return, median(abs(f[1:*] - f[0:n_elements(f)-1])) * 1.48 /sqrt(2)

end



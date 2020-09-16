function polyexpand, x, pars
  y = fltarr(n_elements(x))
  for i = 0, n_elements(pars)-1 do begin
    y += x^i * pars[i]
  end
  return, y
end
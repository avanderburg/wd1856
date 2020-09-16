function maketmodn, t, p, e
  tmodn = (t mod p) - (e mod p)
  tmodn = tmodn + p * (tmodn lt -0.5 * p) - p * (tmodn gt 0.5 * p)
  return, tmodn
end
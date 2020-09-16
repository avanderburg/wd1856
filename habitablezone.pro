pro habitablezone, teff, lum, consinner=consinner, consouter=consouter, optinner=optinner, optouter=optouter, sconsinner=sconsinner, sconsouter=sconsouter, soptinner=soptinner, soptouter=soptouter
; uses ravi kopparapu's habitable zone calculator polynomial coefficients to calculate HZ boundaries given Teff and Luminosity.
; inputs in terms of solar luminosity and temp in kelvin
seffsun = dblarr(8)
a = dblarr(8)
b = dblarr(8)
c = dblarr(8)
d = dblarr(8)

;seffsun(1) = 1.7763
;seffsun(2) = 1.0385
;seffsun(3) = 1.0146
;seffsun(4) = 0.3507
;seffsun(5) = 0.3207
;seffsun(6) = 0.2484
;seffsun(7) = 0.5408
;
;a(1) = 1.4335e-4
;a(2) = 1.2456e-4
;a(3) = 8.1884e-5 
;a(4) = 5.9578e-5
;a(5) = 5.4471e-5
;a(6) = 4.2588e-5
;a(7) = 4.4499e-5
;
;b(1) = 3.3954e-9 
;b(2) = 1.4612e-8
;b(3) = 1.9394e-9
;b(4) = 1.6707e-9
;b(5) = 1.5275e-9
;b(6) = 1.1963e-9
;b(7) = 1.4065e-10
;
;c(1) = -7.6364e-12
;c(2) = -7.6345e-12
;c(3) = -4.3618e-12
;c(4) = -3.0058e-12
;c(5) = -2.7481e-12
;c(6) = -2.1709e-12
;c(7) = -2.2750e-12
;
;d(1) = -1.1950e-15
;d(2) = -1.7511E-15
;d(3) = -6.8260e-16
;d(4) = -5.1925e-16
;d(5) = -4.7474e-16
;d(6) = -3.8282e-16
;d(7) = -3.3509e-16

seffsun(1) = 1.776d  
seffsun(2) = 1.107d
seffsun(3) = 0.356d
seffsun(4) = 0.320d
seffsun(5) = 1.188d
seffsun(6) = 0.99d 

a(1) = 2.136d-4
a(2) = 1.332d-4
a(3) = 6.171d-5
a(4) = 5.547d-5
a(5) = 1.433d-4
a(6) = 1.209d-4

b(1) = 2.533d-8
b(2) = 1.580d-8
b(3) = 1.698d-9
b(4) = 1.526d-9
b(5) = 1.707d-8
b(6) = 1.404d-8

c(1) = -1.332d-11
c(2) = -8.308d-12
c(3) = -3.198d-12
c(4) = -2.874d-12
c(5) = -8.968d-12
c(6) = -7.418d-12

d(1) = -3.097d-15
d(2) = -1.931d-15
d(3) = -5.575d-16
d(4) = -5.011d-16
d(5) = -2.084d-15
d(6) = -1.713d-15



tstar = teff - 5780.0d
seff = dblarr(n_elements(teff), 8)
distance = dblarr(n_elements(teff), 8)
for i = 1,7 do begin
  seff[*,i] = seffsun(i) + a(i) * tstar + b(i) * tstar^2d + c(i) * tstar^3d + d(i) * tstar^4d
  distance[*,i] = sqrt(lum / seff[*,i])
end

;consinner = distance[*,3]
;consouter = distance[*,4]
;optinner = distance[*,1]
;optouter = distance[*,5]

consinner = distance[*,5]
consouter = distance[*,3]
optinner = distance[*,1]
optouter = distance[*,4]

sconsinner = seff[*,5]
sconsouter = seff[*,3]
soptinner = seff[*,1]
soptouter = seff[*,4]


end
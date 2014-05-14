Code to compute redshift space anisotropic P(k) 
-----------------------------------------------

The Kaiser term P_kaiser(k,mu) is 

(1+f mu^2)^2 P_dd (k)   (Linear)
P_dd (k) + 2 f mu^2 P_dt (k) + f^2 mu^4 P_tt (Non-linear)

where f is the growth rate. The redshift space power spectrum is 

P(k,mu) = D_FoG [k mu f sigma_v] P_kaiser(k,mu)

We asume a "finger of God" term 

D_FoG[x] = exp(-x^2)

In addition we include corrections according 

http://jp.arxiv.org/pdf/1006.0699v1

so that 

P(k,mu) = D_FoG [k mu f sigma_v] {P_dd (k) + 2 f mu^2 P_dt (k) + f^2 mu^4 P_tt + A(k,mu,f) + B(k,mu,f)}

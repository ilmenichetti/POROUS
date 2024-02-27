ky=0.8
ko=0.04
  kmix=0.02
  e=0.15
  F_prot=0.1
  phi_mac=0.04
  clay=0.2
  Delta_z_min=14.2
  gamma_o=1.2
  gamma_m=2.7
  proportion=NULL
  phi_min=0.35
  f_text_mic=0.5069
  f_agg=3
  init=c(My_mes=0.132, Mo_mes=0.507, My_mic=0.284, Mo_mic=0.579)

A_a=0.001
k_t=0.1
MY_mes=seq(0,20,0.1)
MO_mes=2
(e*k_t* (ky * (MY_mes/Dz)) + ko * (MO_mes/Dz))

length((1 - A_a/(e*k_t* (ky * (MY_mes/Dz)) + ko * (MO_mes/Dz))))


Dz=20
k_u_mes=pmax(c(0, (1- A_a/(e*k_t* (ky * (MY_mes/Dz)) + ko * (MO_mes/Dz)))))
k_u_mes= (1- A_a/(e*k_t* (ky * (MY_mes/Dz)) + ko * (MO_mes/Dz)))

length(k_u_mes)
length(MY_mes)

plot(MY_mes, k_u_mes, type="l", col=1)

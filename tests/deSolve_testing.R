
test_pars=list(ky=0.8, ko=0.04,
               kmix=0.02,
               e=0.15,
               Im=0.08, Ir=0.048,
               F_prot=0.1,
               phi_mac=0.04,
               clay=0.2,
               Delta_z_min=14.2,
               gamma_o=1.2,
               gamma_m=2.7,
               proportion=NULL,
               phi_min=0.35,
               f_text_mic=0.5069,
               f_agg=3,
               init=c(My_mes=0.132, Mo_mes=0.507, My_mic=0.284, Mo_mic=0.579),
               sim_length=30)


simulation_deSolve<-run_Porous_deSolve(ky=test_pars$ky, ko=test_pars$ko,
                                       kmix=test_pars$kmix,
                                       e=test_pars$e,
                                       Im=test_pars$Im, Ir=test_pars$Ir,
                                       F_prot=test_pars$F_prot,
                                       phi_mac=test_pars$phi_mac,
                                       clay=test_pars$clay,
                                       Delta_z_min=test_pars$Delta_z_min,
                                       gamma_o=test_pars$gamma_o,
                                       gamma_m=test_pars$gamma_m,
                                       proportion=test_pars$proportion,
                                       phi_min=test_pars$phi_min,
                                       f_text_mic=test_pars$f_text_mic,
                                       f_agg=test_pars$f_agg,
                                       init=c(My_mes=as.numeric(test_pars$init[1]), Mo_mes=as.numeric(test_pars$init[2]), My_mic=as.numeric(test_pars$init[3]), Mo_mic=as.numeric(test_pars$init[4])),
                                       sim_length=test_pars$sim_length,
                                       sim_steps=1)

results_deSolve<-simulation_deSolve$results

range_stocks<-range(results_deSolve[,2:5])
plot(results_deSolve$time, results_deSolve$My_mes.stocks, ylim=c(-2,range_stocks[2]), type="l", col=1)
lines(results_deSolve$time, results_deSolve$Mo_mes.stocks, col=2)
lines(results_deSolve$time, results_deSolve$My_mic.stocks, col=3)
lines(results_deSolve$time, results_deSolve$Mo_mic.stocks, col=4)

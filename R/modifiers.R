#' The main accessory function of the model
#'
#' This function  calculates the proportion of inputs in each of the two youg pools deending on the organic matter content
#' @param phi_mac macroporosity, (\eqn{\frac{cm^3 of water}{cm^3 of soil}})
#' @param phi_min minimal porosity, (\eqn{\frac{cm^3 of water}{cm^3 of soil}})
#' @param Delta_z_min minimal soil thickness if no organic matter was present (cm)
#' @param clay fraction of clay content (dimensionless, between 0 and 1)
#' @param gamma_o density of organic matter (\eqn{g cm^{-3})
#' @param My_mic One of the four model pools (they are all summed up in this function for calculating the total)
#' @param Mo_mic One of the four model pools (they are all summed up in this function for calculating the total)
#' @param My_mes One of the four model pools (they are all summed up in this function for calculating the total)
#' @param Mo_mes One of the four model pools (they are all summed up in this function for calculating the total)
#' @return two values, the proportion of input in the mesopore and micropore Y pools
#' @seealso \code{\link{phi_mic}}, \code{\link{phi_mat}}, \code{\link{f_som}}
#' @export

pore_frac<-function(phi_mac,
                    clay,
                    Delta_z_min,
                    gamma_o,
                    My_mic,
                    Mo_mic,
                    My_mes,
                    Mo_mes,
                    phi_min,
                    f_text_mic=NULL, f_agg){

    phi_mic_calc<-phi_mic(My_mic=My_mic, Mo_mic=Mo_mic, My_mes=My_mes, Mo_mes=Mo_mes, gamma_o=gamma_o, clay=clay, Delta_z_min = Delta_z_min, phi_mac=phi_mac, phi_min = phi_min, f_text_mic = f_text_mic, f_agg = f_agg)

    phi_mat_calc<-phi_mat(My_mic=My_mic, Mo_mic=Mo_mic, My_mes=My_mes, Mo_mes=Mo_mes, gamma_o=gamma_o, Delta_z_min = Delta_z_min, phi_mac=phi_mac, phi_min = phi_min, f_agg = f_agg)

    phi_mes_calc=phi_mat_calc-phi_mic_calc

  mes_f=phi_mes_calc/(phi_mes_calc+phi_mic_calc)
  mic_f=phi_mic_calc/(phi_mes_calc+phi_mic_calc)
  return(c(mes_f, mic_f))
}


#' Proportion of micropores
#'
#' This function calculates the proportion of the textural pore space that comprises micropores. It is used in \code{\link{pore_frac}}.
#' This parameter was intended in the original paper (Meurer et al., 2020) as user defined, but its estimation has been developed further by N. Jarvis (personal communication).
#' The method for its estimation is based on a Brooks-Corey soil water retention model: \deqn{ f_{mic_{text}} = \left( \frac{\psi_{mes  \textbackslash mac}}{\psi_{mic  \textbackslash mes}} \right) ^{\lambda_{mat(t)}}}
#' where \eqn{\psi_{mes  \textbackslash mac}} is the pressure head defining the largest mesopore (set to -0.3) and \eqn{\psi_{mic  \textbackslash mes}} is the pressure head defining the largest micropore (set to -0.6). The parameter
#' \eqn{\lambda_{mat(t)}} is in turn estimated as:
#' \deqn{\lambda_{mat(t)}= \frac{log\left(\frac{\theta_w}{\phi_{min}}\right)}  {log\left(\frac{\psi_{mes  \textbackslash mac}}{\psi_{w}}\right)}}
#' where \eqn{psi_w} is the wilting oint pressure head (set to -150 m) and \eqn{\theta_w} is estimated from a pedotransfer function:
#' \deqn{\theta_w=0.004+0.5 \cdot f_{clay}}
#' where \eqn{f_{clay}} is the soil clay content (\ifelse{html}{\out{kg kg<sup>-1</sup>}}{\eqn{kg \: kg^{-1}}}   ).
#' @inheritParams pore_frac
#' @return one single value
#' @export

f_text_mic_func<-function(clay, phi_min){

  psi_w=-150 #meters, the wilting point pressure head

  theta_w=0.004+0.5*clay #psi_w (m) is the wilting point pressure head (= -150 m) and the corresponding water content theta_w (m3 m-3) is estimated from a pedotransfer function (Ostovari et al., 2015)

  lamda_mat=log(theta_w/phi_min)/log(-0.3/psi_w) #lambda_mat is the pore size distribution index

  f_text_mic=(-0.3/-6)^lamda_mat

  return(f_text_mic)
}


#' Variation of the thickness of soil layer
#'
#' @description This function calculates the variation of the thickness of soil layer as a function of organic matter. The parameter \eqn{f_{agg}} should be estimated from data on the relationship between bulk density (or its inverse, the specific volume) and soil organic matter content (see eq. 19 and fig. 4 in Meurer et al., 2020; from this data and other studies, a good average value of fagg should be around 3, which is the default value)
#' @param f_agg an aggregation factor (m3 pore space m-3 organic matter) defined as the slope of the linear relationship assumed between the volume of aggregation pore space \eqn{V_{agg}}, and the volume of organic matter \eqn{V_{s_o}} (dimensionless)
#' @inheritParams pore_frac
#' @return one single value
#' @export

Delta_z<-function(f_agg=f_agg,
                  Delta_z_min,
                  My_mic, Mo_mic, My_mes, Mo_mes,
                  phi_mac,
                  gamma_o){

  Mso=My_mic + Mo_mic + My_mes + Mo_mes
  Delta_z=(((1+f_agg)*(Mso/gamma_o))+Delta_z_min)/(1-phi_mac)

  return(Delta_z)
}



#' Microporosity
#'
#' @description This function calculates the microporosity \eqn{\phi_{mic}} based on the variation of organic matter in the soil. It is used in \code{\link{pore_frac}}
#' @inheritParams pore_frac
#' @inheritParams Delta_z
#' @inheritParams f_text_mic_func
#' @f_text_mic in case the user wants to override the calculation of f_text_mic, a value can be specified here. Otherwise the value is calculated with the function \code{\link{f_text_mic_func}}
#' @return one single value
#' @export

phi_mic<-function(My_mic, Mo_mic, My_mes, Mo_mes,
                  gamma_o, #density of organic matter
                  f_agg=f_agg, #default suggested by Nick Jarvis
                  clay,
                  Delta_z_min, #soil layer thickess (m)
                  phi_min, #minimum matrix porosity, STILL MISSING
                  phi_mac,
                  f_text_mic){
  if(is.null(f_text_mic)){
    f_text_mic_calc=f_text_mic_func(clay=clay, phi_min=phi_min)
  } else {
    f_text_mic_calc=f_text_mic
  }


  Delta_z_calc=Delta_z(f_agg=f_agg,
                       Delta_z_min=Delta_z_min,
                       My_mic, Mo_mic, My_mes, Mo_mes,
                       phi_mac,
                       gamma_o)

  phi_mic=((f_agg*((My_mic+Mo_mic)/gamma_o))+f_text_mic_calc*Delta_z_min*phi_min)/Delta_z_calc
  return(phi_mic)
}




#' Matrix porosity
#'
#' @description This function calculates the matrix porosity \eqn{\phi_{mac}} based on the variation of organic matter in the soil. It is used in \code{\link{pore_frac}} to calculate mesoporosity \eqn{\phi_{mes}=\phi_{mat}-\phi_{mic}}
#' @inheritParams pore_frac
#' @inheritParams Delta_z
#' @return one single value
#' @export
phi_mat<-function(My_mic, Mo_mic, My_mes, Mo_mes,
                  gamma_o,
                  f_agg,
                  Delta_z_min,#soil layer thickess (m)
                  phi_min,
                  phi_mac
                  ){

  Mso=My_mic + Mo_mic + My_mes + Mo_mes

  Delta_z_calc=Delta_z(f_agg,
                       Delta_z_min = Delta_z_min,
                       My_mic, Mo_mic, My_mes, Mo_mes,
                       phi_mac,
                       gamma_o)


  phi_mat=(f_agg*(Mso/gamma_o)+Delta_z_min*phi_min)/Delta_z_calc
  return(phi_mat)
}



#' Mineral mass
#'
#' @description This function calculates the mineral mass
#' @param gamma_m mineral matter density (\eqn{g cm^{-3})
#' @inheritParams Delta_z
#' @return one single value
#' @seealso \code{\link{Delta_z_min}}, \code{\link{gamma_b}}, \code{\link{f_som}}
#' @export
Msm<-function(Delta_z_min,#soil layer thickess (m)
                  phi_min,
                  gamma_m
){

 Msm=Delta_z_min*gamma_m*(1-phi_min)

   return(Msm)
}




#' Soil bulk density
#'
#' @description This function calculates the soil bulk density. It relies on \code{\link{Msm}} to calculate the mineral mass.
#' @inheritParams Msm
#' @inheritParams pore_frac
#' @inheritParams Delta_z
#' @return one single value
#' @seealso \code{\link{Msm}}
#' @export
gamma_b<-function(My_mic, Mo_mic, My_mes, Mo_mes,
                  Delta_z_min,
                  phi_min,
                  gamma_o,
                  gamma_m,
                  f_agg,
                  phi_mac
){

  Delta_z_calc=Delta_z(f_agg=f_agg,
                       Delta_z_min=Delta_z_min,
                       My_mic, Mo_mic, My_mes, Mo_mes,
                       phi_mac,
                       gamma_o)

  Mso_calc = My_mic + Mo_mic + My_mes + Mo_mes
  Msm_calc = Msm(Delta_z_min,
                 phi_min,
                 gamma_m)

  gamma_b= (Mso_calc + Msm_calc)/Delta_z_calc

  return(gamma_b)
}




#' Soil C concentration
#'
#' @description This function calculates the soil C concentration. It relies on \code{\link{Msm}} to calculate the mineral mass.
#' @inheritParams Msm
#' @inheritParams pore_frac
#' @inheritParams Delta_z
#' @return one single value
#' @seealso \code{\link{Msm}}
#' @export
f_som<-function(My_mic, Mo_mic, My_mes, Mo_mes,
                  Delta_z_min,
                  phi_min,
                  gamma_m
                ){

  Mso_calc = My_mic + Mo_mic + My_mes + Mo_mes
  Msm_calc = Msm(Delta_z_min,
                 phi_min,
                 gamma_m)

  f_som = Mso_calc / (Mso_calc+Msm_calc)

  return(f_som)
}







#' Run Porous
#'
#' @description This function is a wrapper for running the main model \code{\link{Porous}}, which is a wrapper for the dual porosity decomposition model (equations 1 to 6) described in Meurer et al. (2020) implemented in  `SoilR`.
#' It then feeds the simulated SOC stocks to the functions \code{\link{gamma_b}} and \code{\link{f_som}} to calculate the variation of bulk density and C concentrations
#' @inheritParams Porous
#' @param init initialization of the four pools. This must be a vector with the names of the four pools in sequence in the format: `c(My_mes=1, Mo_mes=10,My_mic=0.6, Mo_mic=3)`. Units are generally in  (\eqn{g cm^{-2} year^{-1}), but in any case they should match the units of the inputs.
#' @param sim_lenght length of the simulation (years)
#' @param sim_steps steps of the simulation (fraction of one year)
#' @return a data frame with the simulation of C stocks for each of the four pools, soil respiration for each of the four pools, bulk density variation and SOC concentration.
#' @seealso \code{\link{Porous}}, \code{\link{gamma_b}}, \code{\link{f_som}}
#' @export
run_Porous<-function(ky=0.8,
                     ko=0.04,
                     kmix=0.02,
                     e=0.15,
                     Im=0.08,
                     Ir=0.048,
                     F_prot=0.1,
                     phi_mac=0.04,
                     phi_min=0.35,
                     clay=0.2,
                     Delta_z_min=14.2,
                     gamma_o=1.2,
                     gamma_m=2.7,
                     proportion=NULL,
                     f_text_mic=NULL,
                     f_agg=3,
                     init,
                     sim_length,
                     sim_steps
){

  modelObject<-Porous(ky, ko,
                      kmix,
                      e,
                      Im, Ir,
                      F_prot,
                      phi_mac,
                      phi_min,
                      clay,
                      Delta_z_min,
                      gamma_o)



  times<-seq(0,sim_length,by=sim_steps)

  modrun0<-Model_by_PoolNames(smod=modelObject, times=times, initialValues=init)

  Stocks<-getC(modrun0)
  Resp<-getReleaseFlux(modrun0)
  colnames(Stocks)=c("My_mes stocks", "Mo_mes stocks", "My_mic stocks", "Mo_mic stocks")
  colnames(Resp)=c("My_mes resp", "Mo_mes resp", "My_mic resp", "Mo_mic resp")

  f_som_sim<-f_som(My_mic=Stocks[,1],
        Mo_mic=Stocks[,2],
        My_mes=Stocks[,3],
        Mo_mes=Stocks[,4],
        Delta_z_min,
        phi_min,
        gamma_m)

  gamma_b_sim<-gamma_b(My_mic=Stocks[,1],
          Mo_mic=Stocks[,2],
          My_mes=Stocks[,3],
          Mo_mes=Stocks[,4],
          Delta_z_min,
          phi_min,
          gamma_o,
          gamma_m,
          f_agg,
          phi_mac)

  Delta_z_sim<-Delta_z(f_agg=f_agg,
               Delta_z_min=Delta_z_min,
               My_mic=Stocks[,1],
               Mo_mic=Stocks[,2],
               My_mes=Stocks[,3],
               Mo_mes=Stocks[,4],
               phi_mac,
               gamma_o)

  Stocks[,"My_mes stocks"]
  # results<-data.frame(Stocks[,"My_mes stocks"],
  #                       Stocks[,"Mo_mes stocks"],
  #                       Stocks[,"My_mic stocks"],
  #                       Stocks[,"Mo_mic stocks"],
  #                       Resp[,"My_mes resp"],
  #                       Resp[,"Mo_mes resp"],
  #                       Resp[,"My_mic resp"],
  #                       Resp[,"Mo_mic resp"], f_som_sim, gamma_b_sim) #create the data frame with the results
  results<-data.frame(time=times, Stocks,
                      Resp, f_som_sim, gamma_b_sim, Delta_z_sim) #create the data frame with the results

    return(results)

}




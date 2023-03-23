#' The main accessory function of the model
#'
#' This function  calculates the proportion of inputs in each of the two youg pools deending on the organic matter content
#' @param phi_mac macroporosity
#' @param Delta_z_min minimal soil thickness if no organic matter was present
#' @param clay fraction of clay content
#' @param gamma_o density of organic matter
#' @param My_mic One of the four model pools (they are all summed up in this function for calculating the total)
#' @param Mo_mic One of the four model pools (they are all summed up in this function for calculating the total)
#' @param My_mes One of the four model pools (they are all summed up in this function for calculating the total)
#' @param Mo_mes One of the four model pools (they are all summed up in this function for calculating the total)
#' @return two values, the proportion of input in the mesopore and micropore Y pools
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
                    f_text_mic=NULL){

    phi_mic_calc<-phi_mic(My_mic=My_mic, Mo_mic=Mo_mic, My_mes=My_mes, Mo_mes=Mo_mes, gamma_o=gamma_o, clay=clay, Delta_z_min = Delta_z_min, phi_mac=phi_mac, phi_min = phi_min, f_text_mic = f_text_mic)

    phi_mat_calc<-phi_mat(My_mic=My_mic, Mo_mic=Mo_mic, My_mes=My_mes, Mo_mes=Mo_mes, gamma_o=gamma_o, Delta_z_min = Delta_z_min, phi_mac=phi_mac, phi_min = phi_min)

    phi_mes_calc=phi_mat_calc-phi_mic_calc

  mes_f=phi_mes_calc/(phi_mes_calc+phi_mic_calc)
  mic_f=phi_mic_calc/(phi_mes_calc+phi_mic_calc)
  return(c(mes_f, mic_f))
}


#' Proportion of micropores
#'
#' This function calculates the proportion of the textural pore space that comprises micropores. It is used in \code{\link{pore_frac}}
#' @param clay soil clay fraction
#' @param phi_min minimal porosity, user defined
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
#' This function calculates the variation of the thickness of soil layer as a function of organic matter
#' @param f_agg an aggregation factor (m3 pore space m-3 organic matter) defined as the slope of the linear relationship assumed between the volume of aggregation pore space \eqn{V_{agg}}, and the volume of organic matter \eqn{V_{s_o}}
#' @param Delta_z_min minimal soil thickness if no organic matter was present
#' @inheritParams pore_frac
#' @param gamma_o density of organic matter
#' @return one single value
#' @export

Delta_z<-function(f_agg=3,
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
#' This function calculates the microporosity \eqn{\phi_{mic}} based on the variation of organic matter in the soil. It is used in \code{\link{pore_frac}}
#' @inheritParams pore_frac
#' @inheritParams Delta_z
#' @inheritParams f_text_mic_func
#' @f_text_mic in case the user wants to override the calculation of f_text_mic, a value can be specified here. Otherwise the value is calculated with the function \code{\link{f_text_mic_func}}
#' @return one single value
#' @export

phi_mic<-function(My_mic, Mo_mic, My_mes, Mo_mes,
                  gamma_o, #density of organic matter
                  f_agg=3, #default suggested by Nick Jarvis
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


  Delta_z_calc=Delta_z(f_agg=3,
                       Delta_z_min=Delta_z_min,
                       My_mic, Mo_mic, My_mes, Mo_mes,
                       phi_mac,
                       gamma_o)

  phi_mic=((f_agg*((My_mic+Mo_mic)/gamma_o))+f_text_mic_calc*Delta_z_min*phi_min)/Delta_z_calc
  return(phi_mic)
}




#' Matrix porosity
#'
#' This function calculates the matrix porosity \eqn{\phi_{mac}} based on the variation of organic matter in the soil. It is used in \code{\link{pore_frac}} to calculate mesoporosity \eqn{\phi_{mes}=\phi_{mat}-\phi_{mic}}
#' @inheritParams pore_frac
#' @inheritParams Delta_z
#' @return one single value
#' @export
phi_mat<-function(My_mic, Mo_mic, My_mes, Mo_mes,
                  gamma_o,
                  f_agg=3,
                  Delta_z_min,#soil layer thickess (m)
                  phi_min,
                  phi_mac
                  ){

  Mso=My_mic + Mo_mic + My_mes + Mo_mes

  Delta_z_calc=Delta_z(f_agg=3,
                       Delta_z_min = Delta_z_min,
                       My_mic, Mo_mic, My_mes, Mo_mes,
                       phi_mac,
                       gamma_o)


  phi_mat=(f_agg*(Mso/gamma_o)+Delta_z_min*phi_min)/Delta_z_calc
  return(phi_mat)
}

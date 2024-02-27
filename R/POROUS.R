#' The SOC decomposition model **STILL TESTING; NOT YET WORKING, use the deSolve function instead**
#'
#' This function implements with the SoilR model development framework the dual porosity model described in Meurer et al. (2020). THis function provides just the model definition, that needs then to be run with specific SoilR functions.
#' A more comfortable wrapper for doing this automatically is provided in the package (\code{\link{run_Porous}}).
#' The model is an evolution of a two-pool linear SOC model, with two pools (young and old material) running in parallel for micro- and mesopores.
#' While aboveground inputs are rooted in the mesopores, root inputs are distributed between micro and mesopores depending on porosity, which is in turn influenced by organic matter. This makes the model nonlinear, although it still behaves similarly to a linear model within a reasonable calibration range. The model is described by a series of four equations: \cr
#' \cr
#' \eqn{\frac{dM_{Y_{(mes)}}}{dt} = I_m + \left( \frac{\phi_{mes}}{\phi_{mes}+\phi_{mic}}\right) \cdot I_r - k_Y \cdot M_{Y_{(mes)}}+ T_Y } \cr
#' \eqn{ \frac{dM_{O_{(mes)}}}{dt} = \left( \epsilon \cdot k_Y \cdot M_{Y_{(mes)}} \right) - \left( (1- \epsilon) \cdot k_O \cdot M_{O_{(mes)}} \right) + T_O} \cr
#' \eqn{ \frac{dM_{Y_{(mic)}}}{dt} = \left( \frac{\phi_{mic}}{\phi_{mes}+\phi_{mic}}\right) \cdot I_r - k_Y \cdot F_{prot} \cdot M_{Y_{(mes)}}- T_Y } \cr
#' \eqn{ \frac{dM_{O_{(mic)}}}{dt} = \left( \epsilon \cdot k_Y \cdot F_{prot} \cdot M_{Y_{(mes)}} \right) - \left( (1- \epsilon) \cdot k_O \cdot F_{prot} \cdot M_{O_{(mes)}} \right) - T_O } \cr
#'  \cr
#' Please refer to the original paper for more details.  \cr
#' The terms \eqn{T_Y} and \eqn{T_Y} are calculated in the original paper as  \eqn{k_{mix} \cdot \frac{(My_{mic}-My_{mes})}{2}} and \eqn{k_{mix} \cdot \frac{(Mo_{mic}-Mo_{mes})}{2}}, but in a later model development the term 2 disappeared and are now calculaged as
#' \eqn{T_Y = k_{mix} \cdot (My_{mic}-My_{mes})} and \eqn{T_O = k_{mix} \cdot (Mo_{mic}-Mo_{mes})}. \cr
#' \cr
#' The two porosity terms, \eqn{\phi_{mes} = f(M_{Y_{(mes)}}, M_{O_{(mes)}},M_{Y_{(mic)}}, M_{O_{(mic)}})} and \eqn{\phi_{mic} = f(M_{Y_{(mic)}}, M_{O_{(mic)}})}, are dependent on the variation of the different C pools and everything is variable over time, introducing a nonlinearity in the system and defining the biggest peculiarity of this model.  \cr
#' After substituting the terms \eqn{\left( \frac{\phi_{mes}(t)}{\phi_{mes}(t)+\phi_{mic}(t)}\right) = \varphi_{mes}} and \eqn{\left( \frac{\phi_{mic}(t)}{\phi_{mes}(t)+\phi_{mic}(t)}\right) = \varphi_{mic}},
#' The model can be rewritten in matrix form as :  \cr
#' \cr
#' \eqn{I_m(t) + I_r(t) \cdot N(C,t) + A(t) \cdot P(t) \cdot C(t)} \cr
# THE FOLLOWING EQUATION IS COMMENTED BECAUSE IT GAVE ERROR...
# Or, more explicitly:  \cr
#  \deqn{
#  \frac{dC}{dt}=\begin{bmatrix} I_m\\  0\\  0\\ 0\\ \end{bmatrix} +
#    \begin{bmatrix} I_r\\ 0\\ I_r\\ 0\\ \end{bmatrix} \cdot
#    \begin{bmatrix} \varphi_{mes} & 0 & 0 & 0\\ 0 & 1 & 0 & 0\\ 0 & 0 & \varphi_{mic} & 0\\ 0 & 0 & 0 & 1\\ \end{bmatrix}+
#    \begin{bmatrix} -k_y & \epsilon & 0 & 0\\ 0 & -k_o & 0 & 0\\ T_Y & 0 & -k_y & \epsilon\\ 0 & T_O & 0 & -k_o\\ \end{bmatrix} \cdot
#    \begin{bmatrix} 1 & 0 & 0 & 0\\ 0 & 1 & 0 & 0\\ 0 & 0 & F_{prot} & 0\\ 0 & 0 & 0 & F_{prot}\\ \end{bmatrix} \cdot
#    \begin{bmatrix} M_{Y_{mes}}\\ M_{O_{mes}}\\ M_{Y_{mic}}\\ M_{O_{mic}} \end{bmatrix}
#         }
#' \cr
#' The model is implemented with the \code{SoilR} package but it is relying on a more conventional ODE definition (not its matrix form).
#'
#' @param ky decomposition constant of the Young pool \eqn{frac{1}{year}}
#' @param ko decomposition constant of the Old pool \eqn{frac{1}{year}}
#' @param kmix mixing rate \eqn{frac{1}{year}}
#' @param e efficiency, which is the transfer term between the pools and corresponds to the term h in the ICBM model in Kätterer et al. (2001) (dimensionless)
#' @param Im Inputs from aboveground.  Units are generally in (\eqn{g cm^{-2} year^{-1}), but in any case they should match the units of the initialization.
#' @param Ir Inputs from roots, which is partitioned between micropore and mesopores with the function \code{\link{pore_frac}}. Units are generally in (\eqn{g cm^{-2} year^{-1}), but in any case they should match the units of the initialization.
#' @param F_prot protection provided by the micropore space (dimensionless)
#' @param proportion this is a linearization term to make the proportion of the inputs between micro- and mesopores constant. If NULL (or not specified, since default is NULL) then the model is running as nonlinear, as in the original paper. If specified (must be between 0 and 1) then the model is linearized adopting this value as fixed proportion of inputs from roots going into the mesopore space (and its reciprocal into the micropore)
#' @inheritParams pore_frac
#' @inheritParams Delta_z
#' @inheritParams f_text_mic_func
#' @f_text_mic if this value is user defined, the function \code{\link{f_text_mic_func}} is overridden, otherwise it is used to calculate it.
#' @return two values, the proportion (between 0 and 1) of input in the mesopore and micropore Y pools
#' @seealso \code{\link{run_Porous}}
#' @references Meurer, Katharina Hildegard Elisabeth, Claire Chenu, Elsa Coucheney, Anke Marianne Herrmann, Thomas Keller, Thomas Kätterer, David Nimblad Svensson, and Nicholas Jarvis. “Modelling Dynamic Interactions between Soil Structure and the Storage and Turnover of Soil Organic Matter.” Biogeosciences 17, no. 20 (October 19, 2020): 5025–42. https://doi.org/10.5194/bg-17-5025-2020. \cr
#' Kätterer, Thomas, and Olof Andrén. “The ICBM Family of Analytically Solved Models of Soil Carbon, Nitrogen and Microbial Biomass Dynamics — Descriptions and Application Examples.” Ecological Modelling 136, no. 2–3 (January 2001): 191–207. https://doi.org/10.1016/S0304-3800(00)00420-8.
#' @export
#'
Porous<-function(ky,
                 ko,
                 kmix,
                 e,
                 Im,
                 Ir,
                 F_prot,
                 phi_mac,
                 phi_min,
                 clay,
                 Delta_z_min,
                 gamma_o,
                 proportion=NULL,
                 f_text_mic=NULL,
                 f_agg){

  time_symbol='t'



  ##### IN
  ifs=SoilR:::InFluxList_by_PoolName(
    c(
      # #Inputs from shoots, direcly in Y mesopores
      # SoilR:::InFlux_by_PoolName(
      #   destinationName='My_mes',
      #   func=function(t){
      #     Im
      #   }
      # ),
      #Inputs from roots, partitioned based on either a pfixed partitioning factor or a calculated one
      # Inputs to Y mesopores from roots
      if(is.null(proportion)){ #if cycle for inputs going to My_mes, if linear or not. If "proportion" is missing then nonlinear
        SoilR:::InFlux_by_PoolName(
          destinationName='My_mes',
          func=function(t, My_mes, Mo_mes, My_mic, Mo_mic){
            Im+(Ir*pore_frac(phi_mac, clay, Delta_z_min, gamma_o, My_mes, Mo_mes, My_mic, Mo_mic, phi_min, f_text_mic, f_agg)[1])
          }
        )
      } else{ #... else use the proportion
        SoilR:::InFlux_by_PoolName(
          destinationName='My_mes',
          func=function(t, My_mes, Mo_mes, My_mic, Mo_mic){
            Im+(Ir*proportion)
          }
      )
      },
      # Inputs to Y micropores from roots
      if(is.null(proportion)){#if cycle for inputs going to My_mic, if linear or not. If "proportion" is missing then nonlinear
      SoilR:::InFlux_by_PoolName(
        destinationName='My_mic',
        func=function(t, My_mes, Mo_mes, My_mic, Mo_mic){
          Ir*pore_frac(phi_mac, clay, Delta_z_min, gamma_o, My_mes, Mo_mes, My_mic, Mo_mic, phi_min, f_text_mic, f_agg)[2]
        }
      )
    } else {#... else use the reciprocal of the proportion
      SoilR:::InFlux_by_PoolName(
        destinationName='My_mic',
        func=function(t, My_mes, Mo_mes, My_mic, Mo_mic){
          Ir*(1-proportion)
        }
      )
    }
  )
)
  ##### OUT
  ofs=SoilR:::OutFluxList_by_PoolName(
    c(
      SoilR:::OutFlux_by_PoolName(
        sourceName='My_mes',
        func=function(My_mes){
          (1-e)*ky*My_mes
        }
      )
      ,
      SoilR:::OutFlux_by_PoolName(
        sourceName='Mo_mes',
        func=function(Mo_mes){
         (1-e)*ko*Mo_mes
        }
      ),
      SoilR:::OutFlux_by_PoolName(
        sourceName='My_mic',
        func=function(My_mic){
          (1-e)*(ky*F_prot*My_mic)
        }
      ),
      SoilR:::OutFlux_by_PoolName(
        sourceName='Mo_mic',
        func=function(Mo_mic){
          ((1-e)*ko*F_prot*Mo_mic)
        }
      )
    )
  )

  ##### INT
  intfs=SoilR:::InternalFluxList_by_PoolName(
    list(
      #humification fluxes
      SoilR:::InternalFlux_by_PoolName(
        sourceName='My_mes',
        destinationName='Mo_mes',
        func=function(My_mes){
          e*ky*My_mes
        }
      ),
      SoilR:::InternalFlux_by_PoolName(
        sourceName='My_mic',
        destinationName='Mo_mic',
        func=function(My_mic){
          e*ky*F_prot*My_mic
        }
      )
      ,
      #bioturbation fluxes
      #Ty
      SoilR:::InternalFlux_by_PoolName(
       sourceName='My_mic',
       destinationName='My_mes',
       func=function(My_mic, My_mes){
         kmix*(My_mic-My_mes)
       }
      ),
      #To
      SoilR:::InternalFlux_by_PoolName(
       sourceName='Mo_mic',
       destinationName='Mo_mes',
       func=function(Mo_mic, Mo_mes){
         kmix*(Mo_mic-Mo_mes)
       }
      )
    )
  )

  smod <- SoilR:::SymbolicModel_by_PoolNames(
    in_fluxes=ifs,
    internal_fluxes=intfs,
    out_fluxes=ofs,
    timeSymbol=time_symbol
  )
  smod
}



#' Run Porous **STILL TESTING; NOT YET WORKING, use the deSolve function instead**
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
run_Porous<-function(ky,
                     ko,
                     kmix,
                     e,
                     Im,
                     Ir,
                     F_prot,
                     phi_mac,
                     phi_min,
                     clay,
                     Delta_z_min,
                     gamma_o,
                     gamma_m,
                     proportion=NULL,
                     f_text_mic=NULL,
                     f_agg,
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
                      gamma_o,
                      f_agg,
                      proportion)



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

  ### Mass balance check
  #SOC difference
  SOC_diff<-rowSums(Stocks)[length(times)]-rowSums(Stocks)[1]

  #Respiration (cumulated)
  RESP_tot<-sum(rowSums(Resp*sim_steps)[-1]) #removing the first respiration

  #Inputs total
  Input_total<-sum(rep(Im, sim_length)+ rep(Ir, sim_length))
  # check mass balance
  Simulated_total=(SOC_diff+RESP_tot)

  mass_balance=c(Input_total, Simulated_total)

  convergence=all.equal(mass_balance[1], mass_balance[2])


  results<-data.frame(time=times, Stocks,
                      Resp, f_som_sim, gamma_b_sim, Delta_z_sim) #create the data frame with the results


  return(list(results=results, "mass balance"=mass_balance, "mass balance convergence"=convergence))

}






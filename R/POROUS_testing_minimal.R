#' Purely for testing purposes
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
Porous_test_min<-function(ky,
                 ko,
                 e,
                 Im,
                 Ir,
                 F_prot,
                 proportion){

  time_symbol='t'



  ##### IN
  ifs=SoilR:::InFluxList_by_PoolName(
    c(
        SoilR:::InFlux_by_PoolName(
          destinationName='My_mes',
          func=function(t, My_mes, Mo_mes, My_mic, Mo_mic){
            Im+(Ir*proportion)
          }
        ),
        SoilR:::InFlux_by_PoolName(
          destinationName='My_mic',
          func=function(t, My_mes, Mo_mes, My_mic, Mo_mic){
            Ir*(1-proportion)
          }
       )
    )
  )
  ##### OUT
  ofs=SoilR:::OutFluxList_by_PoolName(
    c(
      SoilR:::OutFlux_by_PoolName(
        sourceName='My_mes',
        func=function(My_mes){
          ky*My_mes
        }
      )
      ,
      SoilR:::OutFlux_by_PoolName(
        sourceName='Mo_mes',
        func=function(Mo_mes){
          ko*Mo_mes
        }
      ),
      SoilR:::OutFlux_by_PoolName(
        sourceName='My_mic',
        func=function(My_mic){
          (ky*My_mic)
        }
      ),
      SoilR:::OutFlux_by_PoolName(
        sourceName='Mo_mic',
        func=function(Mo_mic){
          (ko*Mo_mic)
        }
      )
    )
  )

  ##### INT
  intfs=SoilR:::InternalFluxList_by_PoolName(
    list(
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
          e*ky*My_mic
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




run_Porous_test_min<-function(ky,
                     ko,
                     e,
                     Im,
                     Ir,
                     F_prot,
                     proportion,
                     init,
                     sim_length,
                     sim_steps
){

  modelObject_test_min<-Porous_test_min(ky, ko,
                      e,
                      proportion,
                      F_prot,
                      Im, Ir)



  times<-seq(0,sim_length,by=sim_steps)

  modrun0_test<-Model_by_PoolNames(smod=modelObject_test_min, times=times, initialValues=init)

  Stocks<-getC(modrun0_test)
  Resp<-getReleaseFlux(modrun0_test)
  colnames(Stocks)=c("My_mes stocks", "Mo_mes stocks", "My_mic stocks", "Mo_mic stocks")
  colnames(Resp)=c("My_mes resp", "Mo_mes resp", "My_mic resp", "Mo_mic resp")

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
                      Resp) #create the data frame with the results


  return(list(results=results, "mass balance"=mass_balance, "mass balance convergence"=convergence))

}



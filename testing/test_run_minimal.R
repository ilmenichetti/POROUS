

library(SoilR)
#library(Porous)
#library(devtools)
#devtools::install_github("ilmenichetti/Porous")

#building the manual as a site
#pkgdown::build_site()

#building the manual in pdf
#devtools::build_manual(, path = "../Porous")

input_m=5
input_r=2

Porous_test<-function(ky,
                      ko,
                      e,
                      Im,
                      Ir,
                      proportion){

  time_symbol='t'
  ##### IN
  ifs=SoilR:::InFluxList_by_PoolName(
    c(
      #Inputs from shoots, direcly in Y mesopores
      SoilR:::InFlux_by_PoolName(
        destinationName='My_mes',
        func=function(t){
          Im+(Ir*proportion)
        }
      )
      ,
      # #Inputs from roots, partitioned based on either a fixed partitioning factor or a calculated one
      # # Inputs to Y mesopores from roots
      # SoilR:::InFlux_by_PoolName(
      #   destinationName='My_mes',
      #   func=function(t){
      #     Ir*proportion
      #   }
      # )
      # ,
      # Inputs to Y micropores from roots

      SoilR:::InFlux_by_PoolName(
        destinationName='My_mic',
        func=function(t){
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
          ky*My_mic
        }
      ),
      SoilR:::OutFlux_by_PoolName(
        sourceName='Mo_mic',
        func=function(Mo_mic){
          ko*Mo_mic
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






modelObject<-Porous_test(ky=0.8, ko=0.00605,
                         e=0.15,
                         Im=input_m, Ir=input_r,
                         proportion=1)
plotPoolGraph(modelObject)

iv<-c(My_mes=1.32, Mo_mes=5.07, My_mic=2.84, Mo_mic=5.79)

time_step=1
duration=20

times<-seq(0,duration,by=time_step)
modrun0<-Model_by_PoolNames(smod=modelObject, times=times, initialValues=iv)

Ct0<-getC(modrun0)
Rt0<-getReleaseFlux(modrun0)




# checking the results (mass balance)
### mass balance

#SOC difference
SOC_diff<-rowSums(Ct0)[length(times)]-rowSums(Ct0)[1]

#Respiration (cumulated)
RESP_tot<-sum(rowSums(Rt0*time_step)[-1])

#Inputs total
Input_tot<-sum(rep(input_m+input_r, duration))
# check mass balance
Input_tot==(SOC_diff+RESP_tot)
Input_tot
(SOC_diff+RESP_tot)


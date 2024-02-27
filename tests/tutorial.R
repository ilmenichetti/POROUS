library(SoilR)

Iy_value=1.6

TwopoolNonlinearInput<-function(ky=0.8, ko=0.00605, h=0.13, Iy=Iy_value, proportion){
  time_symbol='t'

  ifs=SoilR:::InFluxList_by_PoolName(
    c(
      SoilR:::InFlux_by_PoolName(
        destinationName='Cy',
        func=function(t){
          Iy*proportion
        }
      ),
      SoilR:::InFlux_by_PoolName(
        destinationName='Co',
        func=function(t){
          Iy*(1-proportion)
        }
      )
    )
  )
  ofs=SoilR:::OutFluxList_by_PoolName(
    c(
      SoilR:::OutFlux_by_PoolName(
        sourceName='Cy',
        func=function(Cy){
          ky*Cy
        }
      )
      ,
      SoilR:::OutFlux_by_PoolName(
        sourceName='Co',
        func=function(Co){
          ko*Co
        }
      )
    )
  )
  intfs=SoilR:::InternalFluxList_by_PoolName(
    list(
      SoilR:::InternalFlux_by_PoolName(
        sourceName='Cy',
        destinationName='Co',
        func=function(Cy){
          h*ky*Cy
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

modelObject<-TwopoolNonlinearInput(proportion=0.3)
plotPoolGraph(modelObject)

iv<-c(Cy=1, Co=10)
duration=20
time_step=0.1

times<-seq(0,duration,by=time_step)
modrun0<-Model_by_PoolNames(smod=modelObject, times=times, initialValues=iv)
Ct0<-getC(modrun0)
Rt0<-getReleaseFlux(modrun0)

plot(rowSums(Ct0), type="l")
plot(rowSums(Rt0), type="l")


### mass balance ###

#SOC difference
SOC_diff<-rowSums(Ct0)[length(times)]-rowSums(Ct0)[1]

#Respiration (cumulated)
plot(rowSums(Rt0))
RESP_tot<-sum(rowSums(Rt0*time_step))

#Inputs total
Input_tot<-sum(rep(Iy_value, duration))

# check mass balance
Input_tot==(SOC_diff+RESP_tot)

Input_tot
(SOC_diff+RESP_tot)



#### ICBM
iv<-c(Cy=1, Co=10)
duration=20
time_step=0.1

times<-seq(0,duration,by=time_step)
ICBM_run=ICBMModel(t=times, h=0.250, r=1.10, c0=iv, In=Iy_value+Io_value) #Manure
Ct0<-getC(ICBM_run)
Rt0<-getReleaseFlux(ICBM_run)

plot(rowSums(Ct0), type="l")

colnames(Ct0)=c("Y stocks", "O stocks")
colnames(Rt0)=c("Y resp", "O resp")

### mass balance ###

#SOC difference
SOC_diff<-rowSums(Ct0)[length(times)]-rowSums(Ct0)[1]

#Respiration (cumulated)
plot(rowSums(Rt0))
RESP_tot<-sum(rowSums(Rt0*time_step)[-1])

#Inputs total
Input_tot<-sum(rep(Iy_value, duration))

# check mass balance
Input_tot==(SOC_diff+RESP_tot)
Input_tot
(SOC_diff+RESP_tot)


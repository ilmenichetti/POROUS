## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(Porous)
library(SoilR) #this should load automatically anyway since it is required by Porous

## -----------------------------------------------------------------------------
 modelObject<-Porous()

## -----------------------------------------------------------------------------
modelObject<-Porous(ky=0.8, ko=0.00605,
                    kmix=0.9,
                    e=0.13,
                    Im=1.1, Ir=0.5,
                    F_prot=0.0,
                    phi_mac=0.2,
                    clay=0.2,
                    Delta_z_min=20,
                    gamma_o=1.2)


## -----------------------------------------------------------------------------
plotPoolGraph(modelObject)

## -----------------------------------------------------------------------------
init<-c(My_mes=1, Mo_mes=10,My_mic=0.6, Mo_mic=3)
times<-seq(0,20,by=0.1)

## -----------------------------------------------------------------------------
modrun0<-Model_by_PoolNames(smod=modelObject, times=times, initialValues=init)

## -----------------------------------------------------------------------------
Stocks<-getC(modrun0)
Resp<-getReleaseFlux(modrun0)

## -----------------------------------------------------------------------------
head(Resp)

## ----fig1, fig.height = 8, fig.width = 6--------------------------------------
par(mfrow=c(2,1), mar=c(4,4,0,1))
matplot(times, Stocks, type="l", lty=1, col=1:4, xlab=" ", ylab="Pool contents", bty="n")
legend("bottomright", c("My_mes", "Mo_mes", "My_mic", "Mo_mic"), lty=1, col=1:4, bty="n")
matplot(times, Resp,  type="l", lty=1, col=1:2, xlab="Time", ylab="Respiration", bty="n")


#
# modelObject<-Porous()
# plotPoolGraph(modelObject)
#
# iv<-c(My_mes=1, Mo_mes=10,My_mic=0.6, Mo_mic=3)
# times<-seq(0,20,by=0.1)
# modrun0<-Model_by_PoolNames(smod=modelObject, times=times, initialValues=iv)
#
# Ct0<-getC(modrun0)
# Rt0<-getReleaseFlux(modrun0)
#
#
# #par(mfrow=c(2,1), mar=c(4,4,0,1))
# matplot(times, Ct0, type="l", lty=1, col=1:4, xlab=" ", ylab="Pool contents", bty="n")
# legend("topleft", c("My_mes", "Mo_mes", "My_mic", "Mo_mic"), lty=1, col=1:4, bty="n")
# #matplot(times, Rt0,  type="l", lty=1, col=1:2, xlab="Time", ylab="Respiration", bty="n")

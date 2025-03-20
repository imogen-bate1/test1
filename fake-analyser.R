densities <- c(0,100,500,1000,3000)

load(paste("steady-CNE",densities[1],"pm2.r",sep=""))
OC <- SV[,"OC.f"] + SV[,"OC.s"] 
FeOOH.tot <- SV[,"FeOOH"] + SV[,"FeOOH.s"] 
SV <- cbind(SV,OC,FeOOH.tot)
CNE.SV <- array(dim=c(dim(SV),length(densities)),dimnames=list(dimnames(SV)[[1]],dimnames(SV)[[2]],densities))

load(paste("steady-PNE",densities[1],"pm2.r",sep=""))
OC <- SV[,"OC.f"] + SV[,"OC.s"] 
OC.b <- SV[,"OC.f.b"] + SV[,"OC.s.b"] 
FeOOH.tot <- SV[,"FeOOH"] + SV[,"FeOOH.s"] 
FeOOH.tot.b <- SV[,"FeOOH.b"] + SV[,"FeOOH.s.b"] 
SV <- cbind(SV,OC,OC.b,FeOOH.tot,FeOOH.tot.b)
PNE.SV <- array(dim=c(dim(SV),length(densities)),dimnames=list(dimnames(SV)[[1]],dimnames(SV)[[2]],densities))


for (i in 1:length(densities)){
  load(paste("steady-PNE",densities[i],"pm2.r",sep=""))
  OC <- SV[,"OC.f"] + SV[,"OC.s"] 
  OC.b <- SV[,"OC.f.b"] + SV[,"OC.s.b"] 
  FeOOH.tot <- SV[,"FeOOH"] + SV[,"FeOOH.s"] 
  FeOOH.tot.b <- SV[,"FeOOH.b"] + SV[,"FeOOH.s.b"] 
  SV <- cbind(SV,OC,OC.b,FeOOH.tot,FeOOH.tot.b)
  PNE.SV[,,i] <- SV
  load(paste("steady-CNE",densities[i],"pm2.r",sep=""))
  OC <- SV[,"OC.f"] + SV[,"OC.s"] 
  FeOOH.tot <- SV[,"FeOOH"] + SV[,"FeOOH.s"] 
  SV <- cbind(SV,OC,FeOOH.tot)
  CNE.SV[,,i] <- SV
}

spec.list <- c("OC", "O2", "FeOOH", "SO4")
spec.list.b <- paste(spec.list,".b",sep="")
spec.names <-c("OC^'tot'~mM",
               "O[2]~mM",
               "FeOOH^italic(a)~mM",
               "SO[4]^'2-'~mM"
               )

PNE.SV[,spec.list.b,"0"] <- NA

load("PNE 1-0-pm2.Rdata")
grid <- simulation.list[[1]]$PL$grid$x.mid

width <- length(spec.list)
height <- length(densities)

cex <- 0.6
x11(width=1.2*width+0.5, height=1.2*height)
par(mfcol=c(height,width), mar=c(1,1,.2,.2), oma=c(3,8,0,0),cex=cex)

for(i in spec.list) {
  j <- paste(i,".b",sep="")
  
  xmax <- max(CNE.SV[,i,],PNE.SV[,i,],PNE.SV[,j,],PNE.SV[,i,]*0.94+PNE.SV[,j,]*0.06,na.rm =T)
  xmin <- min(CNE.SV[,i,],PNE.SV[,i,],PNE.SV[,j,],na.rm =T)
  if(xmin < 0.5 * xmax) xmin <- 0  
  
  for(k in 1:length(densities)){
    
    load(paste("PNE 1",densities[k],"pm2.Rdata",sep="-"))
    burrow.coverage <- simulation.list[[1]]$PL$burrow.coverage
    grid <- simulation.list[[1]]$PL$grid$x.mid
    
    plot(CNE.SV[,i,k],grid,col="black",xlim=c(xmin,xmax),ylim=c(30,0),type="l",
         xlab="",ylab="",xaxt="n",yaxt="n", lwd=2)
    
    if(i == spec.list[2] && k == 1) {
      legend("bottomright",legend=c("CNE","PNE-bulk","PNE-burrow","PNE-total"),
             lty=c(1,2,3,1),col=c("black","gray","red","blue"),bty="n",
             lwd=c(2,2,2,2))
    }
    
    if(i == spec.list[1]){
      axis(2)
      mtext("depth (cm)",2,3,cex=cex)
      lab.text <- paste("theta==",densities[k],"~m^'-2'",sep="")
      mtext(parse(text=lab.text),2,4.5,cex=cex,las=1)
    } else {
      axis(2,labels=FALSE)
    }
    
    if(k == length(densities)){
      axis(1, las=2)
      mtext(parse(text=spec.names[match(i,spec.list)]),1,3,cex=cex)
    } else {
      axis(1,labels=FALSE)
    }
    

    lines(PNE.SV[,i,k]*(1-burrow.coverage)+PNE.SV[,j,k]*burrow.coverage,grid,col="blue", lty=1, lwd=2) 
    lines(PNE.SV[,i,k],grid,col="gray",lty=2, lwd=2) 
    lines(PNE.SV[,j,k],grid,col="red", lty=3, lwd=2)  
  }

}











spec.list <- c("Fe", "HS", "FeS", "FeS2")
spec.list.b <- paste(spec.list,".b",sep="")
spec.names <-c("Fe^'2+'~mM",
               "H[2]*S~mM",
               "FeS~mM",
               "FeS[2]~mM"
)

PNE.SV[,spec.list.b,"0"] <- NA

load("PNE 1-0-pm2.Rdata")
grid <- simulation.list[[1]]$PL$grid$x.mid

width <- length(spec.list)
height <- length(densities)

cex <- 0.6
x11(width=1.2*width+0.5, height=1.2*height)
par(mfcol=c(height,width), mar=c(1,1,.2,.2), oma=c(3,8,0,0),cex=cex)

for(i in spec.list) {
  j <- paste(i,".b",sep="")
  
  xmax <- max(CNE.SV[,i,],PNE.SV[,i,],PNE.SV[,j,],PNE.SV[,i,]*0.94+PNE.SV[,j,]*0.06,na.rm =T)
  xmin <- min(CNE.SV[,i,],PNE.SV[,i,],PNE.SV[,j,],na.rm =T)
  if(xmin < 0.5 * xmax) xmin <- 0  
  
  for(k in 1:length(densities)){
    
    load(paste("PNE 1",densities[k],"pm2.Rdata",sep="-"))
    burrow.coverage <- simulation.list[[1]]$PL$burrow.coverage
    grid <- simulation.list[[1]]$PL$grid$x.mid
    
    plot(CNE.SV[,i,k],grid,col="black",xlim=c(xmin,xmax),ylim=c(30,0),type="l",
         xlab="",ylab="",xaxt="n",yaxt="n", lwd=2)
    
    if(i == spec.list[1] && k == 1) {
      legend("bottomright",legend=c("CNE","PNE-bulk","PNE-burrow","PNE-total"),
             lty=c(1,2,3,1),col=c("black","gray","red","blue"),bty="n",
             lwd=c(2,2,2,2))
    }
    
    if(i == spec.list[1]){
      axis(2)
      mtext("depth (cm)",2,3,cex=cex)
      lab.text <- paste("theta==",densities[k],"~m^'-2'",sep="")
      mtext(parse(text=lab.text),2,4.5,cex=cex,las=1)
    } else {
      axis(2,labels=FALSE)
    }
    
    if(k == length(densities)){
      axis(1, las=2)
      mtext(parse(text=spec.names[match(i,spec.list)]),1,3,cex=cex)
    } else {
      axis(1,labels=FALSE)
    }
    
    
    lines(PNE.SV[,i,k]*(1-burrow.coverage)+PNE.SV[,j,k]*burrow.coverage,grid,col="blue", lty=1, lwd=2) 
    lines(PNE.SV[,i,k],grid,col="gray",lty=2, lwd=2) 
    lines(PNE.SV[,j,k],grid,col="red", lty=3, lwd=2)  
  }
  
}
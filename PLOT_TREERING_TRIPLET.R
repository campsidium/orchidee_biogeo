# SCRIPT. Read tree-ring width and isotope data to build the tree-ring triplet surface response.
#
# Dependencies: openxlsx and mgcv packages
# 
#
# @ Dr. Jonathan Barichivich, LSCE, Paris.
# @ oreobulus@gmail.com
# @ 19 June 2021

# --- Set parameters
saveto='./'              # path to save figure
savefigs = 0             # 1 = save figure, 0 = don't save figure, just plot
data_period=c(1953,2000) # period to use for the tree-ring triplet


# --- Load observations simulations
require(openxlsx)
infile='./DATA_bg-2020-446.xlsx'
OBS=read.xlsx(infile,sheet=2, startRow = 2)
OBS_header=read.xlsx(infile,sheet=2, startRow = 1, colNames=F, rows = c(1,2) )
ORCHIDEE=read.xlsx(infile,sheet=3, startRow = 2)
ORCHIDEE_header=read.xlsx(infile,sheet=3, startRow = 1, colNames=F, rows = c(1,2) )


# --- Build observation triplet
idx=which(OBS[,1]>=data_period[1] & OBS[,1]<=data_period[2])
yi=as.numeric(OBS[idx,6])
xi=as.numeric(OBS[idx,12])
zi=as.numeric(OBS[idx,18])
year=as.numeric(OBS[idx,1])
Xobs = as.data.frame(cbind(yy,xi,yi,zi)) 
colnames(Xobs)=c('year','d18O','d13C','TRW')


# --- Build simulated triplet
idx=which(ORCHIDEE[,1]>=data_period[1] & ORCHIDEE[,1]<=data_period[2])
yi=as.numeric(ORCHIDEE[idx,6])
xi=as.numeric(ORCHIDEE[idx,12])
zi=as.numeric(ORCHIDEE[idx,18])
year=as.numeric(ORCHIDEE[idx,1])
Xsim = as.data.frame(cbind(yy,xi,yi,zi)) 
colnames(Xsim)=c('year','d18O','d13C','TRW')


# --- Reset mean of simulations to compare variability
mu0=mean(Xobs$TRW,na.rm=T)
mu1=mean(Xsim$TRW,na.rm=T)
Xsim$TRW = (Xsim$TRW-mu1) + mu0
#
mu0=mean(Xobs$d13C,na.rm=T)
mu1=mean(Xsim$d13C,na.rm=T)
Xsim$d13C = (Xsim$d13C-mu1) + mu0
#
mu0=mean(Xobs$d18O,na.rm=T)
mu1=mean(Xsim$d18O,na.rm=T)
Xsim$d18O = (Xsim$d18O-mu1) + mu0


########################
# FIGURE TREE-RING TRIPLET
########################
library(mgcv)
fout=paste(saveto, "fig01_FONTAINEBLEAU_triplet_surface_response.png", sep="")
if (savefigs==1) png(fout, width=10, height=10, units="in", res=360) # save as png
#if (savefigs==1) cairo_pdf(filename=fout, height = 10, width = 10) # save as eps to edit vectors elsewhere
#if (savefigs==1) setEPS(); postscript(fout,width=10, height=10)
par(fig=c(0.0,0.50,0.,1),cex.main=0.9, mar=c(4,4,1,1), mgp=c(2,0.5,0), tcl=0.2, font=2)
X=Xsim
model=gam(data=X, TRW~t2(d18O, d13C,k=3), method = 'REML',family = scat())
#gam.chack(model)
x <-range(X$d18O)
x <- seq(x[1], x[2], length.out=20)    
y <- range(X$d13C)
y <- seq(y[1], y[2], length.out=20)
z <- outer(x,y, 
           function(d18O,d13C)
             predict(model, data.frame(d18O,d13C)))
p <- persp(x,y,z, theta=40, phi=30,    # phi is angle for vertical
           col="white",expand = 1.,shade = 0.1,
           xlab="Oxygen", ylab="Carbon", zlab="TRW",ticktype="detailed",zlim=c(0.,1.5),xlim=c(29,33),ylim=c(15,19), main='ORCHIDEE')
obs <- trans3d(X$d18O, X$d13C,X$TRW,p)
pred <- trans3d(X$d18O, X$d13C,fitted(model),p)
points(obs, col="gray10",pch=16)
segments(obs$x, obs$y, pred$x, pred$y,col='gray10',lty=1)
#
# overplot 1976 extreme 
idx=which(X$year==1976)
obs <- trans3d(X$d18O[idx], X$d13C[idx],X$TRW[idx],p)
points(obs, col="red",pch=16)

par(new=TRUE)
par(fig=c(0.50,1,0.,1),cex.main=0.9, mar=c(4,4,1,1), mgp=c(2,0.5,0), tcl=0.2, font=2)
X=Xobs
model=gam(data=X, TRW~t2(d18O, d13C,k=3), method = 'REML',family = scat())
#gam.chack(model)
x <-range(X$d18O)
x <- seq(x[1], x[2], length.out=20)    
y <- range(X$d13C)
y <- seq(y[1], y[2], length.out=20)
z <- outer(x,y, 
           function(d18O,d13C)
             predict(model, data.frame(d18O,d13C)))
p <- persp(x,y,z, theta=40, phi=30,    # phi is angle for vertical
           col="white",expand = 1.,shade = 0.1,
           xlab="Oxygen", ylab="Carbon", zlab="TRW",ticktype="detailed",zlim=c(0.,1.5),xlim=c(29,33),ylim=c(15,19), main='OBSERVATIONS')
obs <- trans3d(X$d18O, X$d13C,X$TRW,p)
pred <- trans3d(X$d18O, X$d13C,fitted(model),p)
points(obs, col="gray10",pch=16)
segments(obs$x, obs$y, pred$x, pred$y,col='gray10',lty=1)
#
# overplot 1976 extreme 
idx=which(X$year==1976)
obs <- trans3d(X$d18O[idx], X$d13C[idx],X$TRW[idx],p)
points(obs, col="red",pch=16)
if (savefigs==1) dev.off() 

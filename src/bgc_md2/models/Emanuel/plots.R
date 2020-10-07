library(RColorBrewer)
library(lattice)

pal=brewer.pal(5,"Dark2")

stocks=read.csv("stocks.csv")
meanAges=read.csv("mean_ages.csv")
medianAges=read.csv("median_ages.csv")
ftt=read.csv("forward_transit_time.csv")
agedensEq=read.csv("start_age_dens.csv")
agedens=read.csv("Emmanuel_age_dens.csv")
bttdens=read.csv('backward_transit_time_densities.csv')
bttmeans=read.csv('backward_transit_time_mean_and_median.csv')

pool_names=c("Non-woody tree parts", "Woody tree parts", "Ground vegetation", "Detritus/Decomposers", "Soil" )

calyr=stocks[,6]

pdf("Figures/poolStocks.pdf")
matplot(calyr,stocks[,-6], type="l", lty=1, col=pal, ylim=c(0,1500),xlab="Calendar year", ylab="Carbon stocks (PgC)",bty="n")
legend("topright", pool_names,col=pal,lty=1,bty="n")
dev.off()

TCS=rowSums(stocks[,-6])
difTCS=max(TCS)-min(TCS)

pdf("Figures/totalStocks.pdf")
par(mar=c(4,4,1,0))
plot(calyr,TCS, type="l", lwd=2, ylim=c(1500,2000), xlim=c(1850, 2020),xlab="Calendar year", ylab="Carbon stock (PgC)", bty="n")
arrows(2016, head(TCS, 1), 2016, tail(TCS,1), angle=90, code=3, length=0.05)
text(2020,min(TCS)+(difTCS/2), labels=round(difTCS, 1))
dev.off()

pdf("Figures/meanPoolAges.pdf")
matplot(calyr, meanAges[,-6], type="l", col=pal, lty=1, ylim=c(0,120),bty="n", xlab="Calendar year", ylab="Mean ages (yr)")
legend(1950,80, pool_names,col=pal,lty=1,bty="n")
dev.off()

pdf("Figures/medianPoolAges.pdf")
matplot(calyr, medianAges[,-6], type="l", col=pal, lty=1, bty="n", xlab="Calendar year", ylab="Median ages (yr)")
legend(1850,60, pool_names,col=pal,lty=1,bty="n")
dev.off()

pdf("Figures/systemAges.pdf")
par(mar=c(4,4,1,0))
plot(calyr, meanAges$System.Age,type="l", col=2,bty="n", ylim=c(0,80),xlab="Calendar year", ylab="Carbon age (yr)")
lines(calyr, medianAges$SystemA.Age, col=4)
legend("bottomleft", c("Mean system age", "Median system age"), lty=1, col=c(2,4), bty="n")
dev.off()

ages=seq(0,nrow(ftt)-1, by=1)
years=seq(1851,2000)

pdf("Figures/forwardTransitTime3d.pdf")
wireframe(t(as.matrix(ftt)), shade=TRUE, ylim=c(0,50), xlab="Calendar year", ylab="Age (yr)", zlab="C stock \n (PgC)", aspect = c(61/87, 0.4),
          main="Forward transit time")
dev.off()

pdf("Figures/backwardTransitTime3d.pdf")
wireframe(t(as.matrix(bttdens)), shade=TRUE, ylim=c(0,50), xlab="Calendar year", ylab="Age (yr)", zlab="C stock \n (PgC)", aspect = c(61/87, 0.4),
          main="Backward transit time")
dev.off()

pdf('Figures/meanBackwardTransitTime.pdf')
plot(bttmeans$Year, bttmeans$meanBTT, type="l", xlab="Calenday year", ylab="Mean backward transit time (yr)", bty="n")
dev.off()

pdf('Figures/medianBackwardTransitTime.pdf')
plot(bttmeans$Year, bttmeans$medianBTT, type="l", xlab="Calenday year", ylab="Median backward transit time (yr)", bty="n")
dev.off()

pdf("Figures/poolAgesEq.pdf")
par(mar=c(4,4,1,0))
matplot(ages,agedensEq, type="l", lty=1, col=pal, bty="n", xlab="Ages (yr)", ylab="Carbon stock (PgC)")
legend("topright", pool_names, lty=1, col=pal, bty="n")
dev.off()

x1=subset(agedens,pool==0, select=-pool)
x2=subset(agedens,pool==1, select=-pool)
x3=subset(agedens,pool==2, select=-pool)
x4=subset(agedens,pool==3, select=-pool)
x5=subset(agedens,pool==4, select=-pool)
sysAge=subset(agedens,pool==-1, select=-pool)

pdf("Figures/nonwoodyTreeParts3d.pdf")
wireframe(value~time*age, data=x1, shade=TRUE, xlab="Calendar year", ylab="Age (yr)", zlab="C stock \n (PgC)", ylim=c(0,10), 
          scales = list(arrows = FALSE), main=pool_names[1])
dev.off()

pdf("Figures/woodyTreeParts3d.pdf")
wireframe(value~time*age, data=x2, shade=TRUE, xlab="Calendar year", ylab="Age (yr)", zlab="C stock \n (PgC)", ylim=c(0,50), 
          scales = list(arrows = FALSE), main=pool_names[2])
dev.off()

pdf("Figures/groundVeg3d.pdf")
wireframe(value~time*age, data=x3, shade=TRUE, xlab="Calendar year", ylab="Age (yr)", zlab="C stock \n (PgC)", ylim=c(0,30), 
          scales = list(arrows = FALSE), main=pool_names[3])
dev.off()

pdf("Figures/detritus3d.pdf")
wireframe(value~time*age, data=x4, shade=TRUE, xlab="Calendar year", ylab="Age (yr)", zlab="C stock \n (PgC)", ylim=c(0,30), 
          scales = list(arrows = FALSE), main=pool_names[4])
dev.off()

pdf("Figures/soil3d.pdf")
wireframe(value~time*age, data=x5, shade=TRUE, xlab="Calendar year", ylab="Age (yr)", zlab="C stock \n (PgC)", ylim=c(0,150), 
          scales = list(arrows = FALSE), main=pool_names[5])
dev.off()

pdf("Figures/systemAge3d.pdf")
wireframe(value~time*age, data=sysAge, shade=TRUE, xlab="Calendar year", ylab="Age (yr)", zlab="C stock \n (PgC)", ylim=c(0,100), 
          scales = list(arrows = FALSE), main="System Age")
dev.off()

longftt=as.vector(ftt)          

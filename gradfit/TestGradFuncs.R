<<<<<<< HEAD
####  	Example of gcv.rss and predict.rss 
# simulated data (x,y) and true mean response  

x<-runif(100)*3-0.1; x<-x[order(x)];
ybar<-5*(exp(-x)-4*exp(-2*x)+3*exp(-3*x)); 
y<-ybar+0.3*rnorm(100);

# fit & plot with penalty=2
fit2<-gcv.rss(x,y,hmin=-5,hmax=5,intercept=T,penalty=2,draw=T);
points(x,ybar,type="l",lty=2,lwd=2);

# fit & don't plot with penalty=1
fit1<-gcv.rss(x,y,hmin=-6,hmax=1,intercept=T,penalty=1,draw=F); 

# predict & plot new values from the model with penalty=1
newx<-range(x)[1]+0.01*c(0:100)*(range(x)[2]-range(x)[1]); 
yhat1<-predict.rss(fit1$model,newx);
points(yhat1$x,yhat1$y,type="l",lty=4,lwd=1); 
legend(0.25,1,c("true","penalty=2","penalty=1"),lty=c(2,1,4),lwd=c(2,1,1),cex=0.75); 

##############

############## Example of gcv.rssM 
x<-2*runif(250); x<-x[order(x)]; ybar<-(x^20)/(1+x^20); 
y<-ybar+0.5*rnorm(250);
fitM<-gcv.rssM(x,y,hmin=-12,hmax=0,intercept=F,penalty=1,draw=T); 
points(x,ybar,type="l",lty=2,lwd=2);

# fit & plot the same data without monotonicity constraint
fit<-gcv.rss(x,y,hmin=-12,hmax=0,intercept=F,penalty=1,draw=F); 
points(x,fit$fitted.values,type="l",lty=4,lwd=1); 
legend(0,2.2,c("true","monotonic","unconstrained"),lty=c(2,1,4),lwd=c(2,1,1),cex=0.75); 

##############


################# Example of bcgrad and gcv.rssadd2 on Nicholson data
library(mgcv); par(mfcol=c(2,2)); 

#load data and estimate the gradient 
x<-as.matrix(read.table("c:\\GRADFuncs\\NichAdults.txt")); 
ad<-x/1000; tvals<-2*(c(1:length(ad))-1); 

dxhat<-bcgrad(tvals,ad,hder=2); 
plot(tvals,1000*ad,type="p",xlab="Time (d)",ylab="Adults"); 
points(tvals,1000*dxhat$smooth,type="l",lty=1); title("(a) Time series & smooth"); 
zhat<-dxhat$unbiased; zraw<-dxhat$raw; 

# plot gradient bias correction curve 
plot(1000*zraw,1000*zhat,type="p",xlab="Raw gradient estimate",ylab="Bias corrected estimate");
ez<-order(zraw); 
points(1000*zraw[ez],1000*zraw[ez],type="l",lty=1);
title(" (b) Gradient bias correction"); 

## Create data vectors: xt1=current density, xt0=lagged density, dadt=current gradient
xt1<-trimr(ad,9,2); xt0<-trimr(ad,2,9); dadt<-trimr(zhat,9,2); 

## Unconstrained regression spline fit using gam() from mgcv library
gamfit<-gam(dadt~s(xt1,20)+s(xt0,20));
newd1<-data.frame(xt1=(0:120)/20,xt0=0*(0:120)/20);
f1<- -predict.gam(gamfit,newd1); f1<-f1-f1[1];
newd2<-data.frame(xt1=0*(0:120)/20,xt0=(0:120)/20);
f0<-predict.gam(gamfit,newd2); f0<-f0-f0[1];
matplot(1000*(0:120)/20, 1000*(f1%~%f0),type="l",lty=1,xlab="Number of adults",
ylab="Fitted B and D rates"); title("(c) Unconstrained GAM fit");


## Sign-constrained regression spline fit using gcv.rssadd2
## Note: gam() from mgcv is MUCH FASTER because we use a general-purpose
## optimizer (R's optim() with method="Nelder-Mead") for minimizing GCV(alpha). 
##  mgcv is built for speed, this set of functions was built speedily. 
rssfit<-gcv.rssadd2(x=xt1%~%xt0,y=dadt,Boundary.knots1=c(0,range(xt1)[2]),
Boundary.knots2=c(0,range(xt0)[2]), signs=c(-1,1),intercept1=F,
intercept2=F,penalty=2);

f1c<- -1000*rssfit$fitted.values[,1];
f0c<- 1000*rssfit$fitted.values[,2]; 
e0<-order(xt0); 
plot(1000*c(0,xt0[e0]),c(0,f0c[e0]),type="l",lty=1,xlab="Number of adults x 1000",
ylab="Fitted B and D rates",xlim=c(0,6000)); title("(d) Sign-constrained fit");
e1<-order(xt1); 
points(1000*c(0,xt1[e1]),c(0,f1c[e1]),type="l",lty=1)


## SIMEX example ########################################################################
## This uses some of the 'utility' routines that are not described in the documentation,
## but the comments in the source code should be sufficient. 
##	rssfit() is similar to gcv.rss except that the smoothing parameter is user-supplied
##	rssbasis() generates the spline basis functions at a set of points, given the knots
##	  and the left endpoint of the interval ("origin" in the call). 

# the "data"; x and ybar are the unobserved truth, xobs and y are the error-corrupted observations.
x<-5*runif(200)-2.5; x<-x[order(x)]; 
ybar<-exp(-(x^2)); y<-ybar+0.01*rnorm(200); 
sigme<-0.35;
xobs<-x+sigme*rnorm(200);    
plot(xobs,y,xlim=c(-3,3), ylim=c(0,1.2),xlab="Observed x values",ylab="Observed y values"); 
points(x,ybar,type="l",lty=1,lwd=2); 

# fit using the observed data
rawfit<-gcv.rss(xobs,y,hmin=-6,hmax=6,intercept=T,penalty=1,draw=F);
exobs<-order(xobs); 
points(xobs[exobs],rawfit$fitted.values[exobs],type="l",lty=2,lwd=2); 

# "freeze" the smoothing parameter at the optimum for the observed data
hopt<-rawfit$hopt;

# extract other parameters of the fit to the observed data
nknots<-length(rawfit$model$knots);
Boundary.knots<-rawfit$model$Boundary.knots;

# refit 1000 times with one added dose of measurement error
fitted2<-0*rawfit$model$coef; fitted3<-fitted2;
for (irep in 1:1000) {
	ei<-sigme*rnorm(200);
	x2<-xobs+ei;
	fit2<-rssfit(x2,y,nknots=nknots,Boundary.knots=Boundary.knots,alpha=hopt,intercept=T); 
	fitted2<-fitted2+fit2$coef;
	x3<-xobs+sqrt(2)*ei; 
	fit3<-rssfit(x3,y,nknots=nknots,Boundary.knots=Boundary.knots,alpha=hopt,intercept=T); 
	fitted3<-fitted3+fit3$coef;
	cat(c(irep, "\n")); 
}
fitted2<-fitted2/1000; fitted3<-fitted3/1000; 

SimexCoef<-6*rawfit$model$coef + 4*fitted3 - 9*fitted2;

xplot<- xobs[order(xobs)]; knots<- rawfit$model$knots; origin<-rawfit$model$origin; 
X<-rssbasis(xplot,knots=knots,origin=origin,intercept=T); 
points(xplot,X%*%SimexCoef,type="l",lty=4,lwd=3); 
legend(-3,1,c("True","Raw fit","SIMEX"),lty=c(1,2,4),lwd=c(2,2,3)); 
#####################################################################

=======
####  	Example of gcv.rss and predict.rss 
# simulated data (x,y) and true mean response  

x<-runif(100)*3-0.1; x<-x[order(x)];
ybar<-5*(exp(-x)-4*exp(-2*x)+3*exp(-3*x)); 
y<-ybar+0.3*rnorm(100);

# fit & plot with penalty=2
fit2<-gcv.rss(x,y,hmin=-5,hmax=5,intercept=T,penalty=2,draw=T);
points(x,ybar,type="l",lty=2,lwd=2);

# fit & don't plot with penalty=1
fit1<-gcv.rss(x,y,hmin=-6,hmax=1,intercept=T,penalty=1,draw=F); 

# predict & plot new values from the model with penalty=1
newx<-range(x)[1]+0.01*c(0:100)*(range(x)[2]-range(x)[1]); 
yhat1<-predict.rss(fit1$model,newx);
points(yhat1$x,yhat1$y,type="l",lty=4,lwd=1); 
legend(0.25,1,c("true","penalty=2","penalty=1"),lty=c(2,1,4),lwd=c(2,1,1),cex=0.75); 

##############

############## Example of gcv.rssM 
x<-2*runif(250); x<-x[order(x)]; ybar<-(x^20)/(1+x^20); 
y<-ybar+0.5*rnorm(250);
fitM<-gcv.rssM(x,y,hmin=-12,hmax=0,intercept=F,penalty=1,draw=T); 
points(x,ybar,type="l",lty=2,lwd=2);

# fit & plot the same data without monotonicity constraint
fit<-gcv.rss(x,y,hmin=-12,hmax=0,intercept=F,penalty=1,draw=F); 
points(x,fit$fitted.values,type="l",lty=4,lwd=1); 
legend(0,2.2,c("true","monotonic","unconstrained"),lty=c(2,1,4),lwd=c(2,1,1),cex=0.75); 

##############


################# Example of bcgrad and gcv.rssadd2 on Nicholson data
library(mgcv); par(mfcol=c(2,2)); 

#load data and estimate the gradient 
x<-as.matrix(read.table("c:\\GRADFuncs\\NichAdults.txt")); 
ad<-x/1000; tvals<-2*(c(1:length(ad))-1); 

dxhat<-bcgrad(tvals,ad,hder=2); 
plot(tvals,1000*ad,type="p",xlab="Time (d)",ylab="Adults"); 
points(tvals,1000*dxhat$smooth,type="l",lty=1); title("(a) Time series & smooth"); 
zhat<-dxhat$unbiased; zraw<-dxhat$raw; 

# plot gradient bias correction curve 
plot(1000*zraw,1000*zhat,type="p",xlab="Raw gradient estimate",ylab="Bias corrected estimate");
ez<-order(zraw); 
points(1000*zraw[ez],1000*zraw[ez],type="l",lty=1);
title(" (b) Gradient bias correction"); 

## Create data vectors: xt1=current density, xt0=lagged density, dadt=current gradient
xt1<-trimr(ad,9,2); xt0<-trimr(ad,2,9); dadt<-trimr(zhat,9,2); 

## Unconstrained regression spline fit using gam() from mgcv library
gamfit<-gam(dadt~s(xt1,20)+s(xt0,20));
newd1<-data.frame(xt1=(0:120)/20,xt0=0*(0:120)/20);
f1<- -predict.gam(gamfit,newd1); f1<-f1-f1[1];
newd2<-data.frame(xt1=0*(0:120)/20,xt0=(0:120)/20);
f0<-predict.gam(gamfit,newd2); f0<-f0-f0[1];
matplot(1000*(0:120)/20, 1000*(f1%~%f0),type="l",lty=1,xlab="Number of adults",
ylab="Fitted B and D rates"); title("(c) Unconstrained GAM fit");


## Sign-constrained regression spline fit using gcv.rssadd2
## Note: gam() from mgcv is MUCH FASTER because we use a general-purpose
## optimizer (R's optim() with method="Nelder-Mead") for minimizing GCV(alpha). 
##  mgcv is built for speed, this set of functions was built speedily. 
rssfit<-gcv.rssadd2(x=xt1%~%xt0,y=dadt,Boundary.knots1=c(0,range(xt1)[2]),
Boundary.knots2=c(0,range(xt0)[2]), signs=c(-1,1),intercept1=F,
intercept2=F,penalty=2);

f1c<- -1000*rssfit$fitted.values[,1];
f0c<- 1000*rssfit$fitted.values[,2]; 
e0<-order(xt0); 
plot(1000*c(0,xt0[e0]),c(0,f0c[e0]),type="l",lty=1,xlab="Number of adults x 1000",
ylab="Fitted B and D rates",xlim=c(0,6000)); title("(d) Sign-constrained fit");
e1<-order(xt1); 
points(1000*c(0,xt1[e1]),c(0,f1c[e1]),type="l",lty=1)


## SIMEX example ########################################################################
## This uses some of the 'utility' routines that are not described in the documentation,
## but the comments in the source code should be sufficient. 
##	rssfit() is similar to gcv.rss except that the smoothing parameter is user-supplied
##	rssbasis() generates the spline basis functions at a set of points, given the knots
##	  and the left endpoint of the interval ("origin" in the call). 

# the "data"; x and ybar are the unobserved truth, xobs and y are the error-corrupted observations.
x<-5*runif(200)-2.5; x<-x[order(x)]; 
ybar<-exp(-(x^2)); y<-ybar+0.01*rnorm(200); 
sigme<-0.35;
xobs<-x+sigme*rnorm(200);    
plot(xobs,y,xlim=c(-3,3), ylim=c(0,1.2),xlab="Observed x values",ylab="Observed y values"); 
points(x,ybar,type="l",lty=1,lwd=2); 

# fit using the observed data
rawfit<-gcv.rss(xobs,y,hmin=-6,hmax=6,intercept=T,penalty=1,draw=F);
exobs<-order(xobs); 
points(xobs[exobs],rawfit$fitted.values[exobs],type="l",lty=2,lwd=2); 

# "freeze" the smoothing parameter at the optimum for the observed data
hopt<-rawfit$hopt;

# extract other parameters of the fit to the observed data
nknots<-length(rawfit$model$knots);
Boundary.knots<-rawfit$model$Boundary.knots;

# refit 1000 times with one added dose of measurement error
fitted2<-0*rawfit$model$coef; fitted3<-fitted2;
for (irep in 1:1000) {
	ei<-sigme*rnorm(200);
	x2<-xobs+ei;
	fit2<-rssfit(x2,y,nknots=nknots,Boundary.knots=Boundary.knots,alpha=hopt,intercept=T); 
	fitted2<-fitted2+fit2$coef;
	x3<-xobs+sqrt(2)*ei; 
	fit3<-rssfit(x3,y,nknots=nknots,Boundary.knots=Boundary.knots,alpha=hopt,intercept=T); 
	fitted3<-fitted3+fit3$coef;
	cat(c(irep, "\n")); 
}
fitted2<-fitted2/1000; fitted3<-fitted3/1000; 

SimexCoef<-6*rawfit$model$coef + 4*fitted3 - 9*fitted2;

xplot<- xobs[order(xobs)]; knots<- rawfit$model$knots; origin<-rawfit$model$origin; 
X<-rssbasis(xplot,knots=knots,origin=origin,intercept=T); 
points(xplot,X%*%SimexCoef,type="l",lty=4,lwd=3); 
legend(-3,1,c("True","Raw fit","SIMEX"),lty=c(1,2,4),lwd=c(2,2,3)); 
#####################################################################

>>>>>>> 7aa694dd646fa18774378eaa5eea76f04c87ed74

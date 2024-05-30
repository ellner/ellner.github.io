#### ----------------------     Introduction 
## R functions for 
## (0) making R understand a bit of GAUSS. Some of this code was ported from GAUSS.
## (1) Estimating the gradient dx/dt for a time series x(t) with observations
##     at times t0<t1<t2<... 
## (2) Fitting penalized regression splines of several types:
##     univariate, univariate monotonic increasing, and additive with sign
##     constraints. Fitting criterion is GCV with user-specified penalty parameter 
##     as in smooth.spline (the default, penalty=1, gives the usual GCV criterion)

## To run under Splus, uncomment the following line: 
## 	library(Matrix); solve<-solve.Matrix;

## Not everything runs under Splus, though adapting them to do so should 
## not be too difficult. The main R-only feature is the use of 
## optim() with method = "Nelder-Mead" to fit the two smoothing
## parameters for the ridge functions in additive models; a call to
## to nlmin() could be substituted, but only if you're willing to trust 
## results from nlmin. 


#### -------------------------   CODE STARTS HERE

library(quadprog); 
library(MASS); 

## Utility routines betraying that this code was ported from GAUSS.

"%~%"<-function(a,b) {cbind(a,b)}
"%|%"<-function(a,b) {rbind(a,b)}
minindc<-function(x) {match(range(x)[1],x) };

stack<-function(xx, yy){
	nx <- length(xx)
	ny <- length(yy)
	z <- matrix(0, nx + ny, 1)
	z[1:nx] <- xx[1:nx]
	z[(nx + 1):(nx + ny)] <- yy[1:ny]
	return(z)
}

trimr<-function (a,n1,n2) {
	da<-dim(a); 
	if(is.null(da)) {a[(n1+1):(length(a)-n2)]}
	else {a[(n1+1):(da[1]-n2),]};
};

zeros<-function(nrow,ncol) {matrix(0,nrow,ncol)};
ones<-function(nrow,ncol) {matrix(1,nrow,ncol)};
rndn<-function(nrow,ncol) {matrix(rnorm(nrow*ncol),nrow,ncol)};
rndu<-function(nrow,ncol) {matrix(runif(nrow*ncol),nrow,ncol)};
rows<-function(a) {
	da<-dim(a)
	if(is.null(da)) {length(a)} else {da[1]}
}
cols<-function(a) {
	da<-dim(a)
	if(is.null(da)) {1} else {da[2]}
}

xy<-function(x,y,...) {
	e<-order(x);
	plot(x[e],y[e],...);
}

xypoints<-function(x,y,...) {
	e<-order(x);
	points(x[e],y[e],...);
}

## Routines for gradient estimation via local polynomial regression, 
## Gaussian kernel truncated at (-10*h,10*h) 

# Utility routine for locregG. Does the fit at a single point.		
locreg1G<-function(xvals,yvals,x,h,degree=2) {
	x0<-(xvals-x);
	e<-(abs(x0)<=10*h);
	yvals<-yvals[e]; x0<-x0[e];
	X<-matrix(0,length(x0),degree); 
	for (i in 1:degree) {
		X[,i]<-x0^i;
	}
	w<-exp(-((x0/h)^2));
	betahat<-lsfit(X,yvals,wt=w,intercept=T)$coef;
	return(as.matrix(betahat));
}

# locregG: given data (x,y), produce a fitted value at a vector
# of locations xnew by local regression smoothing with bandwidth h 
# and selected degree (default 2).
# RETURNS a vector ysmu whose ith column is the ith-degree coefficient
# in the local polynomial fit centered at the data points. In particular
# ysmu[,1] is the fitted value, ysmu[,2] is the estimated derivative.
locregG<-function(x,y,xnew,h,degree=2) {
	ysmu<-matrix(0,length(xnew),degree+1); 
	for (i in 1:length(xnew)) {
		xi<-xnew[i]
		bi<-locreg1G(xvals=x,yvals=y,x=xi,h=h,degree=degree)
		ysmu[i,]<-matrix(bi,1,degree+1)
	}
	return(ysmu)
}


## Routines for fitting penalized cubic regression splines by GCV
## "Penalty" is the cost per effective model df, as in smooth.spline(); 

# The "plus" regression spline basis
# x are the valuation points, knots are the knots, origin is the left
# endpoint of the interval corresponding to Boundary.knots[1] in the
# functions that use the basis, and intercept is logical for inclusion
# of a constant term in the basis
rssbasis<-function(x, knots, origin, intercept) {
    	xs<-(x-origin);
	p<-cbind(xs,xs^2,xs^3); 
	if(intercept) p<-cbind( matrix(1,length(x),1),p); 
	nknots<-length(knots); 
	for (i in 1:nknots) {
	  si<-pmax(x-knots[i],0);
	  p<-cbind(p,si^3);
	}
	return(as.matrix(p))
}

# Fit Penalized Regression Spline with user-supplied smoothing parameter
#	NOTE alpha is the log10 of the smoothing parameter
#	Other parameters are the same as in gcv.rss 
rssfit<-function(x,y,nknots,Boundary.knots,intercept=T,alpha,penalty=1)	{
	origin<-Boundary.knots[1]; 
	fac<-(Boundary.knots[2]-Boundary.knots[1])/(nknots+1)
	knots<- origin+c(1:nknots)*fac; 
	X<-rssbasis(x,knots,origin,intercept);
	XTX<-t(X)%*%X;
	d<-rep(1,dim(XTX)[1]); 
	d[1:3]<-0; if(intercept) d[4]<-0; 
	D<-diag(d); 
	betahat<-solve(XTX+(10^alpha)*D,t(X)%*%y,tol=1.e-16);
	yhat<-X%*%betahat;
	
	msr<-mean((yhat-y)^2); 
	hat<-X%*%solve(XTX+(10^alpha)*D,tol=1.e-16)%*%t(X);
	d1<-sum(diag(hat))/length(y);
	gcv<-msr/((1-penalty*d1)^2);
	if(penalty*d1>1) gcv<-1.e12; 
	return(list(x=x,knots=knots,origin=origin,alpha=alpha,Boundary.knots=Boundary.knots,
		intercept=intercept, coef=betahat, fitted.values=yhat,gcv=gcv,enp=sum(diag(hat))))
}

predict.rss<-function(model,newx){
	X<-rssbasis(newx,model$knots,model$origin,model$intercept);
	yhat<-X%*%(model$coef);
	return(list(x=newx,y=yhat))
}


# Fit Penalized Regression Spline with smoothing parameter chosen by GCV(p)
# Hmin and Hmax are user-supplied brackets on log10(optimal smoothing parameter) 
# Uses golden section search so hmin and hmax MUST bracket a minimum or the
# routine will return garbage. Use draw=T (the default) to make sure you have a bracket. 
gcv.rss<-function(xvals,yvals,hmin,hmax,nknots=NA,Boundary.knots=NA,intercept=T,penalty=1,tol=0.01,draw=T) {
	if(is.na(Boundary.knots)) Boundary.knots<-range(xvals)
	if(is.na(nknots)) nknots<-20;
	C<-0.61803399; R<-1-C;
	dh<-hmax-hmin;
	hvals<-hmin+0.1*c(0:10)*dh;
	Gvals<-matrix(0,11,1);
	for (i in 1:11) {
		Gvals[i]<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,
			hvals[i],penalty=penalty)$gcv;
		print(c(i,hvals[i],Gvals[i]));
        }
	imin<-minindc(Gvals);
	if((imin<=1)||(imin>=11)) {hopt<-hvals[imin]} else
	{ 
	ax<-hvals[imin+1]; bx<-hvals[imin]; cx<-hvals[imin-1]; 
	x0<-ax; x3<-cx; 
	if(abs(cx-bx)>abs(bx-ax)) 
		{x1<-bx; x2<-bx+C*(cx-bx)} 
	else {x2<-bx; x1<-bx-C*(bx-ax)}
	f1<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x1)$gcv;	
	f2<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x2)$gcv;
	while(abs(x3-x0)>tol*(abs(x1)+abs(x2))) {
		if(f2<f1) {x0<-x1; x1<-x2; x2<-R*x1+C*x3; f1<-f2; 
			f2<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x2)$gcv;}
		else {x3<-x2; x2<-x1; x1<-R*x2+C*x0; f2<-f1;
			f1<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x1)$gcv;}
	}		
	if(f1<f2) {hopt<-x1; fmin<-f1} else {hopt<-x2; fmin<-f2}
	print(c(hopt,fmin));
	}

	bestfit<-rssfit(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=hopt);	

	if (draw==T) {
	win.graph(); par(mfrow=c(2,1));
	plot(hvals,log10(Gvals),type="l",xlab="log10(h)",
		ylab="GCV");
	title("Bandwidth selection curve");
	ox<-order(xvals);
	plot(xvals,yvals); points(xvals[ox],bestfit$fitted.values[ox],type="l");
	title("Data and P-spline regression curve");
	}

	return(list(hvals=hvals,Gvals=Gvals,hopt=hopt,fitted.values=bestfit$fitted.values,
		penalty=penalty,gcv=bestfit$gcv,model=bestfit))
       				
}

# Use Qprog to Fit Penalized Regression Spline with user-supplied 
# smoothing parameter, constrained to be monotonic at the knots
# NOTE x values are sorted into increasing order on return. If x is unsorted
# on input, plotting fitted values versus input x will give garbage. 
rssfitM<-function(x,y,nknots,Boundary.knots,intercept=T,penalty=1,alpha,ncongrid)	{
	ex<-order(x); x<-x[ex]; y<-y[ex]; 
	origin<-Boundary.knots[1]; 
	fac<-(Boundary.knots[2]-Boundary.knots[1])/(nknots+1)
	knots<- origin+c(1:nknots)*fac; 
	X<-rssbasis(x,knots,origin,intercept);
	XTX<-t(X)%*%X;

	d<-rep(1,dim(XTX)[1]); 
	d[1:3]<-0; if(intercept) d[4]<-0;
	D<-diag(d); 

# Quadratic programming objective function
	qmat<-XTX+(10^alpha)*D; rmat<-as.vector(t(X)%*%y);  

# Impose monotonicity at all constraining grid points
	cfac<-(Boundary.knots[2]-Boundary.knots[1])/ncongrid;
	aknots<- as.vector(Boundary.knots[1]+c(0:ncongrid)*cfac); 
	cmat<-rssbasis(aknots,knots,origin,intercept); 	
	nrc<-dim(cmat)[1];
	zmat<-cmat[2:nrc,]-cmat[1:(nrc-1),];
	qpfit<-solve.QP(Dmat=qmat,dvec=rmat,Amat=t(zmat),meq=0,factorized=F);
	betahat<-qpfit$solution;
	yhat<-X%*%betahat;
	msr<-mean((yhat-y)^2); 
	iact<-qpfit$iact; nact<-sum(iact>0);
	if(nact<1) {
	hat<-X%*%solve(XTX+(10^alpha)*D,tol=1.e-16)%*%t(X); Z<-1;
	}
	else {		
		iact<-iact[iact>0]; 
		Cact<-zmat[iact,];
		if(nact<2) Cact<-matrix(Cact,1,cols(zmat)) 
		Z<-Null(t(Cact));
		ZT<-t(Z);
		hat<-X%*%Z%*%solve(ZT%*%(XTX+(10^alpha)*D)%*%Z,tol=1.e-16)%*%ZT%*%t(X);
	}
	enp<-sum(diag(hat));
	d1<-sum(diag(hat))/length(y);
	gcv<-msr/((1-penalty*d1)^2);
	if(penalty*d1>1) gcv<-1.e12; 
 
	return(list(x=x,knots=knots,origin=origin,alpha=alpha,
		intercept=intercept, coef=betahat, fitted.values=yhat,gcv=gcv,
		penalty=penalty,Z=Z,enp=enp))
}

# Fit Monotone Penalized Regression Spline with smoothing parameter chosen by GCV(p)
# Hmin and Hmax are user-supplied brackets on the log10(optimal smoothing parameter) 
# Uses golden section search; hmin and hmax must bracket a minimum
gcv.rssM<-function(xvals,yvals,hmin,hmax,nknots=NA,Boundary.knots=NA,intercept=T,penalty=1,ncongrid=NA,tol=0.01,draw=T) {
	if(is.na(Boundary.knots)) Boundary.knots<-range(xvals);
	if(is.na(nknots)) nknots<-20;
	if(is.na(ncongrid)) ncongrid<-50; 

	ex<-order(xvals); xvals<-xvals[ex]; yvals<-yvals[ex]; 

	C<-0.61803399; R<-1-C;
	hvals<-hmin+0.1*c(0:10)*(hmax-hmin);
	Gvals<-matrix(0,11,1);
	for (i in 1:11) {
		Gvals[i]<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,
			alpha=hvals[i],penalty=penalty,ncongrid=ncongrid)$gcv;
		print(c(i,hvals[i],Gvals[i]));
        }
	imin<-minindc(Gvals);
	if((imin<=1)||(imin>=11)) {hopt<-hvals[imin]} else
	{ 
	ax<-hvals[imin-1]; bx<-hvals[imin]; cx<-hvals[imin+1]; 
	x0<-ax; x3<-cx; 
	if(abs(cx-bx)>abs(bx-ax)) 
		{x1<-bx; x2<-bx+C*(cx-bx)} 
	else {x2<-bx; x1<-bx-C*(bx-ax)}
	f1<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x1,ncongrid=ncongrid)$gcv;	
	f2<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x2,ncongrid=ncongrid)$gcv;
	while(abs(x3-x0)>tol*(abs(x1)+abs(x2))) {
		if(f2<f1) {x0<-x1; x1<-x2; x2<-R*x1+C*x3; f1<-f2; 
			f2<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x2,ncongrid=ncongrid)$gcv;}
		else {x3<-x2; x2<-x1; x1<-R*x2+C*x0; f2<-f1;
			f1<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=x1,ncongrid=ncongrid)$gcv;}
	}		
	if(f1<f2) {hopt<-x1; fmin<-f1} else {hopt<-x2; fmin<-f2}
	print(c(hopt,fmin));
	}

	bestfit<-rssfitM(xvals,yvals,nknots,Boundary.knots,intercept=intercept,penalty=penalty,alpha=hopt,ncongrid=ncongrid);	

	if (draw==T) {
	win.graph(); par(mfrow=c(2,1));
	plot(hvals,log10(Gvals),type="l",xlab="log10(alpha)",
		ylab="GCV");
	title("Bandwidth selection curve");
	plot(xvals,yvals); points(xvals,bestfit$fitted.values,type="l");
	title("Data and Monotone P-spline fit");
	}

	return(list(hvals=hvals,Gvals=Gvals,hopt=hopt,x=xvals,fitted.values=bestfit$fitted.values,
		penalty=penalty,gcv=bestfit$gcv,model=bestfit))
       				
}

bcgrad<-function(tvals,x,hder=NA) {
	xrange<-range(x)[2]-range(x)[1]; xmin<-range(x)[1];
	x<-(x-range(x)[1])/xrange;
	if(is.na(hder)) {hder<-(tvals[2]-tvals[1])}

	betahat<-locregG(tvals,x,tvals,hder,4); 
	xsmu<-betahat[,1]; xdot1<-betahat[,2];

	xfit<-xsmu;
	nt<-length(tvals); 
	tsim<-0.5*(tvals[2:nt]+tvals[1:(nt-1)]);
	betahat<-locregG(tvals,xfit,tsim,hder,4); 
	dxhat<-betahat[,2];

	f<-splinefun(tvals,xfit);
	dx<-0.05*(tvals[2]-tvals[1]);
	xplus<-f(tsim+dx); xminus<-f(tsim-dx);
	dxsim<-(xplus-xminus)/(2*dx);

	dxhat<-dxhat[5:(nt-4)];
	dxsim<-dxsim[5:(nt-4)]; 

	e<-order(dxhat);
	xvals<-dxhat[e]; yvals<-dxsim[e];

	fit<-gcv.rssM(xvals,yvals,hmin=-8,hmax=5,nknots=20,Boundary.knots=range(dxhat),
		intercept=T,penalty=1,ncongrid=50,tol=0.01,draw=F);
	
	f<-splinefun(xvals,fit$fitted.values,method="natural"); 
	yhat<-f(xdot1); 
	return(list(raw=xrange*xdot1,unbiased=xrange*yhat,smooth=xmin+xrange*xfit,hder=hder))
}
 
rssadd2<-function(x,y,nknots1=20,nknots2=20,Boundary.knots1,Boundary.knots2,intercept1=T,intercept2=T,
	penalty=1, alpha,signs=c(1,1)) {
#	Additive model, y=f_1(x[,1])+f_2(x[,2]) + e, cubic penalized regression splines for f1 & f2
# 	Input variable 'signs' indicates sign constraints: 
#		f_i >= 0 or <=0 according to whether signs[i] is +ive or -ive. 
# 	Uses solve.QP to apply sign constraints at the knots
#
	x1<-x[,1]; x2<-x[,2]; signs<-signs/abs(signs); 
	origin1<-Boundary.knots1[1];
	fac1<-(Boundary.knots1[2]-Boundary.knots1[1])/(nknots1+1);
	knots1<-origin1+c(1:nknots1)*fac1;
	x1mat<-rssbasis(x1,knots1,origin1,intercept1); 
	origin2<-Boundary.knots2[1];
	fac2<-(Boundary.knots2[2]-Boundary.knots2[1])/(nknots2+1);
	knots2<-origin2+c(1:nknots2)*fac2;
	x2mat<-rssbasis(x2,knots2,origin2,intercept2); 
	xmat<-x1mat%~%x2mat;
	XTX<-t(xmat)%*%xmat;

	d1<-rep(1,dim(x1mat)[2]); 
	d1[1:3]<-0; if(intercept1) d1[4]<-0;
	d2<-rep(1,dim(x2mat)[2]); 
	d2[1:3]<-0; if(intercept2) d2[4]<-0;
	D<-c((10^alpha[1])*d1,(10^alpha[2])*d2); D<-diag(D);

	# qprog matrices
	qmat<-XTX+D; rmat<-as.vector(t(xmat)%*%y);  	

	# constrain the sign of each term #
	aknots<-c(Boundary.knots1[1],knots1,Boundary.knots1[2]); 
	cmat1<- signs[1]*rssbasis(aknots,knots1,origin1,intercept1); 

	aknots<-c(Boundary.knots2[1],knots2,Boundary.knots2[2]); 
	cmat2<- signs[2]*rssbasis(aknots,knots2,origin2,intercept2); 

	cmat<-zeros(rows(cmat1)+rows(cmat2),cols(cmat1)+cols(cmat2));
	cmat[1:rows(cmat1),1:cols(cmat1)]<-cmat1;
	cmat[(1+rows(cmat1)):rows(cmat),(1+cols(cmat1)):cols(cmat)]<-cmat2;

	qpfit<-solve.QP(Dmat=qmat,dvec=rmat,Amat=t(cmat),meq=0,factorized=F);
	betahat<-qpfit$solution;
	yhat<-xmat%*%betahat;
	yhat1<-x1mat%*%betahat[1:cols(x1mat)]; 
	yhat2<-x2mat%*%betahat[(1+cols(x1mat)):(cols(x1mat)+cols(x2mat))];
	msr<-mean((yhat-y)^2);


	iact<-qpfit$iact; nact<-sum(iact>0);
	if(nact<1) {
	hat<-xmat%*%solve(XTX+D,tol=1.e-16)%*%t(xmat); Z<-1;
	}
	else {		
		iact<-iact[iact>0]; 
		Cact<-cmat[iact,];
		if(nact<2) Cact<-matrix(Cact,1,cols(cmat)) 
		Z<-Null(t(Cact)); 
		ZT<-t(Z);
		hat<-xmat%*%Z%*%solve(ZT%*%(XTX+D)%*%Z,tol=1.e-16)%*%ZT%*%t(xmat);
	}
	enp<-sum(diag(hat));
	d1<-sum(diag(hat))/length(y);
	if(penalty*d1>1) {gcv<-1.0e12} else {gcv<-msr/((1-penalty*d1)^2)};
 
	return(list(x=x,knots1=knots1,knots2=knots2,origin1=origin1, origin2=origin2,
		alpha=alpha, 
		intercept1=intercept1, intercept2=intercept2,coef=betahat,fitted.values=yhat1%~%yhat2,gcv=gcv,
		penalty=penalty,enp=enp))
}

gcvscore.rssadd2<-function(par,x,y,nknots1=20,nknots2=20,Boundary.knots1,Boundary.knots2,
	intercept1=T,intercept2=T,penalty=1,signs=c(1,1)) 
{
	fit<-rssadd2(x,y,nknots1=nknots1,nknots2=nknots2,Boundary.knots1=Boundary.knots1, 
	Boundary.knots2=Boundary.knots2,intercept1=intercept1,intercept2=intercept2,
	penalty=penalty,alpha=par,signs=signs)
	return(log10(fit$gcv))
}

gcv.rssadd2<-function(x,y,nknots1=20,nknots2=20,Boundary.knots1=NA,Boundary.knots2=NA,
	intercept1=T,intercept2=T,penalty=1,signs=c(1,1),par=c(0,0)) 
{
	if(is.na(Boundary.knots1)) Boundary.knots1<-range(x[,1]);
	if(is.na(Boundary.knots2)) Boundary.knots2<-range(x[,2]);
	bestfit<-optim(par,gcvscore.rssadd2,gr=NULL, method="Nelder-Mead", lower=-Inf, upper=+Inf,
	control=list(maxit=1000,trace=3), hessian=F,
	x,y,nknots1,nknots2,Boundary.knots1,
	Boundary.knots2,intercept1=intercept1,intercept2=intercept2,penalty,signs);
	hopt<-bestfit$par; 
	bestfit<-optim(hopt,gcvscore.rssadd2,gr=NULL, method="Nelder-Mead", lower=-Inf, upper=+Inf,
	control=list(maxit=1000,trace=3), hessian=F,
	x,y,nknots1,nknots2,Boundary.knots1,
	Boundary.knots2,intercept1,intercept2,penalty,signs);
	hopt<-bestfit$par; 

	rssfit<-rssadd2(x,y,nknots1,nknots2,Boundary.knots1,Boundary.knots2,intercept1=intercept1,
	intercept2=intercept2,penalty=penalty,alpha=bestfit$par,signs=signs);

	return(list(hopt=bestfit$par,gcv=rssfit$gcv, fitted.values=rssfit$fitted.values,model=rssfit));
}

predict.rssadd2<-function(model, newx) {
	x1<-newx[,1]; x2<-newx[,2];  
	x1mat<-rssbasis(x1,model$knots1,model$origin1,model$intercept1); 
	x2mat<-rssbasis(x2,model$knots2,model$origin2,model$intercept2); 
	betahat<-model$coef; 
	yhat1<-x1mat%*%betahat[1:cols(x1mat)]; 
	yhat2<-x2mat%*%betahat[(1+cols(x1mat)):(cols(x1mat)+cols(x2mat))];
	return(list(x=newx,y=yhat1%~%yhat2));
}





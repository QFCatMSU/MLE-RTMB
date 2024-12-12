#  PDF for normal with value of 2:  1/(sqrt(2*pi)*sig)*exp(-(2-mu)^2/(2*sig^2))
# exercise of finding joint pdf and demo dnorm, pnorm, qnorm, rnorm

#Set some values
xval=2.5
mu=2.5
var=3.5


#Illustrating finding density for normal
dnorm(x=xval,mean=mu,sd=sqrt(var));

#demonstration that dnorm just calculates the normal pdf

#write my own dnorm function and use it
mydnorm<-function(x,mu,var){
  1/sqrt(2*3.14156*var)*exp(-(xval-mu)^2/(2*var));
}

mydnorm(xval,mu,var)

#demo of calculating joint pdf of normal sample
#I did not bother with function given dnorm is vectorized and I used prod
xvec=c(10.72,7.23,10.07,8.62,8.55);
prod(dnorm(xvec,mean=10,sd=sqrt(2)))

# or the hard way:->
den<-vector("numeric",length=length(xvec))

for (i in 1:length(xvec)){
  den[i]<-dnorm(xvec[i],mean=10,sd=sqrt(2))
}
den[1]*den[2]*den[3]*den[4]*den[5]

# making the hard way a function just to illustrate this kind of coding
jnt_den=function(xvec,mu,var){
prod_den=1;
for (i in 1:length(xvec)){
  prod_den = prod_den*dnorm(xvec[i],mean=mu,sd=sqrt(var));
}
prod_den;
}

jnt_den(xvec,mu=10,var=2)

#we used dnorm, families of R probability functions include "p" and "q" versions
#By default pnorm produces the CDF
pnorm(2,mean=2, sd=sqrt(2))
curve(pnorm(x,mean=2,sd=sqrt(2)), from=-5, to = 5)

#using pnorm to find the probability of being in an interval
pnorm(xval+1,mean=3,sqrt(2))-pnorm(xval,mean=3,sqrt(2))

#qnorm gets quantiles
# i.e., the value for which Pr(value<x)=p
qnorm(p=0.5,mean=3,sqrt(2))


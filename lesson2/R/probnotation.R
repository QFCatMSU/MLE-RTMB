#Simple demo of prob notation
library(RTMB)

#The explicit way we have and will do it in class
tst1 = function(x,y){
  NLL= -dnorm(x,0,1,log=TRUE)
  NLL = NLL - dnorm(y,0,1,log=TRUE)
  NLL
}

#probalistic notation way
tst2<-function(x,y) {
x %~% dnorm(0,1)
y %~% dnorm(0,1)
}

tst1(1,2)
tst2(1,2)


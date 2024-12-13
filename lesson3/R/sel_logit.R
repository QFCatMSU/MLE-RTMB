logit=function(x){
  log(x/(1-x))
}

invlogit=function(x){
  1/(1+exp(-x))
}

# tr_adjust might be starting values of transformed adjusts
tradjust=rep(logit(.9),9)

# code like this would go inside NLL function
adjust = invlogit(tradjust);
sel=rep(NA,10);
sel[10]=1;
for(i in 1:9){
  sel[10-i]=adjust[10-i]*sel[10-i+1]
}
sel



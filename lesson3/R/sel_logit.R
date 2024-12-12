logit=function(x){
  log(x/(1-x))
}

inv_logit=function(x){
  1/(1+exp(-x))
}

# tr_adjust might be starting values of transformed adjusts
tr_adjust=rep(logit(.9),9)

# code like this would go inside NLL function
adjust = inv_logit(tr_adjust);
sel=rep(NA,10);
sel[10]=1;
for(i in 1:9){
  sel[10-i]=adjust[10-i]*sel[10-i+1]
}
sel



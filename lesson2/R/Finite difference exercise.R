
#Finite difference derivatives

#We will find derivative of function wrt to x

# function we will find derivative of
myfun = function(x,a,b){
  a+b*x+sin(x);
}

#centered finite difference deriv function
fin_diff = function(x,a,b,relh){
  h=relh*x;
  (myfun(x+h/2,a,b)-myfun(x-h/2,a,b))/h;
}

fin_diff(x=1,a=2,b=0.5,relh=0.001);
fin_diff(x=2,a=1,b=1,relh=0.001);
#compare with analytical answer b+cos(x)
0.5+cos(1);
1+cos(2);

#what about bigger relh?
fin_diff(x=1,a=2,b=0.5,relh=0.1);
fin_diff(x=2,a=1,b=1,relh=0.1);

#second derivative function
sec_der<-function(x,a,b,relh){
  h=relh*x;
  (fin_diff(x+h/2,a,b,relh)-fin_diff(x-h/2,a,b,relh))/h
}

sec_der(x=1,a=2,b=0.5,relh=0.01)

#analytical check -sin(x)
-sin(1)

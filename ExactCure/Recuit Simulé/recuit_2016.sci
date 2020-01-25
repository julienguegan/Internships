function y=f(x)
    y=x.^2-cos(2*%pi*x)+ones(x)
endfunction
//
x=-20+30*rand();  // point initial
Niter=2000;alpha=0.5;Ytot=[]
for i=1:Niter
  y1=f(x);Ytot=[Ytot,y1];  
  xtilde=x+(-alpha+2*alpha*rand())  
  y2=f(xtilde)
  p=exp(-(y2-y1)/(1/log(i+1)));
  if (rand()<p) then 
    x=xtilde;  // acceptation  
  end              
end
disp('valeur finale obtenue pour x:')
disp(x)
figure(1)
clf()
plot2d(Ytot) // tracé de (Niter,f(xn))
//
figure(2)    // tracé de la fonction à minimiser
xp=-20:0.01:10;
yp=f(xp);
plot2d(xp,yp)

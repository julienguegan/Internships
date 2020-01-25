// Methode d'optimisation sans gradient
//
//   Methode de recuit simule
//----------------------------------------
// definition de la fonction a minimiser
//
function y=f(x)  // the function to optimize
n=prod(size(x));
y=n+sum(x.^2-cos(2*%pi*x));
endfunction
//
///////////////////////////////////////////////////////////////
//
n=evstr(x_dialog('nombre de parametres de Rastrigin','2'));
Nmax=evstr(x_dialog('nombre maximal d''évaluations','2000'));
//
disp('Niter  Temperature  sigma');
u=rand(1,n);
xmin=-5*ones(1,n);
xmax=5*ones(1,n);
x0=xmin+(xmax-xmin).*u;  // random initialisation 
Neval=[];Fbest=[];
y0 = f(x0);
x = x0; y = y0; alpha = 0.92; T = 1; s = 0.2;icount=1;
Neval=[];
Fbest=[];
for k = 1 : 20
   T = alpha * T; 
   accept = 0;
   for l = 1 : Nmax/20
      x1 = x0 + s*(0.5*ones(x0) - rand (x0));  // uniform
//     x1= x0+s*rand(x0,'normal');              // gaussian
      y1 = f(x1);
      icount=icount+1;
      dy = y1 - y0;
      if rand() < exp (- dy / T) then
         x0 = x1; y0 = y1; accept = accept + 1;
      end
      if y0 < y then x = x0; y = y0; end
   Neval=[Neval;icount];
   Fbest=[Fbest;f(x0)];
    end
  if accept < 25/100*Nmax/20 then s = s / 2; end
  if accept > 75/100*Nmax/20 then s = 2 * s; end
  disp([icount,T,s]);
end
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
disp('minimum obtained after SA:');disp(x);
disp('corresponding value by f:');disp(y);
disp('function evaluation number by SA:');disp(Neval($));
xset('window',0)
clf();
plot2d(Neval,Fbest)
//
/////////////////////////////////////////////////////////////

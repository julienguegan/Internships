clear
exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\affichage.sce')
function alphak = linesearch(f,d,x,rho,c)
// d: The search direction
// rho :- The backtrack step between (0,1) usually 1/2
// c: parameter between 0 and 1 , usually 10^{-4}
alphak = 1;
fk = f(x);
gk = numderivative(f,x)
xx = x;
disp(x);disp(alphak);disp(d)
x = x + alphak*d;
fk1 = f(x);
while fk1 > fk + c*alphak*(gk'*d),
  alphak = alphak*rho;
  x = xx + alphak*d;
  fk1 = f(x);
end
endfunction

function x = newton( f, x, tol, maxit)
fx = f(x)
[gx, Hx] = numderivative(f,x)
it = 1;
while( it < maxit & norm(gx,2) > tol )
   s  = -Hx\gx;
   alphak = linesearch(f,s,x,0.5,1e-4);
   x = x + alphak*s;
   it = it+1;
   [fx,gx,Hx] = feval(f,x);
end
endfunction

function z = rozenbrock(x)
    z =  10*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction

affiche(rozenbrock,2,-2,2,-2,'contour')
f = rozenbrock
x0 = [-1;1]
[gx, Hx] = numderivative(f,x0)
s = Hx\gx/*
tol = 10^-3
itermax = 100

x = newton(f,x0,tol,itermax)
disp(x);disp(n)

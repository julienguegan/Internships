exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\gradient\recherche lineaire.sce',-1)
exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\affichage.sce')

function [x,n] = newton(f,tol,itermax,x0)
    x = x0
   [Jac, Hess] = numderivative(f, x0)
    direction = -Jac
    n = 1
    xold = 0
    while (tol<(norm(xold-x)))&(n<itermax) then
        n = n+1 
        xold = x
        alpha = linesearch(f,x,direction,Jac)
        x = x+alpha*direction 
       [Jac, Hess] = numderivative(f, x)
        direction = inv(Hess)*Jac
        plot(x(1),x(2),'k.')
    end
endfunction

function z = rozenbrock(x)
    z =  10*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction
affiche(rozenbrock,2,-2,2,-2,'contour')
f = rozenbrock
tol = 10^-3
itermax = 100
x0 = [-1;1]
[x n] = newton(f,tol,itermax,x0)
disp(x);disp(n)

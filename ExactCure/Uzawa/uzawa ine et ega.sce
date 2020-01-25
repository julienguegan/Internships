clear

exec('C:\Users\Julien Guégan\Desktop\PFE\fonctions test.sce',-1)
exec('C:\Users\Julien Guégan\Desktop\PFE\affichage.sce',-1)
/*function alpha = backtracking(f,x,d,grad)
    alpha = 1
    w = 10^-4
    cpt = 0
    disp(f(x+alpha*d))
    disp('heoh')
    fx = f(x) //une seule evaluation
    while(f(x+alpha*d)>(fx+alpha*w*(grad'*d)))
        alpha = alpha/2 
        cpt = cpt+1
    end
endfunction
function y=Lx(x),y=lagrangien(f,h,g,lambdan,mun,x),endfunction
        function y=Lm(m),y=-lagrangien(f,h,g,lambdan,m,xnp1),endfunction
        function y=Ll(l),y=-lagrangien(f,h,g,l,mun,xnp1),endfunction
        rhox = backtracking(Lx,xn,-dLx,dLx)
        rhol = backtracking(Ll,lambdan,h(xnp1),h(xnp1))
        rhom = backtracking(Lm,munp1,g(xnp1),g(xnp1))
*/
function L = deriveelagrangien(f,h,g,lambda,mu)
    S = 0
    H = h
    G = g
    for i = 1:size(H,2)
        S = S + lambda(i).*H(i)
    end
    for j = 1:size(G,2)
        S = S + mu(j).*G(j)
    end
    L = f+S
endfunction

function [L,dLx] = lagrangien(f,h,g,lambda,mu,x)
    S = 0
    dS = 0
    H = h(x)
    G = g(x)
    F = f(x)
    df = numderivative(f,x)'
    dh = numderivative(h,x)'
    dg = numderivative(g,x)'
    for i = 1:size(H,2)
        S = S + lambda(i).*H(i)
    end
        dS = dS + lambda.*dh
        dS = dS + mu.*dg
    for j = 1:size(G,2)
        S = S + mu(j).*G(j)
    end
    L = F + S
    dLx = df + dS
endfunction


function [xn,n,stock] = uzawa(f,tol,itermax,x0,h,g)
    xn = x0
    stock(:,1) = x0
    for i = 1:length(h(x0))
        lambdan(i) = 1
    end
    for j = 1:length(g(x0))
        mun(j) = 1
    end
    n = 1
    xnp1 = x0

   while (((n<itermax)&(norm(xn-xnp1)>tol))| (n == 1)) then
        xn = xnp1
        stock(:,n)=xnp1
        [L dLx] = lagrangien(f,h,g,lambdan,mun,xn)
        
        rhog = 0.2
        rhoh = 0.6
        rhom = 0.6
        xnp1 = xn - rhox*dLx
        lambdanp1 = lambdan + rhoh*h(xnp1)
        munp1 = max(0, mun + rhog*g(xnp1))
        
        lambdan = lambdanp1
        mun = munp1
        /*plot([xnp1(1) xn(1)],[xnp1(2) xn(2)],'k-')
        plot(xn(1),xn(2),'k.','markersize',3)*/
        n = n+1
    end

endfunction

function z = cout(x)
    z = 0.05*x(1)^4+ 0.1*x(2)^4 + 10*x(1)*x(2)+20*x(1)//x(1)^3+x(2)
endfunction
function eq = egalite(x)
    eq = (x(1)+2)^2+(x(2)+2)^2/9-1
endfunction
function in = inegalite(x)
    in = (-x(1)/2+x(2)+1)
endfunction

clf()
affiche(cout,10,-10,10,-10,'contour')

exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\cercle.sce',-1)
ellipse(-2,-2,1,3,0)
x=-3.5:0.1:-0.5
plot(x,x/2-1,'b-')
legend(["$h = 0$";"$g < 0$"] ,-1, %f)

 
tol = 0.0001
itermax = 1000
x0 = [0;6]
fonction = cout

tic()
[sol n stock] = uzawa(fonction,tol,itermax,x0,egalite,inegalite)
time = toc()

disp('le minimum est x = ')
disp(sol)
disp(' au bout de '+string(n)+' iterations, '+string(time)+' secondes')

tracer = stock'
for i=1:size(tracer,1)-1
    plot(tracer(i,1),tracer(i,2),'k-')//'k.','markersize',3)
    plot([tracer(i,1) tracer(i+1,1)],[tracer(i,2) tracer(i+1,2)],'k-')
end
plot(sol(1),sol(2),'r.')
xstring(sol(1),-1,'$min$')
gce().font_size=3
plot(x0(1),x0(2),'k.')
xstring(x0(1),x0(2),'$x_0$')
gce().font_size=3

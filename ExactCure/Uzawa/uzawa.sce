clear
exec('C:\Users\Julien Guégan\Desktop\PFE\fonctions test.sce',-1)
exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\Uzawa\recherche lineaire.sce',-1)
exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\Gradient\Gradient.sce',-1)
function alpha = backtracking(f,x,d,grad)
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
        stock(:,n)= xn

        [L dLx] = lagrangien(f,h,g,lambdan,mun,xn)
        function y=Lx(x)
            y=lagrangien(f,h,g,lambdan,mun,x)
        endfunction
        /* function y=Lm(m),y=lagrangien(f,h,g,lambdan,m,xn),endfunction
        function y=Ll(l),y=lagrangien(f,h,g,l,mun,xn),endfunction*/
        rhox = 0.01
        rhol = 0.01
        rhom = 0.01
       // a = linesearch(Lx,xn,-dLx,dLx)
        xnp1 = xn - rhox*dLx
        lambdanp1 = lambdan + rhol*h(xnp1)
        munp1 = max(0, mun + rhom*g(xnp1)) 

        lambdan = lambdanp1
        mun = munp1

        plot([xnp1(1) xn(1)],[xnp1(2) xn(2)],'k-')
        //plot(xn(1),xn(2),'k.','markersize',2)

        n = n+1
    end

endfunction

function z = cout(x)
    z = x(1)^3+x(2)
endfunction
function eq = egalite(x)
    eq = x(1)^2+2*x(2)^2-1
endfunction
function in = inegalite(x)
    in = -x(1)+x(2)+0.5
endfunction

clf()
x = -2:0.1:2
for i=-2:0.1:2
    ind_i = i*10+21
    for j=-2:0.1:2
        u=[i j]
        ind_j = j*10+21
        z(ind_i,ind_j) = cout(u)
    end
end
clf()
g = gcf()
g.color_map = rainbowcolormap(64)
xset("fpf"," ")
contour2d(x,x,z,40)
xlabel('$x$','fontsize',4)
ylabel('$y$','fontsize',4)
exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\cercle.sce',-1)
ellipse(0,0,1,1/sqrt(2),0)
plot(x,x)

tol = 0.0001
itermax = 1000
x0 = [-1;1]
fonction = cout

tic()
[sol n stock] = uzawa(fonction,tol,itermax,x0,egalite,inegalite)
time = toc()

disp('le minimum est x = ')
disp(sol)
disp(' au bout de '+string(n)+' iterations, '+string(time)+' secondes')
/*
tracer = stock'
for i=1:size(tracer,1)-1
    plot(tracer(i,1),tracer(i,2),'k.','markersize',3)
    plot([tracer(i,1) tracer(i+1,1)],[tracer(i,2) tracer(i+1,2)],'k-')
end*/
plot(sol(1),sol(2),'r.')
xstring(sol(1),sol(2),'$min$')
gce().font_size=3

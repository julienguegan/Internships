clear all

exec('C:\Users\Julien Guégan\Desktop\PFE\affichage.sce',-1)
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
        dS = dS + lambda(i).*dh(i,:)
    end
        /*dS = dS + lambda.*dh
        dS = dS + mu.*dg*/
    for j = 1:size(G,2)
        S = S + mu(j).*G(j)
        dS = dS + mu(j).*dg(:,j)
    end
    
    L = F + S
    dLx = df + dS
endfunction


function xn = uzawa(f,tol,x0,h,g)
    xn = x0
    for i = 1:length(h(x0))
        lambdan(i) = 1
    end
    for j = 1:length(g(x0))
        mun(j) = 1
    end
    n = 1
    cdtarret = %T
    while (cdtarret) then
        [L,dLx] = lagrangien(f,h,g,lambdan,mun,xn)
        rhox = 0.01
        xnp1 = xn - rhox*dLx
        rhol = 10
        lambdanp1 = lambdan + rhol*h(xnp1)
        rhom = 1
        munp1 = max(0, mun + rhom*g(xnp1)) 
        plot([xnp1(1) xn(1)],[xnp1(2) xn(2)],'k-')
        cdtarret = norm(xn-xnp1)>tol
        
        lambdan = lambdanp1
        mun = munp1
        xn = xnp1
        n = n+1
    end

endfunction

function z = cout(x)
    z = 0.05*x(1)^4+ 0.1*x(2)^4 + 10*x(1)*x(2)
endfunction
function h = egalite(x)
    h = x(2)-x(1)
endfunction
function g = inegalite(x)
    g1 = x(1)-8
    g2 = -x(1)-8
    g3 = -x(2)-8
    g4 = x(2)-8
    g = [g1;g2;g3;g4]
endfunction

clf()
affiche(cout,10,-10,10,-10,'contour')
xlabel('$x_1$','fontsize',4)
ylabel('$x_2$','fontsize',4)
//exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\cercle.sce',-1);ellipse(0,0,5,5,0)
x = [-10:0.1:10]'
plot(x,8*ones(length(x),1),'b-.');plot(x,-8*ones(length(x),1),'b-.');plot(8*ones(length(x),1),x,'b-.');plot(-8*ones(length(x),1),x,'b-.')
//legend(["égalité";"inégalité"] ,-1)
plot(x,x)

tol = 0.0001
x0 = [-4;-4]
fonction = cout
plot(x0(1),x0(2),'k.')
sol = uzawa(fonction,tol,x0,egalite,inegalite)
disp('le minimum est x = ')
disp(sol)
plot(sol(1),sol(2),'k.')

clear all

exec('C:\Users\Julien Guégan\Desktop\PFE\affichage.sce',-1)
exec('C:\Users\Julien Guégan\Desktop\PFE\fonctions test.sce',-1)
function alpha = backtracking(f,x,d,grad)//pour newton mais pas pour quasi et CG
    alpha = 10
    w = 0.1
    cpt = 0
    while(f(x+alpha*d)>(f(x)+alpha*w*(grad'*d)))
        alpha = alpha/2 
        cpt = cpt+1
    end
endfunction


function [L,dLx] = lagrangien(f,g,mu,x)
    S = 0
    dS = 0
    G = g(x)
    F = f(x)
    df = numderivative(f,x)
    dg = numderivative(g,x)
    for j = 1:size(G,1)
        S = S + mu(j).*G(j)
        dS = dS + mu(j).*dg(j,:)
    end
    L = F + S
    dLx = df' + dS'
endfunction


function xn = uzawa(f,tol,x0,g)
    xn = x0
    for j = 1:length(g(x0))
        mun(j) = 1
    end
    n = 1
    cdtarret = %T
    while (cdtarret) then
        [L,dLx] = lagrangien(f,g,mun,xn)
        rhox = 0.05
        xnp1 = xn - rhox*dLx
        rhom = 1
        munp1 = max(0, mun + rhom*g(xnp1)) 
        plot([xnp1(1) xn(1)],[xnp1(2) xn(2)],'k-')
        cdtarret = norm(xn-xnp1)>tol

        mun = munp1
        xn = xnp1
        n = n+1
    end

endfunction

function z = cout(x)
    z = 0.05*x(1)^4+ 0.1*x(2)^4 + 10*x(1)*x(2)+20*x(1)
endfunction

function g = inegalite(x) //dlf g(x)<0
    g1 = -x(1)-2
    g2 = x(1)-2
    g3 = -x(2)+2
    g4 = x(2)-6
    g = [g1;g2;g3;g4]
endfunction

clf()
affiche(quadratique,4.9,-4.9,9.1,-3,'contour')
xlabel('$x_1$','fontsize',4)
ylabel('$x_2$','fontsize',4)
x = [-2.5:0.1:2.5]'
x2 = [0.5:0.1:7.5]'
plot(x,6*ones(length(x),1),'b-.');plot(x,2*ones(length(x),1),'b-.');plot(-2*ones(length(x2),1),x2,'b-.');plot(2*ones(length(x2),1),x2,'b-.')
//legend(["égalité";"inégalité"] ,-1)


tol = 0.001
fonction = quadratique
for x01 = -2:2 //samplé plusieurs initial guess
    for x02 = 2:6
        x0 = [x01;x02]
/*plot(x0(1),x0(2),'k.')
xstring(x0(1),x0(2),'$x_0$')
gce().font_size=3*/
sol = uzawa(fonction,tol,x0,inegalite)
//disp('le minimum est x = ')
//disp(sol)
plot(sol(1),sol(2),'k.')
end
end
title('$\rho_x = 0.05;\rho_{\mu} = 1$')

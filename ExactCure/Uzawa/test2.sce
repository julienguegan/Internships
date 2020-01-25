clear all

exec('C:\Users\Julien Guégan\Desktop\PFE\affichage.sce',-1)
function alpha = backtracking(f,x,d,grad)
    alpha = 1
    w = 10^-4
    cpt = 0
    fx = f(x) //une seule evaluation
    while(f(x+alpha*d)>(fx+alpha*w*(grad'*d)))
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
        function y=Lx(x)
            y=lagrangien(f,g,mun,x)
        endfunction
        function y=Ll(l)
            y=lagrangien(f,h,g,l,mun,x)
        endfunction
        //disp(Lx(xn-0.01*dLx));disp((Lx(xn)-(10^-4)*0.01*(dLx'*dLx)))
        //rhox = backtracking(Lx,xn,-dLx,dLx)
        xnp1 = xn - rhox*dLx
        //rhom = backtracking(Ll,xn,g(xnp1),g(xnp1))
        munp1 = max(0, mun + rho*g(xnp1)) 
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

function g = inegalite(x)
    g1 = x(1)-6
    g2 = -x(1)-6
    g3 = -x(2)-6
    g4 = x(2)-4
    g = [g1;g2;g3;g4]
endfunction

clf()
affiche(cout,10,-10,10,-10,'contour')
xlabel('$x_1$','fontsize',4)
ylabel('$x_2$','fontsize',4)
//exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\cercle.sce',-1);ellipse(0,0,5,5,0)
x = [-10:0.1:10]'
plot(x,-6*ones(length(x),1),'b-.');plot(x,4*ones(length(x),1),'b-.');plot(-6*ones(length(x),1),x,'b-.');plot(6*ones(length(x),1),x,'b-.')

tol = 0.0001
x0 = [-2;-6]
fonction = cout
plot(x0(1),x0(2),'k.')
xstring(x0(1),x0(2),'$x_0$')
gce().font_size=3
plot(-7.2630161,5.6626344,'b.')
plot(5.3864445,-5.1256305,'b.')
sol = uzawa(fonction,tol,x0,inegalite)
disp('le minimum est x = ')
disp(sol)
plot(sol(1),sol(2),'k.')

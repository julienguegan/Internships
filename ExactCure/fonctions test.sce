
// ******** fonctions 1D *********//
function z = carre(x)
    z = x^2
endfunction

function z = oscillante(x)
    z = (x*sin(x))^2+0.01*abs(x)*(cos(1000*x)+1)
endfunction

function z = chebfun(x);
    domaine = 0:6
    z = 20*cos(x)*sin(exp(x));
endfunction

function z = rastr(x);
    domaine = -5:5
    z = 10+x^2-10*cos(%pi*2*x);
endfunction

// ******** courbes monstres  *********//

function [z,domaine] = weierstrass(x)
    a = 1/3
    b = 5
    n = 50
    z = 0
    for i = 0:n
        z = z+a^i*cos(b^i*%pi/2*x)
    end
endfunction

function z = blancmanger(x)
    n = 100
    z = 0
    for i = 0:n
        z = z+abs(2^i*x-round(2^i*x))/(2^i)
    end
endfunction

function z = riemann(x)
    n = 100
    z = 0
    for i = 1:n
        z = z+sin(i^2*x)/(i^2)
    end
endfunction
// ******** fonctions diverses *********//

function z = quadratique(x)
    z = 4*x(1).^2+x(2).^2+3*x(1).*x(2)+x(1)+x(2)
endfunction

function z = exponentielle(x)
    z =  5*exp(2*x(1)^2)+2*x(2)^2+x(1)*x(2)
endfunction

function z = pointselle(x)
    z =  (x(1)^2)-(x(2)^2)
endfunction

function z = prodcos(x)
    z =  sin(0.5*x(1)^2-0.25*x(2)^2+3)*cos(2*x(1)+1-exp(x(2)))
endfunction

function z = sphere(x)
    z = 0
    for i = 1:length(x)
        z =  z + x(i)^2
    end
endfunction

// ******** fonctions TEST optimisation 2D *********//
function z = rozenbrock(x)
    z =  10*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction

function z = himmelblau(x)
    z =  (x(1)^2+x(2)-11)^2 + (x(1)+x(2)^2-7)^2;
endfunction

function [z,domaine] = rastrigin(x)
    domaine = -5.12:0.2:5.12
    n = 2
    A = 5
    z =  A*n+x(1)^2-A*cos(2*%pi*x(1))+x(2)^2-A*cos(2*%pi*x(2));
endfunction

function [z,domaine] = ackley(x)
    domaine = -5:0.2:5
    z = -20*exp(-0.2*sqrt(0.5*(x(1)^2+x(2)^2)))-exp(0.5*(cos(2*%pi*x(1))+cos(2*%pi*x(2))))+%e+20
endfunction

function [z,domaine] = beale(x)
    domaine = -4.5:0.2:4.5
    z = (1.5-x(1)+x(1)*x(2))^2+(2.25-x(1)+x(1)*x(2)^2)^2+(2.625-x(1)+x(1)*x(1)*x(2)^3)^2
endfunction

function [z,domaine] = goldsteinprice(x)
    domaine = -2:0.1:2
    z = (1 + ((x(1)+ x(2) + 1).^2) * (19 - (14 * x(1)) + (3 * (x(1)^2)) - 14*x(2) + (6 .* x(1)*x(2)) + (3 * (x(2).^2)))) *(30 + ((2 * x(1)- 3 * x(2)).^2) .* (18 - 32 * x(1)+ 12 * (x(1)^2) + 48 * x(2) - (36 .* x(1)*x(2)) + (27 * (x(2)^2))) );
endfunction

function [z,domaine] = booth(x)
    domaine = -10:0.1:10
    z = (x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;
endfunction

function [z,domainex,domainey] = bukin(x)
    domainex = -15:0.5:-5
    domainey = -3:0.1:3
    z = 100*sqrt(abs(x(2)-0.01*x(1)^2))+0.01*abs(x(1)+10);
endfunction

function [z,domaine] = levi(x)
    domaine = -10:0.5:10
    z = sin(3*%pi*x(1))^2+(x(1)-1)^2*(1+sin(3*%pi*x(2))^2)+(x(2)-1)^2*(1+sin(2*%pi*x(2)));
endfunction

function [z,domaine] = easom(x)
    domaine = -100:100
    z = -cos(x(1))*cos(x(2))*exp(-((x(1)-%pi)^2+(x(2)-%pi)^2));
endfunction

function [z,domaine] = crossintray(x)
    domaine = -10:0.1:10
    z = -0.0001*(abs(sin(x(1)*sin(x(2))*exp(abs(100-sqrt(x(1)^2+x(2)^2)/%pi))))+1)^0.1;
endfunction

function z = keane(x)
    num = (sin(x(1) - x(2)) .^ 2) .* (sin(x(1) + x(2)) .^ 2); 
    den = sqrt(x(1) .^2 + x(2) .^2);
    z = - num ./ den;
end
/*
scf()
domaine=0:0.0005:1
for i=0:0.0005:1
    ind_i=i/0.0005+1
    z(ind_i) = himmelblau(i)
end
plot(domaine',z)*/
// ****  AFFICHAGE *****//
/*
fonction = himmelblau
domainex=-5:0.2:5
domainey=-5:0.2:5
for i=-5:0.2:5
    for j=-5:0.2:5
        u=[i j]
        ind_i=i*5+26
        ind_j=j*5+26
        z(ind_i,ind_j) = fonction(u)
    end
end
f = scf();
plot3d1(domainex,domainey,z)//flag=[0,1,0])
f.color_map = rainbowcolormap(134)
xlabel('$x$','fontsize',4)
ylabel('$y$','fontsize',4)
zlabel('$z$','fontsize',4)

g=scf()
xset("fpf"," ")
contour2d(domainex,domainey,z,40)
g.color_map = rainbowcolormap(70)
x0 = [-1 -1]
[fopt, xopt] = fminsearch(fonction, x0')
plot(fopt(1),fopt(2),'k.','markersize',4)
xstring(fopt(1),fopt(2),'$min$')
gce().font_size=4

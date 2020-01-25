clear

function y = rastrigin(x)
    n = 2
    A = 5
    y =  A*n+x(1)^2-A*cos(2*%pi*x(1))+x(2)^2-A*cos(2*%pi*x(2));
endfunction

function z = rozenbrock(x)
    z =  10*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction

// ********* calcul F *********//
f = rastrigin

xmax = 3;xmin = -3;ymax = 3;ymin = -3;
xstep = (xmax-xmin)/65
axe_x = xmin:xstep:xmax
ystep = (ymax-ymin)/65
axe_y = ymin:ystep:ymax

z=0
for i = xmin:xstep:xmax
    for j = ymin:ystep:ymax
        u = [i j]
        ind_i = (i/xstep)-(xmin/xstep)+1
        ind_j = (j/ystep)-(ymin/ystep)+1
        z(ind_i,ind_j) = f(u)
    end
end

clf()
contour2d(axe_x,axe_y,z,45)
a = gcf()
a.color_map = rainbowcolormap(64)
////////////////////////////////////

tol = 10^-4
itermax = 100
//points initiales
x1 = [-1.3 2.5]
x2 = [-1.4 1]
x3 = [-0.5 1.5]
x = [x1;x2;x3]
A = [f(x1);f(x2);f(x3)]
//initialisation du critere d'arret
xbar = (x3 + x2)./2
xbarold = [0 0]

k = 0
while (k<itermax )//& norm(xbar - xbarold)>tol)|k==1
    xbarold = xbar

    [fsort index] = gsort(A) //tri par ordre décroissant
    xmin = x(index(3),:)
    xmax = x(index(1),:)
    xbar = (xmin + x(index(2),:))./2

    //****  AFFICHAGE  ***//
    // clf()
    // xset("fpf"," ")
    a = gcf()
    a.color_map = rainbowcolormap(64)
    plot([x(1,1) x(2,1)],[x(1,2) x(2,2)],'--','thickness',1)
    plot([x(1,1) x(3,1)],[x(1,2) x(3,2)],'--','thickness',1)
    plot([x(3,1) x(2,1)],[x(3,2) x(2,2)],'--','thickness',1)
    plot(xbar(1),xbar(2),'o')
    
    //reflection
    xrefl = xbar + (xbar - xmax)

    if f(xrefl)<f(xmin) then //expansion
        xe = xbar +2*(xbar - xmax)
        if f(xe)<f(xrefl) then
            xmax = xe
        else
            xmax = xrefl
        end
    elseif f(xrefl)<f(xmax) & f(xrefl)>f(xmin)
        xmax = xrefl
    else //contraction
        xe = xbar - 1/2*(xbar - xmax)
        if f(xe)<f(xmin)  then
            xmax = xe
        else
            x(index(2),:) = xmax + 0.5*(x(index(2),:)-xmin)
        end
    end
    k = k+1
    x = [xmax;x(index(2),:);xmin]
    A = [f(xmax);f(x(index(2),:));f(xmin)]
    // xs2png(gcf(),"image"+string(k))
end

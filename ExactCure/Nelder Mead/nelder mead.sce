clear

function y = f(x)
    n = 2
    A = 5
    y =  A*n+x(1)^2-A*cos(2*%pi*x(1))+x(2)^2-A*cos(2*%pi*x(2));
endfunction


tol = 10^-4
itermax = 1000

// triangle equilateral de barycentre (x0 y0)  //
x0 = 1
y0 = 1
c = 2
x1 = [y0-sqrt(3)/6 , x0-c/2]
x2 = [y0-sqrt(3)/6 , x0+c/2]
x3 = [x0 , y0+sqrt(3)/3]
//                                            //


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

    //reflection
    xrefl = xbar + (xbar - xmax)

    if f(xrefl)<f(xmin) then //expansion
        xexp = xbar +2*(xbar - xmax)
        if f(xe)<f(xrefl) then
            xmax = xexp
        else
            xmax = xrefl
        end
    elseif f(xrefl)<f(xmax) & f(xrefl)>f(xmin)
        xmax = xrefl
    else //contraction
        xcon = xbar - 1/2*(xbar - xmax)
        if f(xcon)<f(xmin)  then
            xmax = xcon
        else
            x(index(2),:) = xmax + 0.5*(x(index(2),:)-xmin)
        end
    end
    k = k+1
    x = [xmax;x(index(2),:);xmin]
    A = [f(xmax);f(x(index(2),:));f(xmin)]
    disp(xbar)
end

disp(xbar')

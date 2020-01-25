exec('C:\Users\Julien Guégan\Desktop\PFE\algorithmes\affichage.sce',-1)


function z = cout(x)
    z = 0.05*x(1)^4+ 0.1*x(2)^4 + 10*x(1)*x(2)+20*x(1)
endfunction

function g = inegalite(x) //dlf g(x)<0
    g = -x(1)-2
endfunction

f = cout
g = inegalite
clf()
affiche(fonction,10,-10,10,-10,'contour')
xlabel('$x_1$','fontsize',4)
ylabel('$x_2$','fontsize',4)

tol = 0.001
x0 = [0 0] 

x = x0
cdtarret = %T
n = 1
L = length(g(x))
M = []
while (cdtarret) then
    n = n+1
    df = numderivative(f,x)
    G = g(x)
    for i = 1:L
        if (G(i) == 0) //contrainte active
            M = [M numderivative(g,x)]
        end
    end
    r = -(eye(n) - M*inv(M'*M)*M')*df
    while (r == 0)
        u = -inv(M'*M)*M'*df
        if (min(u) < 0)
            indice = find(u == min(u))
            M(:,ind) = []
        else
            return x;
        end
    end 
    alpha = 0.01//linesearch(f,x,-df,df)//le pas
    G = g(x + alpha*r)
    for i = 1:L
        while (G(i) < 0) //contrainte inactive
            alpha = alpha /2
        end
    end
    xnp1 = x + alpha*df
    cdtarret = norm(x-xnp1)>tol
    x = xnp1
    n = n+1
end



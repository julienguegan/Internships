clear


function y=frozenbrock(x)
    if size(x,1) < 2
        error('dimension must be greater one');
    end
    N = size(x,1)
    y = 100*sum((x(1:N-1).^2 - x(2:N)).^2) + sum((x(1:N-1)-1).^2);
endfunction

function  [xm,n] = evolutionstrategie(f,xm0,sm0,itermax,tol) // μ < λ
    lambda = 4+floor(3*log(N)); 
    mu = floor(lambda/2); 
    
    dim = length(xm0)
    tau = 1 //facteur d'apprentissage, petit = sigma plus petit
    n = 1 //nbr iteration

    xmp1 = xm0
    xm = xm0
    smp1 = sm0
    rand("normal")
    while (((n<itermax)&(norm(xm-xmp1)>tol))| (n == 1)) then

        xm = xmp1
        sm = smp1

        for i = 1:lambda //creation
            s(i,:) = sm .* exp(tau*rand(1,dim))
            x(i,:) = xm + s(i).*rand(1,dim)
        end

        Z = zeros(lambda)
        for i = 1:lambda
            Z(i) = f(x(i,:))
        end

        xp = zeros(mu,dim) //vecteur des x parents
        sp = zeros(mu,dim) // les ecarts types associés

        for j = 1:mu //selection 
            mini =  find(Z == min(Z))// indice de minf(x) 
            xp(j,:) = x(mini,:)//vecteur des mu minimum
            sp(j,:) = s(mini,:)
            Z(mini,:) = []//on l'enleve pour la prochaine iteration
            x(mini,:) = []
            s(mini,:) = []
        end

        xmp1 = sum(xp,'r')/mu //recombinaison
        smp1 = sum(sp,'r')/mu //autoadaptation
        n = n+1
    end

endfunction
   

/*****************   x0   ************************/
counteval = []
time = []
for i = 1:10:100
    x0 = i*ones(20,1)
    [xmin c t] = purecmaes(frozenbrock, x0, 0.5, 100000, 10^-8)
    counteval = [counteval c]
    time = [time t]
end
clf()
x0 = 1:10:100
plot(x0,time,'b.','thickness',2)
xlabel('$x_0$','fontsize',4)
ylabel('$temps$','fontsize',4)
yi = smooth([x0;time],0.01);
plot(yi(1,:)',yi(2,:)','r-');

/*****************   s0   ************************/
scf()
counteval = []
time = []
for s0 = 0.01:1:10
    [xmin c t] = purecmaes(frozenbrock, 10*ones(20,1), s0, 100000, 10^-8)
    counteval = [counteval c]
    time = [time t]
end
clf()
s0 = .01:1:10
plot(s0,time,'b.','thickness',2)
xlabel('$s_0$','fontsize',4)
ylabel('$temps$','fontsize',4)
yi = smooth([s0;time],0.01);
plot(yi(1,:)',yi(2,:)','r-');


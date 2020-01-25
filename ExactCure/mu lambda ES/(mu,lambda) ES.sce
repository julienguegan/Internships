clear

exec('C:\Users\Julien Guégan\Desktop\PFE\fonctions test.sce',-1)

f = rastrigin
x0 = [5 5]
pas0 = 3
lambda = 30
mu = 5
itermax = 500
tol = 0.0001


function  [xm,n] = evolutionstrategie(f,xm0,sm0,lambda,mu,itermax,tol) // μ < λ
    
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

[sol n] = evolutionstrategie(f,x0,pas0,lambda,mu,itermax,tol)
disp('x = ')
disp(sol)
disp('n = ')
disp(n)



exec('C:\Users\Julien Guégan\Desktop\PFE\fonctions test.sce',-1)

function  [xm,n] = evolutionstrategie(f,population,lambda,mu,itermax,tol,C) // μ < λ
    rand("normal")
    xm0 = sum(population)./length(population)//barycentre
    sm0 =  sqrt(1/length(population)*sum((population-xm0)^2))//xm0 - max(population) //ecart type moyen
    // ou utiliser la fonction nanstdev(x) de scilab
    tau = 1 //facteur d'apprentissage
    n = 1 //nbr iteration

    xmp1 = xm0
    xm = xm0
    smp1 = sm0

    while (((n<itermax)&(norm(xm-xmp1)>tol))| (n == 1)) then

        xm = xmp1
        sm = smp1
        /////*  AFFICHAGE */////
        clf()
        domaine=-20:0.1:20
        for i=-20:0.1:20
            ind_i=i*10+201
            z(ind_i) = f(i)
        end
        plot(domaine,z')
        ///////////////////////////////
        for i = 1:lambda //creation
            s(i) = sm * exp(tau*rand())
            x(i) = xm + s(i)*rand()
            if (-20<x(i)&x(i)<20)
                plot(x(i),f(x(i)),'k.','markersize',4)
            end
        end
        plot([2.9 2.9], [0 300],'g-')
        //disp('xm0 = '+string(xm)+' ; sm0 = '+string(sm))

        Z = 0
        for i = 1:length(x) //Une matrice Z = [x f(x) σ]
            Z(i,1) = x(i)
            Z(i,2) = f(x(i))
            Z(i,3) = s(i)
        end

        xp = 0 //vecteur des x parents
        sp = 0 // les ecarts types associés

        slct =  find(0 < C(Z(:,1)))// ceux qui respectent la contraintes
        W = 0
        for i = 1:length(slct)
            W(i,1) = Z(slct(i),1)
            W(i,2) = Z(slct(i),2)
            W(i,3) = Z(slct(i),3)
        end
        mu2 = mu
        if length(slct) < mu
            mu2 = length(slct)
        end
        
        for j = 1:mu2 //selection 
            mini =  find(W(:,2) == min(W(:,2)))// indice de minf(x) 
            xp(j) = W(mini,1)//vecteur des mu minimum
            sp(j) = W(mini,3)

            W(mini,:) = []//on l'enleve pour la prochaine iteration
        end

        xmp1 = sum(xp)/mu2 //recombinaison
        smp1 = sum(sp)/mu2 //autoadaptation
        plot(xmp1,f(xmp1),'r*')
        xs2png(gcf(),"image"+string(n))
        n = n+1

    end

endfunction

function C = contraintes(x)
    C = x-3
endfunction

//population = -10:10//;-25:0.1:-10]'

lambda = 15
mu = 5
fonction = rastr
itermax = 200
tol = 0.001
//population initiale
rand("uniform")
for i = 1:10
    population(i) = 28*rand()-14
end
plot(population,fonction(population),'k.','markersize',5)


[sol n] = evolutionstrategie(fonction,population,lambda,mu,itermax,tol,contraintes)
disp('le minimum est x = ')
disp(sol)
disp(' au bout de '+string(n)+' iterations')
/*
x=-15:0.1:15
for i=-15:0.1:15
    ind_i = i*10+151
    z(ind_i) = fonction(i)
end

plot(x',z,'b-')

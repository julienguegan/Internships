clear

exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\fonctions test.sce',-1)
exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\affichage.sce',-1)
exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\cercle.sce',-1)

function C = contraintes(x)
    C = (x(1)-2)^2+(x(2)-2)^2-1
endfunction

C = contraintes
f = rastrigin
rand("uniform")
for i =1:10
    population(i,1) = 3*rand()-1
    population(i,2) = 3*rand()-1
end
plot(population(:,1),population(:,2),'k.','markersize',2)

lambda = 40
mu = 10
itermax = 20
tol = 0.001

//function  [xm,n] = evolutionstrategie(f,population,lambda,mu,itermax,tol,C) // μ < λ
    longueur = size(population,1)
    dim = size(population,2)
    rand("normal")
    xm0 = sum(population,'r')./longueur//barycentre
    sm0 =  nanstdev(population,'r')
    tau = 1 //facteur d'apprentissage, petit = sigma plus petit
    n = 1 //nbr iteration

    xmp1 = xm0
    xm = xm0
    smp1 = sm0

    // ********* calcul F *********//
    xmax=4;xmin=-2.5;ymax=4;ymin=-2.5;
    xstep = (xmax-xmin)/115
    axe_x = xmin:xstep:xmax
    ystep = (ymax-ymin)/115
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
    ////////////////////////////////////


    while (n<itermax)//&(norm(xm-xmp1)>tol))| (n == 1)) then

        xm = xmp1
        sm = smp1

        //****  AFFICHAGE  ***//
        clf()
        xset("fpf"," ")
        contour2d(axe_x,axe_y,z,30)
        a = gcf()
        a.color_map = rainbowcolormap(64)
        //*******************************************************//

        for i = 1:lambda //creation
            s(i,:) = sm .* exp(tau*rand(1,dim))
            x(i,:) = xm + s(i).*rand(1,dim)
            if (xmin<x(i,1)&x(i,1)<xmax)&(ymin<x(i,2)&x(i,2)<ymax)
                plot(x(i,1),x(i,2),'k.','markersize',2)
            end
        end

        ellipse(2,2,1,1,0)
        Z = []
        cpt = 1
        slct=[]
        for i = 1:lambda
            Z(i,:) = f(x(i,:))
            if (C(x(i,:))<=0) //slct =  find(0 < C(x))
                slct(cpt) = i
                cpt = cpt +1
            end
        end
        if cpt>1 then //au moins 1 dans la contrainte 
            L = length(slct)
            W = []//zeros(L,1)
            Y = []//zeros(L,dim)
            Q = []//zeros(L,dim)
            for i = 1:L
                W(i) = Z(slct(i),:)
                Y(i,:) = x(slct(i),:)
                Q(i,:) = s(slct(i),:)
            end
            mu2 = mu
            if length(slct) < mu
                mu2 = length(slct)
            end
            xp = []//zeros(mu2,dim) //vecteur des x parents
            sp = []//zeros(mu2,dim)
            for j = 1:mu2 //selection 
                mini =  find(W == min(W))// indice de minf(x) 
                //disp(mini)
                xp(j,:) = Y(mini,:)//vecteur des mu minimum
                sp(j,:) = Q(mini,:)
                W(mini) = []//on l'enleve pour la prochaine iteration
                Y(mini,:) = []
                Q(mini,:) = []
            end
            xmp1 = sum(xp,'r')/mu2
            smp1 = sum(sp,'r')/mu2 
           // plot(xmp1(1),xmp1(2),'b*')
        end
        n = n+1
        xs2png(gcf(),"esconstr"+string(n))
    end

//endfunction

disp('x = ')
disp(xmp1)
disp('n = ')
disp(n)



clear

exec('C:\Users\Julien Guégan\Desktop\PFE\algorithmes\fonctions test.sce',-1)
exec('C:\Users\Julien Guégan\Desktop\PFE\algorithmes\affichage.sce',-1)

function  [xm,n] = evolutionstrategie(f,xm0,sm0,lambda,mu,itermax,tol) // μ < λ
    //longueur = size(population,1)
    //dim = size(population,2)
    longueur = length(xm0)
    dim = length(sm0)
    rand("normal")
    //xm0 = sum(population,'r')./longueur//barycentre
    //sm0 =  nanstdev(population,'r')
    tau = 1 //facteur d'apprentissage, petit = sigma plus petit
    n = 1 //nbr iteration

    xmp1 = xm0
    xm = xm0
    smp1 = sm0

    // ********* calcul F *********//
    xmax=4;xmin=-4;ymax=4;ymin=-4;
    xstep = (xmax-xmin)/75
    axe_x = xmin:xstep:xmax
    ystep = (ymax-ymin)/75
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

    while (n<itermax)&((norm(xm-xmp1)>tol)|n==1)

        xm = xmp1
        sm = smp1

        //****  AFFICHAGE  ***//
        clf()
        xset("fpf"," ")
        contour2d(axe_x,axe_y,z,30)
        a = gcf()
        a.color_map = rainbowcolormap(64)
        xlabel('$x_1$','fontsize',4)
        ylabel('$x_2$','fontsize',4)
        //*******************************************************//

        //affiche(f,4,-4,4,-4,'contour')
        for i = 1:lambda //creation
            s(i,:) = sm .* exp(tau*rand(1,dim))
            x(i,:) = xm + s(i).*rand(1,dim)
            if (xmin<x(i,1)&x(i,1)<xmax)&(ymin<x(i,2)&x(i,2)<ymax)
                plot(x(i,1),x(i,2),'k.','markersize',3)
            end
        end

        //plot(x(:,1),x(:,2),'k.','markersize',2)
        Z = 0
        for i = 1:lambda
            Z(i,:) = f(x(i,:))
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
        plot(xp(:,1),xp(:,2),'b.','markersize',2)//les mu parents retenu

        xmp1 = sum(xp,'r')/mu //recombinaison
        smp1 = sum(sp,'r')/mu //autoadaptation
        plot(xmp1(1),xmp1(2),'b*')
       xs2png(gcf(),"tmp"+string(n))
        n = n+1
        
    end

endfunction

fonction = rastrigin
rand("uniform")
xm0 = [0,0]
sm0 = [1.5,1.5]
lambda = 30
mu = 10
itermax = 100
tol = 0.0001
[sol n] = evolutionstrategie(fonction,xm0,sm0,lambda,mu,itermax,tol)
disp('x = ')
disp(sol)
disp('n = ')
disp(n)



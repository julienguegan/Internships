clear
clf()

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tinfini\Fonctions 2 cohorts - Tinfini.sci',-1)

T = 120
t = 0:T

// CI
M01 = 400 
M02 = 400

n1 = 100
n2 = 1
n = n1+n2

Mee = fsolve(0,eqnM)
disp(" M** = "+string(Mee))

Mσ = 1500//Mσ au debut 

tic()

n1 = 100
n2 = 1

M1σ = Mσ
M2σ = Mσ+500//*(rand()-0.5)


disp("-------------------")
disp("M1σ = "+string(M1σ))
disp("M2σ = "+string(M2σ))

tic()

i = 0 //compteur
norme = 5 //initialisation critere d'arret
while ((i < 100) & (n1 >0.5) & (n2>0.5) & (norme>10^-4))

    tolM = 10^-2
    tolJ = 10^-2

    //resolution ODE 
    M = ode("stiff",[M01;M02;200;200],0,t,[tolM tolM tolJ tolJ],syst_dyn)  
    J1 = M(3,length(t))
    J2 = M(4,length(t)) //derniere valeur du vecteur (au temps T)

    n1old = n1//stockage des valeurs précédentes
    n2old = n2

    //n nouvelle taille de la cohort
    n11 = ((n1*J1)/(n1*J1+n2*J2))*n 
    n21 = ((n2*J2)/(n1*J1+n2*J2))*n
    n1 = n11
    n2 = n21
    norme = sqrt((n1-n1old)^2+(n2-n2old)^2)

    subplot(211)
    plot2d(t,M(1,:),color('red'))
    plot2d(t,M(2,:),color('blue'))
    plot([0,T],[Mee,Mee],'k--')

    xstring(length(t)*0.7, Mee*0.9, "M** = "+string(Mee));
    gce().font_size = 3;
    xstring(length(t)*0.7, Mee*0.7, "M1* = "+string(M1σ));
    gce().font_size = 3;
    gce().font_color = 5; 
    xstring(length(t)*0.7, Mee*0.5, "M2* = "+string(M2σ));
    gce().font_size = 3; 
    gce().font_color = 2; 
    xlabel("t - jours")
    ylabel("M")
    xtitle("Evolution de M1 et M2 au long d''une saison")

    subplot(212)
    a = gca()
    a.data_bounds=[0,0;120,12]
    plot(i,n1,'r.-')
    plot(i,n2,'b.-')
    legend(["n1";"n2"],-1, %f)
    xlabel("i - saisons")
    ylabel("n")
    xtitle("Evolution de n1 et n2 au long de i saisons")

    i = i+1
end

if (10-n1 > 0) then //n1 diminue = strategie perdante 
    Mσ = M2σ //la 2e est prise comme strategie de reference pour la prochaine iteration
else //sinon elle est gagnante 
    Mσ = M1σ //c'est la strategie de ref. pour la proch. 
end
xstring(100, n1-1, "$n_1$")//" = "+string(n1))
gce().font_size = 3; 
gce().font_color = 5;
xstring(100, n2, "$n_2$")//" = "+string(n2))
gce().font_size = 3; 
gce().font_color = 2; 
disp("n1 = "+string(n1))
disp("n2 = "+string(n2))
disp("iterations : "+string(i))
disp("||n|| = "+string(norme))
disp("t = "+string(toc())+" secondes")

disp("TIME = "+string(toc())+" s")

xstring(50, 5, "t = "+string(round(toc()))+" s")

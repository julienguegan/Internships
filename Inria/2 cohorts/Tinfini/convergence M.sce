clear

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tinfini\Fonctions 2 cohorts - Tinfini.sci',-1)

T = 140
t = 0:T

// CI
M01 = 400 
M02 = 400

n1 = 100
n2 = 1
n = n1+n2

Mee = fsolve(0,eqnM)
disp(" M** = "+string(Mee))

Mσ = 1750

gagnants = [] //on stock les Mσ gagnants
perdants = [] //les perdants
egalite = [] //on stock quand il y a egalité
tic()
j = 0
while (j<250) 

    n1 = 100
    n2 = 1

    M1σ = Mσ
    M2σ = Mσ+15*(rand()-0.5)

        atol = 10^-4
        rtol = 10^-4
        M = ode("stiff",[M01;M02;800;800],0,t,[atol atol rtol rtol],syst_dyn)
        J1 = M(3,length(t))
        J2 = M(4,length(t))

        n11 = (n1*J1*n)/(n1*J1+n2*J2)
        n21 = (n2*J2*n)/(n1*J1+n2*J2)

        //mise a jour
        n1 = n11
        n2 = n21

    if ((100-n1) > 0) then //n1 diminue = strategie perdante 
        Mσ = M2σ //la 2e est prise comme strategie de 'reference' pour la prochaine iteration
        perdants =  [perdants;M1σ ]
        gagnants = [gagnants;M2σ ]
    else //sinon M1 augmente donc on le garde
        perdants =  [perdants;M2σ ]
        gagnants = [gagnants;M1σ ]
        Mσ = M1σ
    end
    j = j+1
end

disp("TIME = "+string(round(toc()))+" sec")

scf()
i = 1:size(gagnants,1)
Mee = fsolve(0,eqnM)

plot([0,length(i)],[Mee,Mee],'r--')
xstring(length(i)/2, Mee, "M** = "+string(Mee)); gce().font_color = 5;
gce().font_size = 3;
//plot(i,perdants,'k.')
plot(i,gagnants,'b')
//j = 1:size(egalite,2)
//plot(j,egalite,'g.')
//xstring(30, 3100, "tol = "+string(10^-4));
//xstring(30, 3000, "t = "+string(966)+" s")
xlabel("i - temps évolutif","fontsize", 3)
ylabel("$M_σ$","fontsize", 3)
gca().data_bounds = [0,1700;250,2000]
title(' T = '+string(T),'fontsize',4)

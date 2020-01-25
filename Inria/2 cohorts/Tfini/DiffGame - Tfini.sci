clear
clf()

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tfini\Fonctions 2 cohorts - T fini.sci',-1)

T = 20
t = 0:T

// CI
M01 = 800 
M02 = 800

n1 = 100
n2 = 1
n = n1+n2

Mee = fsolve(0,eqnM)
τee = fsolve(0,list(lambda,Mee))
disp(" M** = "+string(Mee)+" , τ** = "+string(τee))

Mσ = 1700//Mσ au debut 

donnees = x_mdialog("Mσ",["M1σ";"M2σ"],['1250';'1750'])
M1σ = evstr(donnees(1))
M2σ = evstr(donnees(2))

τ1σ =10//fsolve(0,list(lambda,M1σ))
τ2σ = 10//fsolve(0,list(lambda,M2σ))

disp("-------------------")
disp("M1σ = "+string(M1σ)+", τ1σ = "+string(τ1σ))
disp("M2σ = "+string(M2σ)+", τ2σ = "+string(τ2σ))

tic()

affichage_n = [] //stockage de nos résultats à afficher
affichage_M = []

norme = 5 //initialisation critere d'arret
i = 0
while ((i < 150)  )

    affichage_n = [affichage_n;n1 n2]

    tolM = 10^-3
    tolJ = 10^-3

    //resolution ODE 
    M = ode("stiff",[M01;M02;200;200],0,t,[tolM tolM tolJ tolJ], syst_dyn)  
    J1 = M(3,length(t))
    J2 = M(4,length(t)) //derniere valeur du vecteur (au temps T)
    affichage_M = [affichage_M M]
    n1old = n1//stockage des valeurs précédentes
    n2old = n2

    //n nouvelle taille de la cohort
    n11 = ((n1*J1)/(n1*J1+n2*J2))*n 
    n21 = ((n2*J2)/(n1*J1+n2*J2))*n
    n1 = n11
    n2 = n21
    norme = sqrt((n1-n1old)^2+(n2-n2old)^2)

    i = i+1
end
/*
long = size(affichage_M,2)
subplot(211)
plot(t,affichage_M(1:2,long-T:long))
plot([0,T],[Mee,Mee],'k--')
xstring(length(t)*0.8, Mee*0.8, "M** = "+string(round(Mee)));
xlabel("t")
ylabel("M(t)")
xtitle("Evolution des M1 et M2 au long d''une saison")
*/
k = 1:150
//subplot(212)
xlabel("i - saisons",'fontsize',3)
ylabel("n(i) - densité de la lésion",'fontsize',3)
plot(k,affichage_n(:,1),'r')
plot(k,affichage_n(:,2),'b')

xstring(140, n1+6, "$n_1$")
gce().font_size = 5; 
gce().font_color = 5;
xstring(140, n2, "$n_2$")
gce().font_size = 5; 
gce().font_color = 2; 
//legend(['n1';'n2'],-1, %f)
//title("Evolution des n1 et n2 au long de k saisons",'fontsize',3)






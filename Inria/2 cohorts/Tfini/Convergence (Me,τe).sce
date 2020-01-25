clear

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tfini\Fonctions 2 cohorts - T fini.sci',-1)

T = 30
t = 0:T

// CI
M01 = 400 
M02 = 400

n1 = 100
n2 = 1
n = n1+n2

Mee = fsolve(0,eqnM)
disp(" M** = "+string(Mee))
τee = fsolve(0,list(lambda,Mee))
disp('τ* = '+string(τee))

Mσ = 2100
τσ = 12

gagnants = [] //on stock les Mσ gagnants
perdants = [] //les perdants
egalite = [] //on stock quand il y a egalité
tic()
j = 0
while (j<300) 

    n1 = 100
    n2 = 1

    M1σ = Mσ
    M2σ = Mσ+10*(rand()-0.5)
    τ1σ = τσ
    τ2σ = τσ+0.2*(rand()-0.5)

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
        τσ = τ2σ
        perdants =  [perdants;M1σ τ1σ]
        gagnants = [gagnants;M2σ τ2σ]
    else //sinon M1 augmente donc on le garde
        perdants =  [perdants;M2σ τ2σ]
        gagnants = [gagnants;M1σ τ1σ]

        Mσ = M1σ
        τσ = τ1σ
    end
    j = j+1
end

disp("TIME = "+string(round(toc()))+" sec")


Mσgagnants = gagnants(:,1)'
Mσperdants = perdants(:,1)'
τσgagnants = gagnants(:,2)'
τσperdants = perdants(:,2)'

i = 1:size(gagnants,1)
Mee = fsolve(0,eqnM)
τee = fsolve(0,list(lambda,Mee))


subplot(2,2,2)
plot([0,length(i)],[Mee,Mee],'r--')
xstring(length(i)/2, Mee, "M** = "+string(Mee)); gce().font_color = 5;
gce().font_size = 3;
//plot(i,Mσperdants,'k.')
plot(i,Mσgagnants,'b')
//j = 1:size(egalite,2)
//plot(j,egalite,'g.')
//xstring(30, 3100, "tol = "+string(10^-4));
//xstring(30, 3000, "t = "+string(966)+" s")
xlabel("i - temps évolutif")
ylabel("Mσ","fontsize", 4)
title("Evolution des Mσ","fontsize", 3)

subplot(2,2,1)
plot([0,length(i)],[τee,τee],'r--')
xstring(length(i)/2, τee, "τ** = "+string(τee)); gce().font_color = 5;
gce().font_size = 3;
//plot(i,τσperdants,'k.')
plot(i,τσgagnants,'b')
xlabel("i - temps évolutif")
ylabel("τσ","fontsize", 4)
title("Evolution des τσ","fontsize", 3)


subplot(2,2,3)
//plot(τσperdants,Mσperdants,'k.')
plot(Mσgagnants,τσgagnants,'r.')
plot([Mee,Mee],[8,τee],'k--')
xstring(Mee, 8, "M** = "+string(Mee)); 
gce().font_color = 5;
plot([1700,Mee],[τee,τee],'k--')
xstring(1700,τee,"τ** = "+string(τee)); 
gce().font_color = 5;
xlabel("$Mσ$","fontsize", 4)
ylabel("$τσ$","fontsize", 4)
title("Evolution du couple (Mσ,τσ)","fontsize", 3)


subplot(2,2,4)
param3d(i,τσgagnants,Mσgagnants)
e=gce() //the handle on the 3D polyline
e.foreground=color('cyan');
e.polyline_style=4
L = ones(1,200)
param3d(1:200,τee*L,Mee*L)
gce().foreground=color('red');
gce().thickness = 4
zlabel("Mσ","fontsize", 4)
ylabel("τσ","fontsize", 4)
xlabel("i")
title("Evolution du couple (Mσ,τσ)","fontsize", 3)


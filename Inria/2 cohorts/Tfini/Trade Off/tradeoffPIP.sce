clear

//T < +∞

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tfini\Fonctions 2 cohorts - T fini.sci',-1)

T = 40
t = 0:T

// CI
M01 = 400 
M02 = 400

n1 = 100
n2 = 1
n = n1+n2
Mee = fsolve(0,eqnM)
τee = fsolve(0,list(lambda,Mee))
disp(" M** = "+string(Mee)+" , τ** = "+string(τee))

gagnants = [] //on stock les Mσ gagnants
perdants = [] //les perdants
egalité = [] //on stock quand il y a egalité
tic()

stepM = 10 
Mσ =  1700:stepM:2100

load("C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tfini\Trade-Off\tradeoff_save",'tradeoff_τ','tradeoff_M','Metoileμγ','τetoileμγ')
scf()
x = tradeoff_M
plot(x, tradeoff_τ,'b-','thickness',3)
plot(Metoileμγ,τetoileμγ,'k.')
ylabel('$τ*$','fontsize',4)
xlabel('$M*$','fontsize',4)
title('$Compromis\ Evolutif$','fontsize',3)

hm = stepM
for Mx = 1:length(Mσ)
    disp('M1σ = '+string(Mσ(Mx)))
    for My = 1:length(Mσ)
        
        n1 = 100
        n2 = 1
        
        M1σ = Mσ(Mx)
        index = find(x == M1σ)  
        τ1σ = tradeoff_τ(index)
        M2σ = Mσ(My)
        index = find(x == M2σ)
        τ2σ = tradeoff_τ(index) 
        M = ode("stiff",[M01;M02;800;800],0,t, syst_dyn) 
        J1 = M(3,length(t))
        J2 = M(4,length(t))
        n11 = (n1*J1*n)/(n1*J1+n2*J2) 
        n21 = (n2*J2*n)/(n1*J1+n2*J2)
        s(Mx,My) = n21/n2

        if ((100-n11) > 0) then 
            gagnants = [gagnants M2σ]
            perdants = [perdants M1σ]
        else //sinon M1 augmente donc on le garde
            gagnants = [gagnants M1σ]
            perdants = [perdants M2σ]
        end
        
    end
end
disp("TIME = "+string(toc()/60)+" min")

scf()
plot(perdants,gagnants,'r.')
plot(gagnants,perdants,'b.')
contour(Mσ,Mσ,s,20)
plot(Mee,Mee,'k.')

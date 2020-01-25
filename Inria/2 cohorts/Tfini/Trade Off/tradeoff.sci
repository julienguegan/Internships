clear

//T < +∞
load("C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tfini\Trade-Off\tradeoff_save",'tradeoff_τ','tradeoff_M','Metoileμγ','τetoileμγ')
scf()
x = tradeoff_M
plot(x, tradeoff_τ,'b-','thickness',3)
plot(Metoileμγ,τetoileμγ,'k.')
ylabel('$τ*$','fontsize',4)
xlabel('$M*$','fontsize',4)
title('$Compromis\ Evolutif$','fontsize',3)
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

stepM = 20 
Mσ =  1600:stepM:2200
hm = stepM
for Mx = 2:length(Mσ)-1
    disp('M1σ = '+string(Mσ(Mx)))
    for My = 1:3//tribande
        
        n1 = 100
        n2 = 1
        
        tribande = [-1 0 1]
        M1σ = Mσ(Mx)
        index = find(x == M1σ)  
        τ1σ = tradeoff_τ(index)
        M2σ = Mσ(Mx)+tribande(My)
        index = find(x == M2σ)
        τ2σ = tradeoff_τ(index) 
        M = ode("stiff",[M01;M02;800;800],0,t, syst_dyn) 
        J1 = M(3,length(t))
        J2 = M(4,length(t))
        n11 = (n1*J1*n)/(n1*J1+n2*J2) 
        n21 = (n2*J2*n)/(n1*J1+n2*J2)
        s(Mx-1,My) = n21/n2

        if ((95-n1) > 0) then 
            gagnants = [gagnants M2σ]
            perdants = [perdants M1σ]
        else //sinon M1 augmente donc on le garde
            gagnants = [gagnants M1σ]
            perdants = [perdants M2σ]
        end
        
    end
end
disp("TIME = "+string(toc()/60)+" min")

for Mx = 3:length(Mσ)-2
        dsdMy(Mx-2) = (s(Mx,3)-s(Mx,1))/(2*hm) 
end

scf() 
ind = 3:length(Mσ)-2 // ! le 1er et le dernier elmts n'existe plus avec la derivée
plot(Mσ(ind),dsdMy)
plot([Mσ(1),Mσ(length(ind))],[0,0],'r--')
//xarrows([1400:50:1800;1410:50:1850],zeros(2,9),400)
//xarrows([2310:-50:1950;2300:-50:1920],zeros(2,8),400)
plot([Mee,Mee],[0,0],'k.')
plot([Mee,Mee],[dsdMy(1),dsdMy(length(ind))],'r--')
xstring(Mee,dsdMy(length(ind)),"M* = "+string(Mee))
gce().font_color = 5;
Mmax = interp1(dsdMy,Mσ(ind),0) //quand ds/dy(x,x)=0, interpolation de l'indice
plot([Mmax,Mmax],[dsdMy(1),dsdMy(length(ind))],'k--')
xstring(Mmax-150,dsdMy(length(ind)),"Mmax = "+string(Mmax))
//title("$gradient\ de\ selection : \frac{\partial s}{\partial y}(x,y=x) $",'fontsize',3)
xlabel("Mσ")

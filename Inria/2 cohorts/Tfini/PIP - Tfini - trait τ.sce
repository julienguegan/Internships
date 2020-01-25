clear

//T < +∞

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tfini\Fonctions 2 cohorts - T fini.sci',-1)

T = 20
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

step = 0.4
τσ = 7:step:15
gagnants = [] //on stock les Mσ gagnants
perdants = [] //les perdants
egalité = [] //on stock quand il y a egalité
tic()

M1σ = 1800
M2σ = 2000
for r = 1:length(τσ)
    disp(r)
    for m = 1:length(τσ)

        n1 = 100
        n2 = 1

        τ1σ = τσ(r)
        τ2σ = τσ(m)
        /*
        disp("-----------------")
        disp("M1σ = "+string(M1σ))
        disp("M2σ = "+string(M2σ))
        */
        i = 0 //compteur
        norme = 5
        while ((i < 2) & (n1>0.2) & (n2>0.2) & (norme>10^-5))

            tolM = 10^-6
            tolJ = 10^-6
            M = ode("stiff",[M01;M02;800;800],0,t,[tolM tolM tolJ tolJ], syst_dyn) 
            J1 = M(3,length(t))
            J2 = M(4,length(t))

            n1old = n1 //stockage des valeurs précédentes pour le critere d'arret
            n2old = n2

            n11 = (n1.*J1.*n)./(n1.*J1+n2.*J2) //vecteur
            n21 = (n2.*J2.*n)./(n1.*J1+n2.*J2)
            n1 = n11 
            n2 = n21

            i = i+1
            norme = sqrt((n1-n1old)^2+(n2-n2old)^2)

        end

        if ((100-n1) > 0) then 
            gagnants = [gagnants τ2σ]
            perdants = [perdants τ1σ]
        else //sinon M1 augmente donc on le garde
            gagnants = [gagnants τ1σ]
            perdants = [perdants τ2σ]
        end

        s(r,m) = n21/n2old
        
    end
end

disp("TIME = "+string(toc()/60)+" min")
scf()
plot(perdants,gagnants,'r.')
plot(gagnants,perdants,'b.')
contour(τσ,τσ,s,20)
xlabel("τ1σ")
ylabel("τ2σ")

f=scf()
f.color_map = rainbowcolormap(32);
surf(s)
scf() //affichage de s(x,x*) avec ds/dy=0
X = 1:length(τσ)
snglr = find(τσ == 11)
plot(τσ,s(snglr,X))
xlabel("$y$")
ylabel("$s$")
title("$s(x*,x)$",'fontsize',3)
ind = find(s(snglr,X)==max(s(snglr,X))) 
maxi = τσ(ind)
plot([τee,τee],[s(snglr,1),s(snglr,ind)],'r--')
xstring(maxi-1,s(snglr,1),"Mmax = "+string(maxi))
plot([maxi,maxi],[s(snglr,1),s(snglr,ind)],'k--')
xstring(τee,s(snglr,1),"M* = "+string(τee))
gce().font_color = 5;
xstring(10,s(snglr,snglr), "$\frac{\partial s}{\partial y} = 0$")
plot([10.5 12],[s(snglr,ind) s(snglr,ind)],'r--')


h = step
for x = 2:length(τσ)-1
    gradselect(x-1) = (s(x,x+1)-s(x,x-1))/(2*h) //difference centrée
end

scf() //affichage de ds/dy(x,x)
y = 2:length(τσ)-1
plot(τσ(y),gradselect)
plot([τσ(1),τσ(length(y))],[0,0],'r--')
xarrows([9.5:0.2:10.5;9.6:0.2:10.6],zeros(2,6),1.5)
xarrows([12.1:-0.2:11;12:-0.2:11],zeros(2,6),1.5)
plot([τee,τee],[0,0],'k.')
plot([τee,τee],[gradselect(1),gradselect(length(y))],'r--')
xstring(τee,gradselect(1),"τ* = "+string(round(τee)))
gce().font_color = 5;
Mmax = interp1(gradselect,τσ(y),0)
plot([Mmax,Mmax],[gradselect(1),gradselect(length(y))],'k--')
xstring(Mmax-150,gradselect(length(y)),"M = "+string(round(Mmax)))
title("$gradient\ de\ selection : \frac{\partial s}{\partial y}(x,y=x) $",'fontsize',3)
xlabel("τσ")


hesselecty = (s(snglr,snglr+1)- 2*s(snglr,snglr)+s(snglr,snglr-1))./(h^2)
hesselectx = (s(snglr+1,snglr)- 2*s(snglr,snglr)+s(snglr-1,snglr))./(h^2)

if(hesselecty<0) then
    disp('ESS')
end
if (hesselecty<hesselectx) then
    disp('stable par convergence')
end
if (hesselecty>0)&(hesselectx>0) then
    disp('mutuellement invasible')
end

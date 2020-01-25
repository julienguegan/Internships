clear

//T < +∞

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tinfini\Fonctions 2 cohorts - Tinfini.sci',-1)

T = 150
t = 0:T

// CI
M01 = 400 
M02 = 400

n1 = 100
n2 = 1
n = n1+n2

Mee = fsolve(0,eqnM)
disp(" M** = "+string(Mee))

step = 35
Mσ = 1600:step:2200
gagnants =[] //on stock les Mσ gagnants
perdants = [] //les perdants
tic()
for r = 1:length(Mσ)
    disp(r)
    for m = 1:length(Mσ)

        n1 = 100
        n2 = 1

        M1σ = Mσ(r)
        M2σ = Mσ(m)

        i = 0 //compteur
        norme = 5
       // while ((i < 10) & (n1>1) & (n2>1) & (norme>10^-5))

            tolM = 10^-5
            tolJ = 10^-5
            M = ode("stiff",[M01;M02;800;800],0,t,[tolM tolM tolJ tolJ], syst_dyn) 
            J1 = M(3,length(t))
            J2 = M(4,length(t))
            n1old = n1
            n2old = n2
            n11 = (n1*J1*n)/(n1*J1+n2*J2) 
            n21 = (n2*J2*n)/(n1*J1+n2*J2)
            n1 = n11 
            n2 = n21
            norme = sqrt((n1-n1old)^2+(n2-n2old)^2)

            if i==0 then
                s(m,r) = n21/n2old
            end
            i = i+1
       // end
      //  disp("t = "+string(toc())+" s") 
        if ((100-n1) > 0) then 
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
/*subplot(311)
plot(perdants,gagnants,'r.')
subplot(312)
plot(gagnants,perdants,'b.')
subplot(313)
plot(perdants,gagnants,'r.')
plot(gagnants,perdants,'b.')
contour(Mσ,Mσ,s,20)
xlabel("M1σ")
ylabel("M2σ")*/
M = 1600:2200
plot(M,M)
gca().data_bounds=[1500 1500;2300 2300]
plot(perdants,gagnants,'r.')
plot(gagnants,perdants,'b.')
//plot(perdants,perdants)
xstring(1625,1800,'+')
gce().font_size = 5; 
gce().font_color = 5; 
xstring(2175,2000,'+')
gce().font_size = 5; 
gce().font_color = 5; 
xstring(1750,1525,'_')
gce().font_size = 5; 
gce().font_color = 2; 
xstring(1750,2275,'_')
gce().font_size = 5; 
gce().font_color = 2; 
plot(Mee,Mee,'k.')
plot([Mee Mee],[1600 Mee],'k--')
xstring(1900,1600,'M* = '+string(Mee))
xlabel('$M_r$','fontsize',4)
ylabel('$M_m$','fontsize',4)
//contour(Mσ,Mσ,s,7)

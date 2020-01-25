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

Mσ = 1600:50:2200
gagnants = [] //on stock les Mσ gagnants
perdants = [] //les perdants
egalité = [] //on stock quand il y a egalité
tic()
τ1σ = τee
τ2σ = τee
for r = 1:length(Mσ)
    for m = 1:length(Mσ)

        n1 = 100
        n2 = 1

        M1σ = Mσ(r)
        M2σ = Mσ(m)
        /*
        disp("-----------------")
        disp("M1σ = "+string(M1σ))
        disp("M2σ = "+string(M2σ))
        */
        i = 0 //compteur
        norme = 5
        while ((i < 2) & (n1>0.2) & (n2>0.2) & (norme>10^-5))

            tolM = 10^-5
            tolJ = 10^-5
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
            gagnants = [gagnants M2σ]
            perdants = [perdants M1σ]
        else //sinon M1 augmente donc on le garde
            gagnants = [gagnants M1σ]
            perdants = [perdants M2σ]
        end

        s(r,m) = n21/n2old
        
    end
end

disp("TIME = "+string(toc()/60)+" min")
scf()
plot(perdants,gagnants,'r.')
plot(gagnants,perdants,'b.')
contour(Mσ,Mσ,s,20)
xlabel("M1σ")
ylabel("M2σ")


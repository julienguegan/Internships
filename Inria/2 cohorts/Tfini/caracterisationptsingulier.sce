clear

//T < +∞

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tfini\Fonctions 2 cohorts - T fini.sci',-1)

T = 40
t = 0:T

// CI
M01 = 400 
M02 = 400
tol = 10^-8

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


stepτ = 0.1
τσ = 8:stepτ:14
stepM = 20
Mσ =  1700:stepM:2100

hτ = stepτ
hm = stepM

Mr = [Mee-hm Mee Mee+hm]
τr = [τee-hτ τee τee+hτ]
for i = 1:3 //M1σ
    for j = 1:3 //τ1σ
        for v = 1:3 //M2σ
            for w = 1:3 //τ2σ
                n1 = 100
                n2 = 1
                τ1σ = τr(j)
                M1σ = Mr(i)
                M2σ = Mr(v)
                τ2σ = τr(w)
                M = ode("stiff",[M01;M02;800;800],0,t,tol, syst_dyn) 
                J1 = M(3,length(t))
                J2 = M(4,length(t))
                n11 = (n1*J1*n)/(n1*J1+n2*J2) 
                n21 = (n2*J2*n)/(n1*J1+n2*J2)
                s(i,j,v,w) = n21/n2 //$s(M_r,τ_r,M_m,τ_m)$
            end
        end
    end
end
H11 = (s(2,2,3,2)-2*s(2,2,2,2)+s(2,2,1,2))/(hm^2) //Mm²
H12 = (s(2,2,3,3)-s(2,2,1,3)-s(2,2,3,1)+s(2,2,1,1))/(4*hm*hτ) //Mmτm
H22 = (s(2,2,2,3)-2*s(2,2,2,2)+s(2,2,2,1))/(hτ^2) //τm²
H = 10^6*[H11,H12;H12,H22]

J11 = (s(3,2,3,2)-s(1,2,3,2)-s(3,2,1,2)+s(1,2,1,2))/(4*hm*hτ) //MmMr
J12 = (s(2,3,3,2)-s(2,1,3,2)-s(2,3,1,2)+s(2,1,1,2))/(4*hm*hτ) //Mmτr
J21 = (s(3,2,2,3)-s(1,2,2,3)-s(3,2,2,1)+s(1,2,2,1))/(4*hm*hτ) //τmMr
J22 = (s(2,3,2,3)-s(2,1,2,3)-s(2,3,2,1)+s(2,1,2,1))/(4*hm*hτ) // τmτr
J = 10^6*[H11+J11,H12+J12;H12+J21,H22+J22]
disp('H = ')
disp(H)
disp('J = ')
disp(J)

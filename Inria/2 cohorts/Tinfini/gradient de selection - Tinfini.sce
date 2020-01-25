clear

//T < +∞

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tinfini\Fonctions 2 cohorts - Tinfini.sci',-1)

T = 300
t = 0:T

// CI
M01 = 400 
M02 = 400

n1 = 100
n2 = 1
n = n1+n2

Mee = fsolve(0,eqnM)
disp(" M** = "+string(Mee))

step = 20
Mσ = 1400:step:2300
gagnants =[] //on stock les Mσ gagnants
perdants = [] //les perdants

tic()
for r = 1:length(Mσ)
    for m = 1:3 //on regarde une tribande autour de y=x
        n1 = 100
        n2 = 1

        tribande = [-1 0 1]
        M1σ = Mσ(r)
        M2σ = Mσ(r)+tribande(m)

        tol = 10^-7
        M = ode("stiff",[M01;M02;800;800],0,t,tol, syst_dyn) 
        J1 = M(3,length(t))
        J2 = M(4,length(t))
        n1old = n1
        n2old = n2
        n11 = (n1*J1*n)/(n1*J1+n2*J2) 
        n21 = (n2*J2*n)/(n1*J1+n2*J2)
        n1 = n11 
        n2 = n21
        //m+2 en indice car m = -1 0 1
        s(r,m) = n21/n2old
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

h = step
for x = 2:length(Mσ)-1
    gradselect(x-1) = (s(x,3)-s(x,1))/(2*h) //difference centrée
end
//gradselect = 10^6*gradselect
scf() //affichage de ds/dy(x,x)
ind = 2:length(Mσ)-1 // ! le 1er et le dernier elmts n'existe plus avec la derivée
plot(Mσ(ind),gradselect)
plot([Mσ(1),Mσ(length(ind))],[0,0],'r--')
//xarrows([1400:50:1800;1410:50:1850],zeros(2,9),400)
//xarrows([2310:-50:1950;2300:-50:1920],zeros(2,8),400)
plot([Mee,Mee],[0,0],'k.')
plot([Mee,Mee],[gradselect(1),gradselect(length(ind))],'r--')
xstring(Mee,gradselect(length(ind)),"$M^{**} = $"+string(Mee)+"$")
gce().font_color = 5;
gce().font_size = 4;
Mmax = interp1(gradselect,Mσ(ind),0) //quand ds/dy(x,x)=0, interpolation de l'indice
plot([Mmax,Mmax],[gradselect(1),gradselect(length(ind))],'k--')
xstring(Mmax-200,gradselect(length(ind)),"$M_{opt} = $"+string(Mmax)+"$")
gce().font_size = 4;
//title("$gradient\ de\ selection : \frac{\partial s}{\partial y}(x,y=x) $",'fontsize',3)
xlabel("$M_r$",'fontsize',4)
ylabel("$\left.\frac{\partial s}{\partial M_m}\right| _{M_m=M_r}$",'fontsize',3,'rotation',0)

snglr = find(Mσ == 1900)
hesselecty = (s(snglr,3)- 2*s(snglr,2)+s(snglr,1))./(h^2)
hesselectx = (s(snglr+1,2)- 2*s(snglr,2)+s(snglr-1,2))./(h^2)

if(hesselecty<0) then
    disp('ESS')
end
if (hesselecty<hesselectx) then
    disp('stable par convergence')
end
if (hesselecty>0)&(hesselectx>0) then
    disp('mutuellement invasible')
end


//*******EQUATION CANONIQUE***********//
/*
Mxd = Mσ

Mx = 1700:2100

function dy  = canonique(t,y) 
   
    dy =10^7*interp1(Mxd(ind), gradselect,y)
    
endfunction

scf()

temps = 0:100/*marche pas ?
y = 1700:25:2100
fchamp(canonique,0,y,temps)*/

for i=1500:100:2200
trajectoire = ode(i,0,temps,canonique) 
Mpointx = trajectoire(1,:)
plot(temps,Mpointx,"r-")
end

title('$canonical\ equation\ : \dot{M_x} \ =  \frac{\partial s}{\partial M_y}(M_x,τ_x) $','fontsize',3)
xlabel('$temps $','fontsize',4)
ylabel('$M_x $','fontsize',4)
xstring(1970,13,'$\dot{M_x} = 0 $')
gce().font_size = 3; 
gce().font_color = 2; 

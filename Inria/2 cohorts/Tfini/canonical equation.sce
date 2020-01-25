clear

//T < +∞

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tfini\Fonctions 2 cohorts - T fini.sci',-1)

T = 40
t = 0:T

// CI
M01 = 400 
M02 = 400
tol = 10^-5

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


stepτ = 0.5
τσ = 8:stepτ:14
stepM = 20
Mσ =  1700:stepM:2100

hτ = stepτ
hm = stepM
for Mx = 1:length(Mσ)
    disp(Mx)
    for τx = 1:length(τσ)

        n1 = 100
        n2 = 1

        τ1σ = τσ(τx)
        M1σ = Mσ(Mx)

        //cas 1 : My= Mx+h et τy = τx
        M2σ = Mσ(Mx)+hm
        τ2σ = τσ(τx)
        M = ode("stiff",[M01;M02;800;800],0,t,tol, syst_dyn) 
        J1 = M(3,length(t))
        J2 = M(4,length(t))
        n11 = (n1*J1*n)/(n1*J1+n2*J2) 
        n21 = (n2*J2*n)/(n1*J1+n2*J2)
        s1(Mx,τx) = n21/n2
        //cas 2 : My= Mx-h et τy = τx
        M2σ = Mσ(Mx)-hm
        τ2σ = τσ(τx)
        M = ode("stiff",[M01;M02;800;800],0,t,tol, syst_dyn) 
        J1 = M(3,length(t))
        J2 = M(4,length(t))
        n11 = (n1*J1*n)/(n1*J1+n2*J2) 
        n21 = (n2*J2*n)/(n1*J1+n2*J2)
        s2(Mx,τx) = n21/n2
        //cas 3 : My= Mx et τy = τx+h
        M2σ = Mσ(Mx)
        τ2σ = τσ(τx)+hτ
        M = ode("stiff",[M01;M02;800;800],0,t,tol, syst_dyn) 
        J1 = M(3,length(t))
        J2 = M(4,length(t))
        n11 = (n1*J1*n)/(n1*J1+n2*J2) 
        n21 = (n2*J2*n)/(n1*J1+n2*J2)
        s3(Mx,τx) = n21/n2
        //cas 4 : My= Mx et τy = τx-h
        M2σ = Mσ(Mx)
        τ2σ = τσ(τx)-hτ
        M = ode("stiff",[M01;M02;800;800],0,t,tol, syst_dyn) 
        J1 = M(3,length(t))
        J2 = M(4,length(t))
        n11 = (n1*J1*n)/(n1*J1+n2*J2) 
        n21 = (n2*J2*n)/(n1*J1+n2*J2)
        s4(Mx,τx) = n21/n2

    end
end
disp("TIME = "+string(toc()/60)+" min")


for Mx = 1:length(Mσ)
    for τx = 1:length(τσ)
        dsdMy(Mx,τx) = (s1(Mx,τx)-s2(Mx,τx))/(2*hm) 
        dsdτy(Mx,τx) = (s3(Mx,τx)-s4(Mx,τx))/(2*hτ) 
    end
end

//affiche ds/dy = 0
scf() //on cherche coordonées du pt d'intersection des 2 courbes
// contour + interpolation
[xM,yM] = contour2di(Mσ,τσ,dsdMy,[-.01,0,.01]) //recupere les pts de contour
[xτ,yτ] = contour2di(Mσ,τσ,dsdτy,[-.01,0,.01])
plot(xM(2:yM(1)),yM(2:yM(1)),'b-')//seul le niveau 0 nous interesse
plot(xτ(2:yτ(1)),yτ(2:yτ(1)),'c-')
//title('$\left(\begin{array}{c}\left \frac{\partial s}{\partial M_y}\right| _{(My,τy)=(Mx,τx)} \\  \left \frac{\partial s}{\partial τ_y}\right| _{(My,τy)=(Mx,τx) \end{array}\right) =  \left(\begin{array}{c} 0 \\ 0\end{array}\right)$','fontsize',3)
xlabel("$M_r$",'fontsize',4)
ylabel("$\tau_r$",'fontsize',4)
legend(['$\frac{\partial s}{\partial M_m} = 0 $';'$ \frac{\partial s}{\partial τ_m} = 0$'],-1, %f);

xx = linspace(Mee-5,Mee+5,300)
yy = linspace(τee-0.07,τee+0.07,300)
yMitpl = interp1(xM,yM,xx)//interpolation des courbes autour de x*
yτitpl = interp1(xτ,yτ,xx)
M0 = interp1(yMitpl-yτitpl,xx,0) //en quel point elles se coupent 
τ0 = interp1(yMitpl-yτitpl,yy,0) // ? marche moyen ?
plot([M0 M0],[τσ(1) τ0],'k--')
plot([Mσ(1) M0],[τ0 τ0],'k--')
xstring(M0,τσ(1),'$M_r = '+string(M0)+'$')
xstring(Mσ(1),τ0,'$τ_r = '+string(τ0)+'$')
plot(Mee,τee,'k.')



//*******EQUATION CANONIQUE***********//

Mxd = Mσ
τxd = τσ
τx = 8:0.05:14
Mx = 1700:2100
[MX τX] =  ndgrid(Mx,τx)
C1 = splin2d(Mxd,τxd,dsdMy)//necessaire pour utiliser interp2d
C2 = splin2d(Mxd,τxd,dsdτy)


function [dy ] = canonique(t,y) 
    
    dy(1) = 10000*interp2d(y(1),y(2),Mxd,τxd,C1)
    dy(2) = 4*interp2d(y(1),y(2),Mxd,τxd,C2)
    
endfunction

scf()

x = 8:0.5:14
y = 1700:25:2100
fchamp(canonique,0,y,x)

temps = 0:10:10000
trajectoire = ode([2000;14],0,temps,10^-5,canonique) 
Mpointx = trajectoire(1,:)
τpointx = trajectoire(2,:)
plot(Mpointx,τpointx,"r-")

title('$canonical\ equation\ : \left(\begin{array}{c} \dot{M_x} \\ \dot{τ_x} \end{array}\right)\ = \left(\begin{array}{c} \frac{\partial s}{\partial M_y}(M_x,τ_x) \\ \frac{\partial s}{\partial τ_y}(M_x,τ_x) \end{array}\right)$','fontsize',3)
xlabel('$M_x$','fontsize',4)
ylabel('$τ_x$','fontsize',4)

xstring(1970,13,'$\dot{M_r} = 0 $')
gce().font_size = 3; 
gce().font_color = 2; 
xstring(1800,10.2,'$ \dot{τ_r} = 0$')
gce().font_size = 3; 
gce().font_color = 2; 


/*
for i = 1:length(Mσ)
    for j =  1:length(τσ)
        m(i,j) = dsdMy(i,j)^2+dsdτy(i,j)^2
    end
end
contour(Mσ,τσ,m,[-.01,0,.01])
f=scf()
f.color_map = rainbowcolormap(32)
surf(m)
[i,j] = find(m == min(m))
*/

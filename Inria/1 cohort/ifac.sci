clear
//clf()

//T < +∞

T = 20
t = 0:0.01:T

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\1 cohort\Fonctions 1 cohort.sci',-1)

Metoile = fsolve(0,eqnM)
Uetoile = 1-g(Metoile)./f(Metoile)//theorem 1

disp('M* = '+string(Metoile))

// ∃! τ*>0 such that λ(τ∗,M∗) = 0.
τetoile = fsolve(0,list(lambda,Metoile))

disp('τ* = '+string(τetoile))


M = 0:3000



//***************T infini******************************//

tolM = 10^-1/*
//Preparation du graphique
subplot(121)
v = 0:20
plot([0,length(v)],[Metoile,Metoile],'r--','LineWidth',2)
xstring(length(v), 2300, "M*");
gce().font_size=3
gce().font_color=5
xstring(length(v)/2, 1000, "u(M) = 0");
gce().font_size=3
xstring(length(v)/2, 4000, "u(M) = 1");
gce().font_size=3
xlabel('t - temps')
ylabel('M - mycélium')
//RESOLUTION ODE

ms1 = ode("stiff",[0.1;1],0,v,[tolM tolM], list(syst_dyn,"infini"))// ! attention oscillations ! tolérances !
//ms2 = ode("stiff",[5500;5500],0,v,[tolM tolM], list(syst_dyn,"infini"))
ms3 = ode("stiff",[100;100],0,v,[tolM tolM], list(syst_dyn,"infini"))
//ms4 = ode("stiff",[3500;3500],0,v,[tolM tolM], list(syst_dyn,"infini"))
ms5 = ode("stiff",[0.0001;0.0001],0,v,[tolM tolM], list(syst_dyn,"infini"))
plot2d(v,ms1(1,:))
//plot2d(v,ms2(1,:))
plot2d(v,ms3(1,:))
//plot2d(v,ms4(1,:))
plot2d(v,ms5(1,:))
xtitle("Optimal feedback control strategies for n = 10 and T = +∞")
*/

//***************T fini******************************//

//subplot(122)
plot([0,T-τetoile],[Metoile,Metoile],'r--','LineWidth',2)
plot([T-τetoile,T-τetoile],[0,Metoile],'k:')
plot2d([T+0.01,T+0.01],[0,5000],color("green"))
xstring(0, Metoile, "M*");
gce().font_size=3
gce().font_color=5
xstring(T, 0, "T");
gce().font_size=3
gce().font_color=3
xstring(T-τetoile/1.5, Metoile/2, "Γ");
gce().font_size=3
gce().font_color=2
xstring(T-τetoile, 0, "T-τ∗");
e = gce()
e.font_size=3
xstring(length(t)/4, 1000, "u(M) = 0");
gce().font_size=3
xstring(length(t)/4, 4000, "u(M) = 1");
gce().font_size=3
xlabel('t - temps')
ylabel('M - mycélium')
xtitle("Optimal feedback control strategies for n = 10 and T = 20")
//pour tracer sigma = switching surface
xres=[]
for i = 1:floor(Metoile)
    xres = [xres  fsolve(0,list(lambda,M(i)))]
end
plot2d(T-xres,M(1:floor(Metoile)),color('blue'))
//RESOLUTION ODE ! attention oscillations ! tolérances !
ms1 = ode("stiff",[1;1],0,t,[tolM tolM], list(syst_dyn,"fini"))
ms2 = ode("stiff",[5000;100],0,t,[tolM tolM], list(syst_dyn,"fini"))
ms3 = ode("stiff",[100;2000],0,t,[tolM tolM], list(syst_dyn,"fini"))
ms4 = ode("stiff",[3500;3500],0,t,[tolM tolM], list(syst_dyn,"fini"))
ms5 = ode("stiff",[0.01;0.01],0,t,[tolM tolM], list(syst_dyn,"fini"))
ms6 = ode("stiff",[0.0001;0.0001],0,t,[tolM tolM], list(syst_dyn,"fini"))
plot2d(t,ms1(1,:))
plot2d(t,ms2(1,:))
plot2d(t,ms3(1,:))
plot2d(t,ms4(1,:))
plot2d(t,ms5(1,:))
plot2d(t,ms6(1,:))

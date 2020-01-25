clear
clf()

//T < +∞

T = 20
t = 0:T

exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\1 cohort\Fonctions 1 cohort bis.sci',-1)

Metoile = fsolve(0,eqnM)
Uetoile = 1-g(Metoile)./f(Metoile)//theorem 1

disp('M* = '+string(Metoile))

// ∃! τ*>0 such that λ(τ∗,M∗) = 0.
τetoile = fsolve(0,list(lambda,Metoile))

disp('τ* = '+string(τetoile))


M = 0:5000



//***************T infini******************************//

tolM = 10^-1
//Preparation du graphique
subplot(121)
v = 0:20
plot([0,length(v)],[Metoile,Metoile],'r--','LineWidth',2)
xstring(length(v)-1, 2000, "$M^*$");
gce().font_size=3
gce().font_color=5
xstring(length(v)/2, 1000, "$u(M) = 0$");
gce().font_size=3
xstring(length(v)/2, 3000, "$u(M) = 1$");
gce().font_size=3
xlabel('t - temps','fontsize',4)
ylabel('M - mycélium','fontsize',4)
//RESOLUTION ODE

ms1 = ode("stiff",[0.1;1],0,0:0.1:11,[tolM tolM], list(syst_dyn,"infini"))// ! attention oscillations ! tolérances !
ms2 = ode("stiff",[4200;4200],0,0:0.1:12.9,[tolM tolM], list(syst_dyn,"infini"))
ms3 = ode("stiff",[100;100],0,0:0.1:5,[tolM tolM], list(syst_dyn,"infini"))
ms4 = ode("stiff",[3000;3000],0,0:0.1:7.8,[tolM tolM], list(syst_dyn,"infini"))
ms5 = ode("stiff",[0.0001;0.0001],0,0:0.1:17,[tolM tolM], list(syst_dyn,"infini"))
plot2d(0:0.1:11,ms1(1,:))
plot2d(0:0.1:12.9,ms2(1,:))
plot2d(0:0.1:5,ms3(1,:))
plot2d(0:0.1:7.8,ms4(1,:))
plot2d(0:0.1:17,ms5(1,:))
plot2d([5 20],[1920 1920])
title("T = +∞",'fontsize',5)


//***************T fini******************************//

subplot(122)
plot([0,T-τetoile],[Metoile,Metoile],'r--','LineWidth',2)
plot([T-τetoile,T-τetoile],[0,Metoile],'k:')
plot2d([T+0.01,T+0.01],[0,5000],color("green"))
xstring(0, Metoile, "$M^*$");
gce().font_size=3
gce().font_color=5
xstring(T, 0, "$T$");
gce().font_size=3
gce().font_color=3
xstring(14.5, 1000, "$Γ$");
gce().font_size=4
gce().font_color=2
xstring(T-τetoile, 100, "$T-τ^∗$");
e = gce()
e.font_size=3
xstring(7, 800, "$u(M) = 0$");
gce().font_size=3
xstring(7, 3000, "$u(M) = 1$");
gce().font_size=3
xlabel('t - temps','fontsize',4)
ylabel('M - mycélium','fontsize',4)
title(" T < +∞",'fontsize',5)
//pour tracer sigma = switching surface
xres=[]
for i = 1:floor(Metoile)
    xres = [xres  fsolve(0,list(lambda,M(i)))]
end
plot2d(T-xres,M(1:floor(Metoile)),color('blue'))
//RESOLUTION ODE ! attention oscillations ! tolérances !
t = 0:0.01:20
ms1 = ode("stiff",[1;1],0,0:0.1:9,[tolM tolM], list(syst_dyn,"fini"))
ms2 = ode("stiff",[4700;4700],0,t,[tolM tolM], list(syst_dyn,"fini"))
ms3 = ode("stiff",[100;100],0,0:0.1:4.5,[tolM tolM], list(syst_dyn,"fini"))
ms4 = ode("stiff",[2700;2700],0,0:0.1:5.9,[tolM tolM], list(syst_dyn,"fini"))
ms5 = ode("stiff",[0.00001;0.00001],0,t,[tolM tolM], list(syst_dyn,"fini"))
ms6 = ode("stiff",[0.000001;0.000001],0,t,[tolM tolM], list(syst_dyn,"fini"))
ms7 = ode("stiff",[1920;1920],0,10.85:0.1:20,[tolM tolM], list(syst_dyn,"fini"))
plot2d(0:0.1:9,ms1(1,:))
plot2d(t,ms2(1,:))
plot2d(0:0.1:4.5,ms3(1,:))
plot2d(0:0.1:5.9,ms4(1,:))
plot2d(t,ms5(1,:))
plot2d(t,ms6(1,:))
plot2d(10.85:0.1:20,ms7(1,:))
plot2d([4.3 10.85],[1920 1920])


exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tinfini\PIP - Tinfini.sce',-1)

f=scf()
f.color_map = rainbowcolormap(32);
surf(s)
xlabel("$M1σ$")
ylabel("$M2σ$")
zlabel("$s$")
title("$s(r,m) = \lim_{n_2(0) \to 0} \frac{n_2(1)}{n_2(0)}$",'fontsize',3)

scf() //affichage de s(x,x*) avec ds/dy=0
X = 1:length(Mσ)
snglr = find(Mσ == 1900)
plot(Mσ,s(snglr,X))
xlabel("$y$")
ylabel("$s$")
title("$s(x,x*)$",'fontsize',3)
ind = find(s(snglr,X)==max(s(snglr,X))) 
maxi = Mσ(ind)
plot([Mee,Mee],[s(snglr,1),s(snglr,ind)],'r--')
xstring(maxi-200,s(snglr,1),"Mmax = "+string(maxi))
plot([maxi,maxi],[s(snglr,1),s(snglr,ind)],'k--')
xstring(Mee,s(snglr,1),"M* = "+string(Mee))
gce().font_color = 5;
xstring(1900,s(snglr,snglr), "$\frac{\partial s}{\partial y} = 0$")
plot([1500 2200],[s(snglr,ind) s(snglr,ind)],'r--')


h = step
for x = 2:length(Mσ)-1
    gradselect(x-1) = (s(x,x+1)-s(x,x-1))/(2*h) //difference centrée
end

scf() //affichage de ds/dy(x,x)
y = 2:length(Mσ)-1
plot(Mσ(y),gradselect)
plot([Mσ(1),Mσ(length(y))],[0,0],'r--')
xarrows([1400:50:1800;1410:50:1850],zeros(2,9),400)
xarrows([2310:-50:1950;2300:-50:1920],zeros(2,8),400)
plot([Mee,Mee],[0,0],'k.')
plot([Mee,Mee],[gradselect(1),gradselect(length(y))],'r--')
xstring(Mee,gradselect(1),"M* = "+string(round(Mee)))
gce().font_color = 5;
Mmax = interp1(gradselect,Mσ(y),0)
plot([Mmax,Mmax],[gradselect(1),gradselect(length(y))],'k--')
xstring(Mmax-150,gradselect(length(y)),"M = "+string(round(Mmax)))
title("$gradient\ de\ selection : \frac{\partial s}{\partial y}(x,y=x) $",'fontsize',3)
xlabel("Mσ")

/*
title('$canonical\ equation\ :\ \frac{d}{dt}x = \frac{\partial s}{\partial y}(x,x)$','fontsize',3)
*/

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

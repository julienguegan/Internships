clear
exec('C:\Users\Julien Guégan\Desktop\PFE\Algorithmes\affichage.sce',-1)
function z = rozenbrock(x)
    z =  10*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction
function z = cout(x)
    z = 6*x(1).^2+2*x(2).^2+4*x(1).*x(2)+x(1)+x(2)
endfunction
 
function alpha = linearsearch1(f,x,d,grad)
    alpha = 1
    w1 = 0.9
    w2 = 0.1
    cpt = 0
    if(f(x+alpha*d)<(f(x)+alpha*w1*(grad'*d)))
        while(f(x+alpha*d)<(f(x)+alpha*w1*(grad'*d)))
            alpha = alpha*2 
            cpt = cpt+1
        end
    else 
        while(f(x+alpha*d)>(f(x)+alpha*w2*(grad'*d)))
            alpha = alpha/2 
            cpt = cpt+1
        end
    end
endfunction

function pas = linearsearch2(f,x,df)
    pas = 1
    grad = numderivative(f,)
    cpt=0
    while(f(x+pas*df)>(f(x)+0.5*grad'*df))
        cpt=cpt+1
        pas = (-pas*grad*df')/(2*(f(x+pas*df)-f(x)-pas*grad*df'))
    end
endfunction
clf(0)
tol = 0.001
f = cout
affiche(f,0.8,-2,2,-1.1,'contour')
g = gcf();
g.color_map = bonecolormap(50);
x0 =[-1.5 1.5]


///////////////pas variable//////////////////////////
x(1,:) = x0
df = numderivative(f,x0)
n = 1
while (norm(df'*df)>tol)
    pas = linearsearch2(f,x(n,:),df)
    //pas = -linearsearch1(f,x(n,:),-df,df)
    x(n+1,:) = x(n,:)+pas*df
    df = numderivative(f, x(n+1,:))
   // plot(x(n,1),x(n,2),'k.','markersize',3)
   // plot([x(n,1) x(n+1,1)],[x(n,2) x(n+1,2)],'k-')
    n = n+1;
    aploter(n) = f(x(n,:))
end
    for i=1:size(x,1)-1
        plot([x(i,1) x(i+1,1)],[x(i,2) x(i+1,2)],'b-')
    end
    
clf(1)
ite = 1:n
aploter(1) = f(x0)
plot(ite',aploter,'b','markersize',3)
////////////////gradient cpnjugué////////////////////////////
x = []
aploter = []
x(1,:) = x0+0.01
grad0 = numderivative(f,x0)
direction = grad0
k = 1
    while k<n then
        alpha = linearsearch2(f,x(k,:),direction)
        x(k+1,:) = x(k,:)+alpha*direction //le nouvel itéré solution du pb
        grad1 = numderivative(f,x(k+1,:)) //le nouveau gradient au point x(n+1)
        Beta = (grad1*grad1')/(grad0*grad0')// la constante Beta (un scalaire)
        direction = grad1+Beta*direction // la nouvelle direction
        grad0 = grad1//mise a jour pour la prochaine iteration
        k = k+1
        aploter(k) = f(x(k,:))
end
disp(k)
aploter(1) = f(x0)
ite = 1:k
plot(ite',aploter,'r-')
//////////////////////////////////////////////


///////////////pas = 0.01 //////////////////////////
x = []
aploter = []
x(1,:) = x0
df = numderivative(f,x0)
k = 1
while k<n
    pas = -0.05
    x(k+1,:) = x(k,:)+pas*df
    df = numderivative(f, x(k+1,:))
    k = k+1;
    aploter(k) = f(x(k,:))
end
aploter(1) = f(x0)
ite = 1:k
plot(ite',aploter,'b--')
//////////////////////////////////////////////////////


g = gca()
g.data_bounds=[0,0; n,f(x0)];
legend(["$pas\ variable$";;"$GC$";"$pas = 0.05$"] ,-1, %f)
xlabel('$nbr\ d\ évaluation$','fontsize',4)
ylabel('$f(x)$','fontsize',4)

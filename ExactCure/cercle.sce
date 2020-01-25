
function C = cercle(x,y,r)
// trace un cercle de centre (x,y) et de rayon r
    cpt1 = 1
    cpt2 = 1
    D = 2*r
    step = D/100
    for i = x-r:step:x+r
        ind_i = i/step-(x-r)/step+1
        C1(ind_i) = sqrt(r^2-(i-x)^2)+y
        C2(ind_i) = -sqrt(r^2-(i-x)^2)+y
        if imag(C1(ind_i)) == 0 then //on slct les nbr réels ! sqrt
            A(cpt1) = C1(ind_i) 
            cpt1 = cpt1 + 1 
        end
        if imag(C2(ind_i)) == 0 then
            B(cpt2) = C2(ind_i)
            cpt2 = cpt2 + 1 
        end
    end
    C = [A B]
    I = [x-r:step:x+r]
    plot(I',C,'b-','thickness',2)
    
endfunction

function E = ellipse(xC,yC,u,v,a)
// trace une ellipse de centre (xC,yC), de grand axe u et de petit axe v
//$(x-u)^2/a^2+(y-v)^2/b^2=1$
    cpt1 = 1
    cpt2 = 1
   D = 2*u
    step = D/100
    for i = xC-u:step:xC+u
        ind_i =i/step-(xC-u)/step+1
        
        E1(ind_i) = sqrt(v^2*(1-(i-xC)^2/u^2))+yC
        E2(ind_i) = -sqrt(v^2*(1-(i-xC)^2/u^2))+yC
        
        if imag(E1(ind_i)) == 0 then //on slct les nbr réels ! sqrt
            A(cpt1) = E1(ind_i) 
            cpt1 = cpt1 + 1 
        end
        if imag(E2(ind_i)) == 0 then
            B(cpt2) = E2(ind_i)
            cpt2 = cpt2 + 1 
        end
    end
    I = [xC-u:step:xC+u]
    E = [A B]
    
    x1 = I'.*cos(a)-A.*sin(a)
    x2 = I'.*cos(a)-B.*sin(a)
    y1 = I'.*sin(a)+A.*cos(a)
    y2 = I'.*sin(a)+B.*cos(a)
    X = [x1 x2]
    Y = [y1 y2]
    
    plot(X,Y,'b--')
endfunction

function y = ellipse2()

x = 0//upper left x
y = 0 //upper left y 
w = 1 // largeur
h = 0.5 //hauteur
a1 =  100//angle 1
a2 = 360*64 //angle2

arcs=[ x;y;w;h;a1;a2]
xarcs(arcs)
rand("normal")
for i = 1:100
    disp(rand(1,1))
end
endfunction


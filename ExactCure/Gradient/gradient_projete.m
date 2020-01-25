function [] = gradient_projete(f,n,x0)
%Méthode du gradient projeté
% gradient_projete(f,n,x0)
% f : Rn -> R, contrainte F sur x
% n, dimension du problème
% x0, point d'initialisation dans Rn
stop = 0;

if feval(f,x0)>0
    disp('Le point d initialisation doit etre dans le domaine F(x)<0');
    return
end

iter = 0;
itermax = 10000;
epsilon = 1e-6;
xk = x0;
while (stop == 0 && iter<itermax)
    
    iter = iter+1;
    
    % Vérification de la saturation de la contrainte
    fxk = feval(f,xk);
    sat = 0;
    if abs(fxk) <= epsilon
        sat = 1;
    end
    stopc = 0;
    while stopc == 0
        
        dgdx = grad(f,xk,n,epsilon);
        % Calcul de la matrice de projection
        P0 = eye(n);
        if sat == 1
            P0 = eye(n) - dgdx' * 1/(dgdx*dgdx') * dgdx;
        end
        yk = - P0 * xk / norm(xk);
        stopd = 0;
        if abs(yk) < epsilon
            stopd = 1;
        end
        if stopd ~= 1
            options = optimset('LargeScale','off');
            [alphamax,fval] = fmincon(@(alpha) -alpha,1,[],[],[],[],[], [],...
                @(alpha)constr(f,alpha,xk,yk), options);
            
            II
            if alphamax>epsilon
                [alphak,fxk] = fminbnd(@(alpha)norm(xk+alpha*yk),0,alphamax);
            else
                alphak = 0;
                fxk = 0;
            end
            
            % Calcul du nouveau point courant
            xk = xk+alphak*yk;
            stopc = 1;
            if norm(alphak*yk) < epsilon
                stop = 1;
                break
            end
        else
            u = - 1/(dgdx*dgdx')*dgdx*2*xk;
            if u >= -epsilon
                stopc = 1;
                stop = 1;
            else
                sat = 0;
            end
        end
    end
end
x = xk
disp('Solution x calculée par l algorithme de type gradient projeté');
x
disp('Nombre d iterations');
iter
disp('Norme de x');
norm(x)
end
function [c,ceq] = constr(f,alpha,xk,yk)
c = feval(f,xk+alpha*yk);
ceq = [];
end

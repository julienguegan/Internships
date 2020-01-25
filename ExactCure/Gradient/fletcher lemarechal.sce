clear all

exec('C:\Users\Julien Guégan\Desktop\PFE\fonctions test.sce',-1)
/*
* ne fonctionne pas ? donne alpha = 0.5
*/
function alpha = FL(f, x, d, grad)
    alpha_r = 10^10; //alpha_r tend vers l'infini
    alpha_l = 0;
    w1 = 0.1
    w2 = 0.9
    alpha = 1;
    cdt1 = %F //initialisation
    cdt2 = %F
    while ~(cdt1 && cdt2)
        cdt1 = numderivative(f,x + alpha * d) * d > w2*grad*d
        cdt2 = (f(x + alpha*d)) > (f(x) + alpha*w1*grad'*d)
        if (cdt1)
            alpha_r = alpha;
            alpha = (alpha_r + alpha_l) / 2;
        else
            if (cdt2)
                alpha_l = alpha;
                if (alpha_r < 10^10)
                    alpha = (alpha_l + alpha_r) / 2;
                else
                    alpha = 2 * alpha_l;
                end
            end
        end
    end
endfunction

tol = 0.001
itermax = 500
f = quadratique
dim = 2
x0 =[-0.9;0.9]//1*ones(dim,1)
df =  numderivative(f,x0)'
x(:,1) = x0
n = 1
tic()
while (tol<(df'*df))&(n<itermax) then
    df = numderivative(f,x(:,n))'

    //*************  FL **************//

    d = -df
    grad = df
    alpha_r = 10^10; //alpha_r tend vers l'infini
    alpha_l = 0;
    w1 = 0.1
    w2 = 0.9
    alpha = 1;
    cdt1 = %F //initialisation
    cdt2 = %F
    while ~(cdt1 && cdt2)
        cdt1 = numderivative(f,x(:,n) + alpha * d) * d > w2*grad'*d
        cdt2 = (f(x(:,n) + alpha*d)) > (f(x(:,n)) + alpha*w1*grad'*d)
        if (cdt1)
            alpha_r = alpha;
            alpha = (alpha_r + alpha_l) / 2;
        else
            if (cdt2)
                alpha_l = alpha;
                if (alpha_r < 10^10)
                    alpha = (alpha_l + alpha_r) / 2;
                else
                    alpha = 2 * alpha_l;
                end
            end
        end
    end 

    //*****************************************//
    x(:,n+1) = x(:,n) - alpha*df
    n = n+1
    disp('n = '+string(n) +' :  ,  alpha = '+string(alpha))
    disp(df)
    disp(x(:,n))
    plot(x(1,n),x(2,n),'k.')
end

time = toc()
disp('le minimum est x = ')
disp(x(:,n))
disp(' au bout de '+string(n)+' iterations, '+string(time)+' secondes')

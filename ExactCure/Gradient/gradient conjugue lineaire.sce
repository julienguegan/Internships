clear
clf()
exec('C:\Users\Julien Guégan\Desktop\PFE\fonctions test.sce',-1)

function [x,n] = gradientconjuguelineaire(A,b,tol,itermax,x0)
    x(:,1) = x0
    grad(:,1) = b - A*x0
    direction(:,1) = grad(:,1)
    n = 1

    while (tol<abs(grad(:,n)))&(n<itermax) then
        alpha(n) = (grad(:,n)'*grad(:,n))./(direction(:,n)'*A*direction(:,n))
        x(:,n+1) = x(:,n)+alpha(n)*direction(:,n)
        grad(:,n+1) = grad(:,n) - alpha(n)*A*direction(:,n)

        Beta(n) = (grad(:,n+1)'*grad(:,n+1))./(grad(:,n)'*grad(:,n))
        direction(:,n+1) = grad(:,n+1)+Beta(n)*direction(:,n)
        n = n+1
    end
endfunction

A = [4 2;2 2]
b = [-1;-1]
tol = [0.0001;0.0001]
itermax = 10
x0 = [0;0]

[x n] = gradientconjuguelineaire(A,b,tol,itermax,x0)
disp('le minimum est '+string(x(:,size(x,2)))+' au bout de '+string(n)+' iterations')

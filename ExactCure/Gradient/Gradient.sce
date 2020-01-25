
function alpha = FL(f, x, d, grad)
    alpha_r = 10^10; //alpha_r tend vers l'infini
    alpha_l = 0;
    w1 = 0.1
    w2 = 0.9
    alpha = 0.5;
    if ((f(x + alpha*d)) <= (f(x) + alpha*w1*grad'*d)) & (numderivative(f,x + alpha * d) * d > w2*grad*d)
        alpha = alpha;
    else
        if ((f(x + alpha*d)) > (f(x) + alpha*w1*grad'*d))
            alpha = (alpha_r + alpha_l) / 2;
        else
            alpha_l = alpha;
            if (alpha_r < 10^10)
                alpha = (alpha_l + alpha_r) / 2;
            else
                alpha = 2 * alpha_l;
            end
        end
    end
endfunction

function alpha = linesearch(f,x,d,grad)
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

function alpha = backtracking(f,x,d,grad)//pour newton mais pas pour quasi et CG
    alpha = 10
    w = 0.9
    cpt = 0
    while(f(x+alpha*d)>(f(x)+alpha*w*(grad'*d)))
        alpha = alpha/2 
        cpt = cpt+1
    end
endfunction


function [alpha,n] = goldstein(f,x,d,grad)
    alpha = 1
    w1 = 10^-4
    w2 = 1-w1
    cpt = 0
    while ~((f(x+alpha*d)<(f(x)+alpha*w1*(grad'*d)))&(f(x+alpha*d)>(f(x)+alpha*w2*(grad'*d))))
            alpha = alpha/2 
            cpt = cpt+1
    end
endfunction

function [x,n,stock] = gradientpasfixe(f,tol,itermax,x0)
    x = x0
    df =  numderivative(f,x0)'
    stock(:,1) = x0
    n = 1
    while (tol<(df'*df))&(n<itermax) then
        n = n+1
        df = numderivative(f,x)'
        //alpha = FL(f,x,-df,df)
        //alpha = backtracking(f,x,-df,df)
        alpha = linesearch(f,x,-df,df)//le pas
       // alpha = 0.02

        x = x - alpha*df
        stock(:,n)=x 

    end
endfunction

function [sol,n,x] = gradientconjugue(f,tol,itermax,x0,option)
    x(:,1) = x0
    grad0 = numderivative(f,x0)'// calcul gradient
    direction = grad0
    n = 1//nbr iteration

    while /*(tol<(grad0'*grad0))&*/(n<itermax) then
        alpha = linesearch(f,x,direction,grad0)//,'armijo')//le pas
        //alpha = search(f,x,direction,grad0)
        //alpha = backtracking(f,x,direction,grad0)
        // alpha = 0.0005
         
        x(:,n+1) = x(:,n)-alpha*direction //le nouvel itéré solution du pb
        grad1 = numderivative(f,x(:,n+1))' //le nouveau gradient au point x(n+1)
        if option=='FletcherReeves'
            Beta = (grad1'*grad1)/(grad0'*grad0)// la constante Beta (un scalaire)
        elseif option=='PolakRibiere'
            Beta = (grad1'*(grad1-grad0))/(grad0'*grad0)
        elseif option=='LiuStorey'
            Beta = -(grad1'*(grad1-grad0))/(direction'*grad0)
        elseif option=='DaiYuan'
            Beta = (grad0'*grad0)/(direction'*(grad1-grad0))
        elseif option=='HestenesStiefel'
            Beta = (grad0'*(grad1-grad0))/((grad1-grad0)'*direction)
        else 
            error('Veuillez rentrer une option valide')
            break
        end

        direction = grad1+Beta*direction // la nouvelle direction
        grad0 = grad1//mise a jour pour la prochaine iteration
        n = n+1
    end
    sol = x(:,n) //le dernier x est la solution
endfunction


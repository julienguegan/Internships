/*function alpha = linesearch(f,x,d,grad,methode)
    tau = 10^-2
    alpha = 1
    w1 = 10^-4
    w11 = 0.99
    w2 = 0.1
    i = 0
    condition = %F //initialisation condition

    while(~condition)
        cdtA = f(x+alpha*d)<(f(x)+alpha*w1*(grad'*d))
        cdtG = f(x+alpha*d)>(f(x)+alpha*w11*(grad'*d))
        // cdtW = (grad(x+alpha*d)'*d>w2*(grad'*d))

        if methode=='armijo' then
            condition = cdtA
        elseif methode=='goldstein' then
            condition = cdtA & cdtG
        elseif methode=='wolfe' then
            condition =cdtA & cdtW
        end
        alpha = alpha/2 
        i = i+1
    end
endfunction
*/
function alpha = linesearch(f,x,d,grad)
    // fonctionne pour tout alpha si w1=w2=0.9 ?
    // si w1=10^-4 w2 =0.99 ne converge pas toujours ?
    alpha = 1
    w1 = 0.1
    w2 = 0.9
    n = 0
    fx = f(x)
    if(f(x+alpha*d)<(fx+alpha*w1*(grad'*d)))
        while(f(x+alpha*d)<(fx+alpha*w1*(grad'*d)))
            alpha = alpha*2 
            n=n+1
            //disp(f(x+alpha*d))
        end
    else 
        while(f(x+alpha*d)>(fx+alpha*w2*(grad'*d)))
            alpha = alpha/2 
            n = n+1
            //disp(f(x+alpha*d))
        end
    end

endfunction


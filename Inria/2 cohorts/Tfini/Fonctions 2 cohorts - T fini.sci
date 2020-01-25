

//CONSTANTES BIOLOGIQUES
α = 0.2*10^4
k = (1/6)*10^4
β = 10^-5
γ = 0.06
δ = 1 //rendement de la production de spore en comparaison au mycelium
μ = 0.03

function y = f(i,m1,m2) //mycelial growth
    if (i == 1) then
        ρ = α*m1./(m1 + k)
    elseif (i == 2) then
        ρ = α*m2./(m2 + k)
    else 
        error("parameter i sould be 1 ou 2")
    end
    ν = 1./(1+β*(n1*m1+n2*m2))

    y = ν.*ρ
endfunction
function y = f(i,m1,m2,n) //mycelial growth
    if (i == 1) then
        ρ = α*m1./(m1 + k)
    elseif (i == 2) then
        ρ = α*m2./(m2 + k)
    else 
        error("parameter i sould be 1 ou 2")
    end
    ν = 1./(1+β*(n*m1))
    y = ν.*ρ
endfunction

function y = g(M) //mycelial decay
    y = γ*M
endfunction

// M** verifie f'(M**)-g'(M**) = μ  => ν(nM**)ρ'(M**)-g'(M**)=μ
function eqn = eqnM(M)
    eqn = α*k/((1+β*n*M)*(M+k)^2)-μ-γ
endfunction
function eqn = eqnM(M,n)
    eqn = α*k/((1+β*n*M).*(M+k).^2)-μ-γ
endfunction
function i = termintegral(s,M)
    i = exp(-(γ+μ)*s)*α*k/((1+β*n*M*exp(-γ*s))*(M*exp(-γ*s)+k)^2)
endfunction
function l = lambda(tau,m)
    l = 1-intg(0,tau,list(termintegral,m),0.001)
endfunction

//commande u
function U = u(t,M,i)
    if (i == 1) then
        mσ = M1σ
        Ue = 1+(/*0.01*(M1σ-M1)*/-g(M1))./f(1,M1,M2) 
        τe = τ1σ
    else
        mσ = M2σ
        Ue = 1+(/*0.01*(M2σ-M2)*/-g(M2))./f(2,M1,M2) 
        τe = τ2σ
    end
    if (t < T-τe) then
        if (M >= mσ) then
            U = Ue
        else 
            U = 0
        end
    else 
        U = 1
    end

endfunction

function dMJ = syst_dyn(t,M)  
    M1 = M(1)
    M2 = M(2)
    dMJ(1) = (1 - u(t,M1,1)).*f(1,M1,M2) - g(M1) //dM1
    dMJ(2) = (1 - u(t,M2,2)).*f(2,M1,M2) - g(M2) //dM2

    dMJ(3) = u(t,M1,1).*f(1,M1,M2).*exp(-μ*t) //dJ1
    dMJ(4) = u(t,M2,2).*f(2,M1,M2).*exp(-μ*t) //dJ2

endfunction

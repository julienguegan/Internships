

//CONSTANTES BIOLOGIQUES
α = 0.2*10^4
k = (1/6)*10^4
β = 10^-5
n = 10//nombre de mycelium d'un groupe
γ = 0.06
δ = 1//rendement de la production de spore en comparaison au mycelium
μ = 0.03

//croissance du mycélium limitée par la capacité d'absorption du champignon
function f = f(M,n)
    f1 = α.*M./(M+k) //flux de ressources pouvant être obtenu par 1 mycélium
    f2 = 1./(1+β*n*M) //influence négative de compétition entre mycélium sur un même hôte
    f = f1.*f2
endfunction

//mortalité du mycélium
function g = g(M)
    g = γ.*M
endfunction

// M* s.t f'(M*)-g'(M*) = μ  
function eqn = eqnM(M,n)
    eqn = (α*k-α*β*n*M.^2-(γ+μ)*(1+β*n*M).^2.*(M+k).^2)./((1+β*n*M).^2.*(M+k).^2)
endfunction

// ∃! τ*>0 such that λ(τ∗,M∗) = 0.
function y = fprime(M)
    y = α*(-β*n*M^2+k)./((1+β*n*M).^2.*(M+k).^2)
endfunction
function i = termintegral(s,M)
    i = exp(-(γ+μ)*s).*fprime(M.*exp(-γ*s))
endfunction
function l = lambda(tau,m)
    l = -intg(0,tau,list(termintegral,m))+1
endfunction

//commande u d'apres Pontryagin (theorem 4)
function u_opti = u(t,M,cas) //cas = "fini" ou "infini"
    u_opti = [] 
    if (cas == "fini") then
            if ((M < Metoile) & (lambda(T-t,M) >= 0)) then
                u_opti =  1
            elseif ((M < Metoile) & (lambda(T-t,M) < 0)) then 
                u_opti = 0
            elseif (M == Metoile) then
                u_opti = Uetoile
            else 
                u_opti = 1
            end
    elseif (cas == "infini") then
            if (M < Metoile) then
                u_opti = 0
            elseif (M == Metoile) then 
                u_opti =  Uetoile
            else
                u_opti =  1
            end
    end
endfunction

//systeme dynamique mycelium/spores
function dMS = syst_dyn(t,MS,cas) 
    M = MS(1)//mycelium
    S = MS(2)//spore
    dMS(1) = (1-u(t,M,cas)).*f(M)-g(M)
    dMS(2) = δ*u(t,M,cas).*f(M)
endfunction


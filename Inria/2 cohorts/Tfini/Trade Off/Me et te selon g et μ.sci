clear

α = 0.2*10^4
k = (1/6)*10^4
β = 10^-5
δ = 1
n = 101

function f = f(M)
    f1 = α.*M./(M+k)
    f2 = 1./(1+β*n*M)
    f = f1.*f2
endfunction
function g = g(M,γ)
    g = γ.*M
endfunction
function eqn = eqnM1(M,γ,μ)//1 cohort
    eqn = (α*k-α*β*n*M.^2-(γ+μ)*(1+β*n*M).^2.*(M+k).^2)./((1+β*n*M).^2.*(M+k).^2)
endfunction
function eqn = eqnM2(M,γ,μ)//2 cohort
    eqn = α*k/((1+β*n*M).*(M+k).^2)-μ-γ
endfunction
function i = termintegral(s,M,γ,μ)
    i = exp(-(γ+μ)*s)*α*k/((1+β*n*M*exp(-γ*s))*(M*exp(-γ*s)+k)^2)
endfunction
function l = lambda(tau,m,γ,μ)
    l = 1-intg(0,tau,list(termintegral,m,γ,μ))
endfunction
/*
scf()
i=0
for γ = 0.01:0.005:0.1
    i=i+1
    μ = 0.03
    Metoileγ(i) = fsolve(0,list(eqnM2,γ,μ))
    τetoileγ(i) = fsolve(0,list(lambda,Metoileγ(i),γ,μ))
end
j=0
for μ =0.01:0.005:0.1
    j=j+1
    γ = 0.06
    Metoileμ(j) = fsolve(0,list(eqnM2,γ,μ))
    τetoileμ(j) = fsolve(0,list(lambda,Metoileμ(j),γ,μ))
end
γ = 0.01:0.005:0.1
μ = 0.01:0.005:0.1
subplot(221)
plot2d(γ ,τetoileγ,color("blue"))
xlabel('$γ$')
legend("$τ*$",-1, %f)
xgrid()
subplot(222)
plot2d(γ ,Metoileγ,color("red"))
xlabel('$γ$')
legend("$M*$",-1, %f)
xgrid()
subplot(223)
plot2d(μ ,τetoileμ,color("blue"))
xlabel('$μ$')
legend("$τ*$",-1, %f)
xgrid()
subplot(224)
plot2d(μ ,Metoileμ,color("red"))
xlabel('$μ$')
legend("$M*$",-1, %f)
xgrid()

*/
scf()
j=0
μ = 0.01
while (μ < 0.1)
    γ = 0.01
    while  (γ < 0.1)
        j=j+1
        Metoileμγ(j) = fsolve(0,list(eqnM2,γ,μ))
        τetoileμγ(j) = fsolve(0,list(lambda,Metoileμγ(j),γ,μ))
        
        γ = γ*1.2//γ*1.15
    end
    μ= μ*1.2//μ*1.15
end
plot(Metoileμγ,τetoileμγ,'k.')
ylabel('$τ^{**}$','fontsize',4)
xlabel('$M^{**}$','fontsize',4)



//MOINDRES CARRES - quadratique
function y = model(t, x)
    y  = x(1)*t^2+x(2)*t+x(3) //forme quadratique
endfunction

y = τetoileμγ
t = Metoileμγ

// we want to find the parameters x such that 
//  minimize  $f(x) = \sum_{i=1}^N   ( f(t_i,x) - y_i )^2$
function e = aminimiser(x,t, y) 
    e =  model(t, x) - y
endfunction
x0 = [1;1;1];//CI
[f,xopt] = leastsq(list(aminimiser,t,y),x0)

tradeoff_M = linspace(1000,4000,3001)';
tradeoff_τ = model(tradeoff_M, xopt);
plot(tradeoff_M,tradeoff_τ,'b-','thickness',3)
legend(["$0.01<γ<0.1\\ 0.01<μ<0.1$";"$f(x)$"] ,-1, %f)
title('$Moindres\ Carrés:minimiser\ f(x) = \sum_{i=1}^N   ( f(t_i,x) - y_i )^2$')

save("C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tfini\Trade-Off\tradeoff_save",'tradeoff_τ','tradeoff_M','Metoileμγ','τetoileμγ')

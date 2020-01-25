clear

μ = 0.03
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
function eqn = eqnM1(M,γ,α)//1 cohort
    eqn = (α*k-α*β*n*M.^2-(γ+μ)*(1+β*n*M).^2.*(M+k).^2)./((1+β*n*M).^2.*(M+k).^2)
endfunction
function eqn = eqnM2(M,γ,α)//2 cohort
    eqn = α*k/((1+β*n*M).*(M+k).^2)-μ-γ
endfunction
function i = termintegral(s,M,γ,α)
    i = exp(-(γ+μ)*s)*α*k/((1+β*n*M*exp(-γ*s))*(M*exp(-γ*s)+k)^2)
endfunction
function l = lambda(tau,m,γ,α)
    l = 1-intg(0,tau,list(termintegral,m,γ,μ))
endfunction

scf()
i=0
for γ = 0.0001:0.001:0.1
    i=i+1
    α = 2*10^3
    Metoileγ(i) = fsolve(0,list(eqnM2,γ,α))
    τetoileγ(i) = fsolve(0,list(lambda,Metoileγ(i),γ,α))
end
j=0
for α = 500:50:3500
    j=j+1
    γ = 0.06
    Metoileμ(j) = fsolve(0,list(eqnM2,γ,α))
    τetoileμ(j) = fsolve(0,list(lambda,Metoileμ(j),γ,α))
end
γ = 0.0001:0.001:0.1
α = 500:50:3500
subplot(221)
plot2d(γ ,τetoileγ,color("blue"))
xlabel('$γ$','fontsize',3)
legend("$τ*$",-1, %f)
xgrid()
subplot(222)
plot2d(γ ,Metoileγ,color("red"))
xlabel('$γ$','fontsize',3)
legend("$M*$",-1, %f)
xgrid()
subplot(223)
plot2d(α,τetoileμ,color("blue"))
xlabel('$α$','fontsize',3)
legend("$τ*$",-1, %f)
xgrid()
subplot(224)
plot2d(α,Metoileμ,color("red"))
xlabel('$α$','fontsize',3)
legend("$M*$",-1, %f)
xgrid()


scf()
j=0
for α =  1000:1000:10000 //
    
   for  γ = 0.0001:0.0025:0.1
       
        j=j+1
        Metoileαγ(j) = fsolve(0,list(eqnM2,γ,α))
        τetoileαγ(j) = fsolve(0,list(lambda,Metoileαγ(j),γ,α))
       plot(Metoileαγ(j),τetoileαγ(j),'k.')
   end
end

//plot(Metoileαγ,τetoileαγ,'k.')
ylabel('$τ*$','fontsize',4)
xlabel('$M*$','fontsize',4)



//MOINDRES CARRES - quadratique
function y = model(t, x)
   y  = x(1)*t+x(2) //forme quadratique
endfunction

y = τetoileαγ
t = Metoileαγ

// we want to find the parameters x such that 
//  minimize  $f(x) = \sum_{i=1}^N   ( f(t_i,x) - y_i )^2$
function e = aminimiser(x,t, y) 
   e =  model(t, x) - y
endfunction
x0 = [1;1];//CI
[f,xopt] = leastsq(list(aminimiser,t,y),x0)

tradeoff_M = linspace(1000,7000,5001)';
tradeoff_τ = model(tradeoff_M, xopt);
plot(tradeoff_M,tradeoff_τ,  'b-','thickness',3)
legend(["$0.01<γ<0.1\\ 1000<α<10000$";"$f(x)$"] ,-1, %f)
title('$Moindres\ Carrés:minimiser\ f(x) = \sum_{i=1}^N   ( f(t_i,x) - y_i )^2$')

save("C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tfini\Trade-Off\tradeoff_save2",'tradeoff_τ','tradeoff_M','Metoileαγ','τetoileαγ')

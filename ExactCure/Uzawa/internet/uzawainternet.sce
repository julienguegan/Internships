clear

function [x,i,resi]=U1(A,b,C,f,eta,Imax,rho)

    x=zeros(b);
    r=A*x-b;
    nr=norm(r);
    i=0;
    err=eta;
    resi=[];
    lambda=zeros(n,1);
    while (i<Imax & norm(C*x-f)>eta) 
        xold=x;
        i=i+1;
        x=inv(A)*(b-C'*lambda);
        lambda = lambda+rho*(C*x-f);
        r=A*x-b;
        resi=[resi;err];
        if nr>1e10; disp('Explosion!'); break; end;
    end
    if (err<eta); disp('Methode terminee normalement'); end
    if (i==Imax); disp('Imax atteint! Augmenter le nombre d''iterations'); end;
endfunction

n=29;
h=1/(n+1);
n=29;
A=zeros(n,n);
A=zeros(28,28);
A(n,n)=2;
A(n,n-1)=-1;
for i=1:n-1 do
    A(i,i)=2;
    A(i,i+1)=-1;
    A(i+1,i)=-1;
end

b=-20*ones(n,1);

eta=1e-4;
Imax=10000;
S=spec(A);
rho=1000;

t=zeros(n,1);
for i=1:n do
    t(i)=i*h;
end

g=zeros(n,1);
for i=1:n do
    g(i)=-1+max(0,0.565-10*(t(i)-0.4)^2);
end

C=-eye(n,n);
f=-g;

disp('rho :');disp(rho);

[x,i,resi]=U1(A,b,C,f,eta,Imax,rho);
clf;
xgrid(3);
plot2d((1:i),resi,logflag='nl', style=[5]);
xtitle('nombre d itérations = '+string(i));
disp('dernier residu :');disp(resi(i));

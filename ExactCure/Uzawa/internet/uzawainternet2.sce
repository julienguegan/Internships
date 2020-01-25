N=4; 
h = (1/(N-1)); 
A = zeros(N,N); 

for i = 1:N 
for j = 1:N 

if (j==i) then 
A(i,j)=2; 
end 
if (j== i+1) then A(i,j)=-1; end 
if (j== i-1) then A(i,j)=-1; end 

end; 
end; 
disp ('afficher la matrice', A); 
A = (1/h)*A; 
F=h*ones(N,1) 
B=zeros(2,N); 
B(1,1)=1; 
B(2,2)=1; 
B(1,N)=-1; 
B(2,N-1)=-1; 


T=zeros(1,N+2); 


k=zeros(N+2,N+2); 
k=[A,B';B,[0 0;0 0]] 
T=[F;0;0] 
S=linsolve(k,-T) 

rho=10^(-3); 
eps=10^(-6); 
u=rand(zeros(N,1)); 
p=rand(zeros(2,1)); 
nmax=10000 
function[p1,u1,n]=Uzawa(u,p,A,B,F,nmax,eps,rho) 
p1=p; 
for n=1:nmax 
d=F-B'*p1; 
u=inv (A)*d; 
p2=p1+rho*B*u; 
if (norm(p2-p1))<eps then 
return 
end 
p1=p2 
end 
u1=u 
endfunction 

[p1,u1,n]=Uzawa(u,p,A,B,F,nmax,eps,rho); 


disp('la solution 1 est :'); disp(p1)
disp('la solution 2 est '); disp(u1)

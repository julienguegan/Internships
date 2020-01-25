function y=f(x)
    y=n+sum(x.^2-cos(2*%pi*x));
//    y=10*(x(2)-x(1)^2)^2+(x(1)-1)^2;
endfunction
////////parametres/////////////////
Npop=50;// pair
Ngen=100;
n=2;
pc=0.6// probabilité de croisement
pm=0.9// probabilité de mutation
R=0.1  // rayon de la mutation
////////////////////////////////////
//
A=zeros(Npop,n+1); // matrice de population 
// derniere colonne: valeur de la fonction
//
// initialisation aléatoire sur [-5,5]^n
A=10*rand(Npop,n+1)-5*ones(Npop,n+1);
//
for i=1:Ngen
// 
//evaluation
y=[];
for j=1:Npop
    y(j)=f(A(j,1:n));
end  
A(:,n+1)=y; 
//selection
[u,v]=gsort(A(:,n+1));
A=A(v,:); // rangement du plus mauvais au meilleur
//
disp('meilleur'),disp(A(Npop,1:n))
p=(1:Npop)/sum(1:Npop);ps=cumsum(p);
Asel=[]
for i=1:Npop
  u=rand();isel=1;
  while (u>ps(isel)) 
     isel=isel+1;
  end
Asel=[Asel;A(isel,:)];
end  
A=Asel;
//
//croisement
for k=1:Npop/2
  n1=int(Npop*rand())+1;
  n2=int(Npop*rand())+1;
  alpha=rand();
  u1=A(n1,:);u2=A(n2,:);
  if(rand()<pc) then
   A(n1,:)=alpha*u1+(1-alpha)*u2;
   A(n2,:)=(1-alpha)*u1+alpha*u2;
  end
end
// mutation
for k=1:Npop
  if(rand()<pm) then
   A(k,1:n)=A(k,1:n)-R*ones(1,n)+2*R*rand(1,n);   
  end    
end
//
clf()
plot2d(A(:,1),A(:,2),-1,rect=[-5,-5,5,5])
//pause
//
end
////////////////
//



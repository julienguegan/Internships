clear

function y = frozenbrock(x)
    if size(x,1) < 2
        error('dimension must be greater one');
    end
    N = size(x,1)
    y = 100*sum((x(1:N-1).^2 - x(2:N)).^2) + sum((x(1:N-1)-1).^2);
endfunction

function y = frastrigin(x)
    if size(x,1) < 2
        error('dimension must be greater one');
    end
    N = size(x,1)
    y = 10*N + sum(x(1:N).^2-10*cos(2*%pi*x(1:N)));
endfunction

f = frastrigin
x0 = 5*rand(10,1,'uniform')+20
s0 = 0.25
itermax = 1500
tol = 0.00001
lambda = 30
mu = 10


// --------------------  Initialization --------------------------------  
// User defined input parameters (need to be edited)

N = length(x0);               // dimension
xmean = x0;                   // initial point
sigma = s0;                   // step size

// Strategy parameter setting: Selection  
weights = log(mu+1/2)-log(1:mu)';        // muXone array for weighted recombination       
weights = weights/sum(weights);          // normalize recombination weights array
mueff = sum(weights)^2/sum(weights.^2);    // variance-effectiveness of sum w_i x_i

// Strategy parameter setting: Adaptation
cc = (4+mueff/N) / (N+4 + 2*mueff/N);                       // time constant for cumulation for C
cs = (mueff+2) / (N+mueff+5);                               // t-const for cumulation for sigma control
c1 = 2 / ((N+1.3)^2+mueff);                                 // learning rate for rank-one update of C
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));   // and for rank-mu update
damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs;         // damping for sigma usually close to 1

// Initialize dynamic (internal) strategy parameters and constants
pc = zeros(N,1); ps = zeros(N,1);   // evolution paths for C and sigma
B = eye(N,N);                       // B defines the coordinate system
D = ones(N,1);                      // diagonal D defines the scaling
C = B * diag(D.^2) * B';            // covariance matrix C
invsqrtC = B * diag(D.^-1) * B';    // C^-1/2 
eigeneval = 0;                      // track update of B and D
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  // expectation of ||N(0,I)|| == norm(randn(N,1)) 

// -------------------- Generation Loop --------------------------------
counteval = 0;  // the next 40 lines contain the 20 lines of interesting code 
rand("normal")
cpt = 0

while counteval < itermax

    // Generate and evaluate lambda offspring
    for k = 1:lambda
        x(:,k) = xmean + sigma * B * (D .* rand(N,1)); // m + sig * Normal(0,C) 
        fitness(k) = f(x(:,k)); // objective function call
        counteval = counteval+1;
    end
    [fitness, index] = gsort(fitness,'g','i'); // minimization
    xold = xmean;
    xmean = x(:,index(1:mu))*weights;
    ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma; 
    hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4 + 2/(N+1);
    pc = (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;

    // Adapt covariance matrix C
    artmp = (1/sigma) * (x(:,index(1:mu))-repmat(xold,1,mu));
    C = (1-c1-cmu) * C ...                  // regard old matrix  
    + c1 * (pc*pc' ...                      // plus rank one update
    + (1-hsig) * cc*(2-cc) * C) ...         // minor correction if hsig==0
    + cmu * artmp * diag(weights) * artmp'; // plus rank mu update

    // Adapt step size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); 

    if counteval - eigeneval > lambda/(c1+cmu)/N/10 
        eigeneval = counteval;
        C = triu(C) + triu(C,1)';  
        [B,D] = spec(C);
        D = sqrt(diag(D));                          
        invsqrtC = B * diag(D.^-1) * B';
    end

    if fitness(1) <= tol || max(D) > 1e7 * min(D)
        break;
    end
    
    cpt = cpt+1
    aploter(cpt) = f(xmean)
end 

xmin = x(:, index(1));

clf()
g = gca()
g.data_bounds=[0,-10; counteval,f(x0)];
xlabel('$nbr\ d\ évaluation$','fontsize',4)
ylabel('$f(x)$','fontsize',4)


nbeval=lambda*[1:cpt]
plot(nbeval,aploter','b.','markersize',3)
ycmi = smooth([nbeval;aploter'],0.01);
plot(ycmi(1,:)',ycmi(2,:)','b-');





function  [xm,n,aploter] = evolutionstrategie(f,xm0,sm0,lambda,mu,itermax,tol) // μ < λ
    
    dim = length(xm0)
    tau = 1 //facteur d'apprentissage, petit = sigma plus petit
    n = 0 //nbr iteration
    xmp1 = xm0
    xm = xm0
    smp1 = sm0
    rand("normal")
    counteval = 0
    while counteval < itermax

        xm = xmp1
        sm = smp1

        for i = 1:lambda //creation
            s(i,:) = sm .* exp(tau*rand(1,dim))
            x(i,:) = xm + s(i).*rand(1,dim)
        end

        Z = zeros(lambda)
        for i = 1:lambda
            Z(i) = f(x(i,:)')
            counteval = counteval + 1
        end

        xp = zeros(mu,dim) //vecteur des x parents
        sp = zeros(mu,dim) // les ecarts types associés

        for j = 1:mu //selection 
            mini =  find(Z == min(Z))// indice de minf(x) 
            xp(j,:) = x(mini,:)//vecteur des mu minimum
            sp(j,:) = s(mini,:)
            Z(mini,:) = []//on l'enleve pour la prochaine iteration
            x(mini,:) = []
            s(mini,:) = []
        end

        xmp1 = sum(xp,'r')/mu //recombinaison
        smp1 = sum(sp,'r')/mu //autoadaptation
            n = n+1
         aploter(n) = f(xm')
    end

endfunction

[sol n ploter] = evolutionstrategie(f,x0',s0,lambda,mu,itermax,tol)
nbeval=lambda*[1:n]
plot(nbeval,ploter','r.','markersize',3)
yi = smooth([nbeval;ploter'],0.01);
plot(yi(1,:)',yi(2,:)','r-');

legend(["$CMAES$";"$(\mu,\lambda)-ES$"] ,-1, %f)

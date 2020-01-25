clear


function y=frozenbrock(x)
    if size(x,1) < 2
        error('dimension must be greater one');
    end
    N = size(x,1)
    y = 100*sum((x(1:N-1).^2 - x(2:N)).^2) + sum((x(1:N-1)-1).^2);
endfunction



// --------------------  Initialization --------------------------------  

N = 20               
xmean = rand(N,1)                  
sigma = 0.3                   
itermax = 1e3*N^2
tol = 1e-10

// Strategy parameter setting: Selection  
lambda = 4+floor(3*log(N)); 
mu = lambda/2;    
weights = log(mu+1/2)-log(1:mu)';          
mu = floor(mu);
weights = weights/sum(weights);          
mueff = sum(weights)^2/sum(weights.^2);    


// Strategy parameter setting: Adaptation
cc = (4+mueff/N) / (N+4 + 2*mueff/N);                     
cs = (mueff+2) / (N+mueff+5);                             
c1 = 2 / ((N+1.3)^2+mueff);                         
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));   
damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs;        

// Initialize dynamic (internal) strategy parameters and constants
pc = zeros(N,1); ps = zeros(N,1); 
B = eye(N,N);                     
D = ones(N,1);                   
C = B * diag(D.^2) * B';        
invsqrtC = B * diag(D.^-1) * B';    
eigeneval = 0;                    
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  


counteval = 0;  
rand("normal")
while counteval < itermax
    // Generate and evaluate lambda offspring
    for k = 1:lambda
        x(:,k) = xmean + sigma * B * (D .* rand(N,1));  
        fitness(k) = frozenbrock(x(:,k)); 
        counteval = counteval+1;
    end
    // Sort by fitness and compute weighted mean into xmean
    [fitness, index] = gsort(fitness,'g','i'); 
    xold = xmean;
    xmean = x(:,index(1:mu))*weights;   
    // Cumulation: Update evolution paths
    ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma; 
    hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4 + 2/(N+1);
    pc = (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;
    // Adapt covariance matrix C
    artmp = (1/sigma) * (x(:,index(1:mu))-repmat(xold,1,mu));
    C = (1-c1-cmu) * C+ c1 * (pc*pc'+ (1-hsig) * cc*(2-cc) * C)+ cmu * artmp * diag(weights) * artmp'; 
    //Adapt step size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); 
    //Diagonalisation de C
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

end 

xmin = x(:, index(1));



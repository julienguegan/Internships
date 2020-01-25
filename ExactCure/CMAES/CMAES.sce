clear

exec('C:\Users\Julien Guégan\Desktop\PFE\algorithmes\fonctions test.sce',-1)

function [xmin,cpt] = purecmaes(f, x0, s0, lambda, mu, itermax, tol)//CMA-ES
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

    // ********* calcul F *********//
    xmax=1.5;xmin=-1.5;ymax=1.5;ymin=-1.5;
    xstep = (xmax-xmin)/95
    axe_x = xmin:xstep:xmax
    ystep = (ymax-ymin)/95
    axe_y = ymin:ystep:ymax
    z=0
    for i = xmin:xstep:xmax
        for j = ymin:ystep:ymax
            u = [i j]
            ind_i = (i/xstep)-(xmin/xstep)+1
            ind_j = (j/ystep)-(ymin/ystep)+1
            z(ind_i,ind_j) = f(u)
        end
    end
    ////////////////////////////////////
    cpt = 0
    xold = [3;3]
    while cpt < itermax & norm(xmean-xold)>tol
        //****  AFFICHAGE  ***//
        clf()
        xset("fpf"," ")
        contour2d(axe_x,axe_y,z,30)
        a = gcf()
        a.color_map = rainbowcolormap(64)
        xlabel('$x_1$','fontsize',4)
        ylabel('$x_2$','fontsize',4)
        //*******************************************************//


        // Generate and evaluate lambda offspring
        for k = 1:lambda
            x(:,k) = xmean + sigma * B * (D .* rand(N,1)); // m + sig * Normal(0,C) 
            fitness(k) = f(x(:,k)); // objective function call
            counteval = counteval+1;
               if (xmin<x(1,k)&x(1,k)<xmax)&(ymin<x(2,k)&x(2,k)<ymax)
                plot(x(1,k),x(2,k),'k.','markersize',3)
            end
        end
        //plot(x(1,:),x(2,:),'k.','markersize',2)
        // Sort by fitness and compute weighted mean into xmean
        [fitness, index] = gsort(fitness,'g','i'); // minimization
        xold = xmean;
        xmean = x(:,index(1:mu))*weights;   // recombination, new mean value

        // Cumulation: Update evolution paths
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

        // Decomposition of C into B*diag(D.^2)*B' (diagonalization)
        if counteval - eigeneval > lambda/(c1+cmu)/N/10  // to achieve O(N^2)
            eigeneval = counteval;
            C = triu(C) + triu(C,1)';                   // enforce symmetry
            [B,D] = spec(C);//eigs()                             // eigen decomposition, B==normalized eigenvectors
            D = sqrt(diag(D));                          // D is a vector of standard deviations now
            invsqrtC = B * diag(D.^-1) * B';
        end
        // Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable 
        /*if fitness(1) <= tol || max(D) > 1e7 * min(D)
            break;
        end*/
        plot(xmean(1),xmean(2),'b*')
        cpt = cpt+1
        xs2png(gcf(),"cmaes"+string(cpt))
    end // while, end generation loop

    xmin = x(:, index(1)); // Return best point of last iteration.
    // Notice that xmean is expected to be even
    // better.
endfunction

f = rozenbrock
x0 = [0;0]
s0 = 0.1
lambda = 30
mu = 10
itermax = 100
tol = 0.001
[sol n] =  purecmaes(f,x0,s0,lambda,mu,itermax,tol)
disp('x = ')
disp(sol)
disp('n = ')
disp(n)

/*
cpti=0
for i=1:10//plusieurs iterations car stochastiques 
    cptmoy = 0
    cptx0 = 0
    for x0 = [[0;0],[1;1],[5;5],[10;10]]//differents initial guess
        cptevalmoy = 0
        cptx0 = cptx0+1
        for s0 = [0.01,0.1,1,5,10,100]//differents pas
            cpts0 = 0
            cpteval = 0
            for lambda = [1,10,20,30,40,50]//differents parametres lambda
                for mu = [1,5,10,20,30]//differents parametres mu
                    if(mu<lambda)
                        cpts0 = cpts0+1
                        //timer()
                        [sol counteval] = purecmaes(f,x0,s0,lambda,mu,itermax,tol)
                        //disp('l='+string(lambda)+'m='+string(mu)+' : eval = '+string(counteval)+' , time = '+string(timer()))
                        cpteval = cpteval + counteval
                    end
                end
            end
            cptevalmoy = cptevalmoy + cpteval/cpt
            //disp('s0='+string(s0)+' : '+string(cpteval/cpt)+' evaluations en moyenne')
        end
        //disp('x0=('+string(x0(1))+','+string(x0(2))+') : '+string(cptevalmoy/cptx0)+' evaluations en moyenne')
    end

end


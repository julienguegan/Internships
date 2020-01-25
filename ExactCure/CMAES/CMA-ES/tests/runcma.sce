// getf('fitfuns.sci'); // for frosen, felli, etc, should be in the demos dir
if 1 < 3 | ~isdef('felli')
  // getf
  exec(get_absolute_file_path('runcma.sce') + ...
       'fitfuns.sci'); // felli, frosen, etc
end
// getf('cma.sci', 'p'); // p for profiling, plotprofile(cma_tell)

stacksize(30e6); // cma_plot might need a large stack 

N = 10; 
fitfun = frosen; // fsphere felli frastrigin10 
clear p 
p.x0 = '0.1 + 0.0*rand(N,1)-0';  // initial solution point, is evaluated in constructor
p.sigma0 = 0.3;                  // initial standard deviation on x0 
p.opt.scaling_of_variables = []; 2*ones(N,1); // zero excludes variables from actual search space

// p.stop.fitness = 1e-9;            // user defined

p.verb.logmodulo = 10;            // writing of output files every logmodulo iteration
p.verb.plotmodulo = '2e4/lambda'; // uses cma_plot() and output files
p.verb.displaymodulo = 1e2;       // screen display 

p.opt.separable = 0; '30*N*10/lambda'; // in particular for high dimension 

date0 = getdate(); timer(); 

number_of_runs = 1; // > 1 restarts, only in case the target fitness was not reached
for irun = 1:number_of_runs

  es = cma_new(p);  // initialize

  // iteration loop
  while ~ cma_stop(es)
    X = cma_ask(es);   // returns a list with lambda column vectors
    y = [];             // just in case
    for i = 1:length(X) // evaluate all lambda candidate solution
      y(i) = fitfun(X(i)); // must return a scalar
      // re-sampling of infeasible solution if f == %nan (rejection method)
      while isnan(y(i)) 
	X(i) = cma_ask(es, 1); 
	disp(' solution ' + string(i) + ' re-sampled');
	y(i) = fitfun(X(i));
      end

    end // for each new candidate solution
    // do not use modified X here, unless you very well know what you do
    es = cma_tell(es, X, y'); // finish iteration
  end

  // output
  for str = es.out.stopflags
    printf('terminated on ''' + str + '''\n');
  end
  [ybest, xbest] = cma_best(es); 
  printf('  x(1)=%f   y=%f', xbest(1), ybest);

  disp('elapsed CPU time: ' + string(timer()) + ' s, ' ...
      + 'elapsed time: ' + string(etime(getdate(),date0)) + ' s');

  // manage restart with increasing population size, IPOP-, G-CMA-ES
  // does not add up number of f-evals
  if irun < number_of_runs & es.out.stopflags(1) ~= 'fitness' 
    disp('restart');
    p.opt.lambda = 2*es.param.opt.lambda;
    p.verb.append = es.counteval;
    // p.verb.filenameprefix = es.param.verb.filenameprefix + string(irun+1);
  else
    break; // incomment for several runs anyway
  end

end // for irun


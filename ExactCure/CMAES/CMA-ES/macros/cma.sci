// FUNCTIONS IN THIS FILE
//external functions, interface
//  function [es, p] = cma_new(inparam)  // constructor, p is optional the parameter struct
//  function data    = cma_plot(fignb, name_prefix, name_extension, object_variables_name, plotrange)
//  function X       = cma_ask(es, lam)  // delivers lam new points
//  function es      = cma_tell(es, X, arfitness) // update, concludes iteration
//  function flag    = cma_stop(es)     // termination criterion met?
//  function [y, x]  = cma_best(es)     // best-ever solution, see also es.out.solutions
//  function X       = %cma_ask(es)      // delivers const es.param.opt.lambda new points
//  function es      = %cma_tell(es, X, arfitness) // update, concludes iteration
//  function [y, x]  = %cma_best(es)     // best-ever solution, see also es.out.solutions

//potentially useful functions
//  function param   = check_and_merge_param(inparam, defparam) // verifies and merges two parameter lists
//  function y       = cma_genophenotransform(gp, x)   // gp=es.out.genopheno
//  function x       = cma_phenogenotransform(gp, x)   // transform into genotypic
//  function res     = cma_minnan(vec)                 // utility, compute min disregarding %nan
//internal functions
//  function p       = cma_getdefaults()
//  function gp      = cma_tlistgenopheno(typical_x, scales, x0)
//  function scales  = cma_check_initial_sizes(typical_x, x0, scales)
//  function x0      = cma_xstart(typical_x, x0, scales)
//  function           cma_annotate(xfinal, yfinal, ann) // not in use yet
//  function           cma_plotintern(name_prefix, fignb, name_extension, object_variables_name)

// see http://www.cert.fr/dcsd/idco/perso/Magni/s_sym/overload.html
// TODO:
//       o replace D completely with a vector diagD, then
//       o make separable with linear *space* complexity B,D,C
//       o (not trivial) implement cma_ask_array(es,
//         form=['col'|'row') and cma_tell_array(es, form='col') by
//         moving most code of cma_tell() into _cma_tell() leaving
//         %cma_tell(),
//         cma_tell() and cma_tell_array() merely as interface functions.
//       o constraints handling should be copied from cmaes.sci (should be easy)
//       o cma_annotate should be used to allow plotting of the genotype if need be
//       o xml docu: help_skeleton('funnam', '.'); xmltohtml(path, name);
//

// --------------------------------------------
function param = cma(dummy)
//This is a dummy for the cma library/toolbox
//that needs to be referenced once to make all funtions available.
//
// See functions
//
//      cma_new
//      cma_optim
//
// help for getting started.
//
  param = cma_getdefaults();
endfunction

// --------------------------------------------
function p = cma_getdefaults()
// returns struct param
// see function cma_new
//
  p = struct(); // the below might be slightly faster, but produces ugly output

//  p = tlist(['param', 'x0', 'typical_x', 'sigma0', 'opt_scaling_of_variables', ...
//      'stop_tolfun', 'stop_tolfunhist', 'stop_tolx', 'stop_tolupx', ...
//      'stop_fitness', 'stop_maxfunevals', 'stop_maxiter', 'opt_lambda', ...
//      'opt_mu', 'opt_separable', 'verb_logfunvals', 'verb_logmodulo', ...
//      'verb_displaymodulo', 'verb_plotmodulo', 'verb_readfromfile', ...
//      'verb_filenameprefix', 'readme']);

// reason for using stop_... instead of a struct: cma_new() displays everything, otherwise assignment is needed

  p.x0 = [];
  p.typical_x = [];
  p.sigma0 = [];
  p.opt.scaling_of_variables = [];

  p.stop.tolfun = 1e-12;
  p.stop.tolfunhist = 1e-12;
  p.stop.tolx = '1e-11*sigma0';
  p.stop.tolupx = '1e3*sigma0';
  p.stop.fitness = -%inf;
  p.stop.maxfunevals = %inf;
  p.stop.maxiter = '100*N+(N+5)^2*1e4/lambda'; // should depend on restarts as well

  p.opt.lambda = '4+floor(3*log(N))';
  p.opt.mu = 'lambda/2';
  p.opt.separable = 0; '5*N*10/lambda'; // not yet default, but in near future?

  // p.verb.silent = 0; // turn off all verbosity
  p.verb.logmodulo  = 1;        // writing output data
  p.verb.displaymodulo = 100;   // display
  p.verb.plotmodulo = 'max(logmodulo, 500)';   // only if output data are written
  p.verb.logfunvals = [];
  p.verb.readfromfile = 'signals.par'; // not really verbosity
  p.verb.filenameprefix = 'outcmaes';    //
  p.verb.append = 0             // overwrite old output files

  p.readme = [];
  p.readme.x0 = 'initial solution, either x0 or typical_x MUST be provided';
  p.readme.typical_x = 'typical solution, the genotypic zero for more convenient plots in future';
  p.readme.sigma0 = 'initial coordinate-wise standard deviation MUST be provided';
  p.readme.stop = 'termination options, see .stop.readme';
  p.readme.opt = 'more optional parameters, see .opt.readme';
  p.readme.verb = 'verbosity parameters, see .verb.readme';
  p.stop.readme.tolfun = 'stop if within-iteration function value differences are smaller';
  p.stop.readme.tolfunhist = 'stop if function value backward differences are smaller';
  p.stop.readme.tolx = 'stop if x-changes are smaller than tolx*sigma0*scaling';
                               // internally scaling must be omitted as it is part of GP-TF
  p.stop.readme.tolupx = 'stop if x-changes are larger than tolupx*sigma0*scaling';
  p.stop.readme.fitness = 'target objective function value (minimization)';
  p.stop.readme.maxfunevals = 'maximal number of function evaluations';
  p.stop.readme.maxiter = 'maximal number of iterations';
  p.opt.readme.scaling_of_variables = ...
      ['sigma0*scaling is in effect the initial standard deviation.' ...
       + ' scaling == 0 means that the variable does not take part in the optimization'];
  p.opt.readme.lambda = 'population size, number of sampled candidate solutions per iteration';
  p.opt.readme.mu = 'parent number, for experts only. Weighted recombination and mu=lambda/2 are default';
  p.opt.readme.separable = ['number of initial iterations with independent sampling for partly' ...
      + ' separable high-dimensional functions (n>100), 1==always'];
  // p.verb.readme.silent = 'if >0, turn off all verbosity';
  p.verb.readme.logmodulo = 'data logging, write every logmodulo iteration to file';
  p.verb.readme.displaymodulo = 'display every displaymodulo iteration some output to console';
  p.verb.readme.plotmodulo = 'call cma_plot() every plotmodulo iteration, only if logmodulo>0';
  p.verb.readme.logfunvals = 'record whenever these f-values are reached first';
  p.verb.readme.readfromfile = 'filename, ""stop"" as first line leads to manual termination';
  p.verb.readme.filenameprefix = 'output files';
  p.verb.readme.append = 'if true, start with Fevals=append and append files'

endfunction

// --------------------------------------------
function flag = cma_stop(es)
  flag = ~isempty(es.out.stopflags);
endfunction

// --------------------------------------------
function [y, x] = cma_best(es)
  x = es.out.solutions.bestever.x;
  y = es.out.solutions.bestever.f;
endfunction

// --------------------------------------------
function [y, x] = %cma_best(es)
   [y, x] = cma_best(es);
endfunction

// --------------------------------------------
function X = %cma_ask(es)
// returns a list with lambda N-dimensional candidate solutions
   X = cma_ask(es, es.sp.lambda);
endfunction

// --------------------------------------------
function X = cma_ask(es, lambda, mod)
// returns a list with lambda N-dimensional candidate solutions
// cma_tell only accepts a list of param.opts.lambda
// solutions, but parameter lambda is useful to resample
// single solution(s) conveniently. For lambda==1 one
// solution (not a list) is returned. With mod==1 always
// a list is returned.

  if ~isdef('lambda', 'local') | isempty(lambda) | lambda == 0
    lambda = es.sp.lambda
  end
  if ~isdef('mod', 'local') | isempty(mod)
    mod = 0;
  end

  // X = [];
  X = list();
  for k=1:lambda
    if es.countiter <= es.sp.separable // diagonal CMA
      x = es.xmean + es.sigma * (diag(es.D) .* rand(es.N, 1, 'normal'));
    else
      x = es.xmean + es.sigma * (es.BD * rand(es.N, 1, 'normal'));
    end
    // X(:,k) = cma_genophenotransform(es.out.genopheno, x);
    X(k) = cma_genophenotransform(es.out.genopheno, x);
  end
  if lambda == 1 & mod == 0
    X = X(1);
  end
endfunction

////////////////////////////////////////////////////////////
// ------------------- CONSTRUCTOR -------------------------
////////////////////////////////////////////////////////////
function [es, p] = cma_new(inparam)

  // string 'cma' is the object type, mlist because of ??
  es = tlist(['cma', 'stop', 'param', 'out', 'inparam', 'defparam', ...
          'N', 'sp', 'countiter', 'counteval', ...
          'sigma', 'xmean', 'pc', 'ps', 'B', 'D', 'BD', 'C', 'fitnesshist', 'const', ...
          'verb', 'files', 'logfunvals']);
        // list must contain all variables

  es.defparam = cma_getdefaults();
  p = [];

  // return default parameters
  if ~isdef('inparam','local')
    es = es.defparam;
    printf('cma_new has two mandatory fields in its input parameter struct: ' ...
        +'\n   x0 (or typical_x) and sigma0. ' ...
        +'\nA complete parameter struct has been returned.\n')
    return;
  end
  if (typeof(inparam) == 'string' & inparam == 'getdefaults') ...
      | isempty(inparam)
    es = es.defparam;
    return;
  end

  es.inparam = inparam;
  p = check_and_merge_param(inparam, es.defparam);

  // keep empty readfromfile
  if isfield(inparam, 'verb') & isfield(inparam.verb, 'readfromfile')
    p.verb.readfromfile = inparam.verb.readfromfile;
  end

  es.out = [];
  es.out.seed = rand('seed'); // needs to be done before any input parameter is evaluated
  es.out.version = 0.998;


// -------------------- Handle Input Arguments -------------------------------

  // mandatory arguments
  p.x0 = evstr(p('x0'));
    x0 = p.x0; // for historical reasons, might be removed in near future
  p.typical_x = evstr(p('typical_x'));
    typical_x = p.typical_x; // for historical reasons, might be removed in near future
  if (isempty(x0) & isempty(typical_x))
    error('CMA-ES: x0 or typical_x must be defined');
  end

  // set dimension
  p.opt.scaling_of_variables = evstr(p.opt.scaling_of_variables);
  scaling_of_variables = cma_check_initial_sizes(typical_x, x0, p.opt.scaling_of_variables);
  es.N = sum(scaling_of_variables > 0);
  N = es.N; // for evstr()

  // set sigma
  p.sigma0 = evstr(p('sigma0')); // todo: remove this?
  if isempty(p.sigma0)
    error('sigma0 (initial one-sigma) must be defined');
  end
  if ~and(size(p.sigma0) == [1 1])
    error('input parameter sigma must be a scalar, use scaling_of_variables otherwise');
  end
  es.sigma = p.sigma0;
  sigma0 = p.sigma0; // for evstr()
  sigma = p.sigma0; // for evstr()

  //
  // TODO: check: sigma tolx stopping criteria are consistent with settings of sigma
  //       and of the gp-transformation?
  es.out.genopheno = list(); // otherwise the next assignment produces an error
  es.out.genopheno = cma_tlistgenopheno(typical_x, scaling_of_variables, x0);

  // ------ we need to know N and sigma from here on (not xmean)

  // --- optional arguments

  p.opt.lambda = evstr(p.opt.lambda);
  es.sp.lambda = p.opt.lambda;
  lambda = p.opt.lambda; // for evstr()
  if (lambda < 3) // even for 2-D it is unreliable
    error('param.opt.lambda must be larger than 2');
  end
  if (lambda < 4 & N > 4)
    warning('param.opt.lambda set too small (unreliable setting)');
  end

  p.opt.mu = evstr(p.opt.mu);
  es.sp.mu = p.opt.mu; // integer is not enforced, 1:mu is used, better should round up for trunc(mu)>0.5
  if ceil(es.sp.mu) > lambda
    error('param.opt.mu must be smaller than param.opt.lambda');
  end

  // verbosity
  logmodulo = p.verb.logmodulo; // for evstr()
  if p.verb.logmodulo == 0
    p.verb.plotmodulo = 0;
  end
  if ~isempty(p.verb.logfunvals) // array of decreasing function values when to record data
    p.verb.logfunvals = gsort(p.verb.logfunvals);
  end
  es.logfunvals = p.verb.logfunvals; // dynamic list, written elements are removed
  p.verb.plotmodulo = evstr(p.verb.plotmodulo);
  if 11 < 3 // | p.verb.silent
    p.verb.logmodulo = 0;
    p.verb.displaymodulo = 0;
    p.verb.plotmodulo = 0;
    if p.verb.silent > 0
      es.logfunvals = [];
    end
  end

  // termination
  p.stop.maxfunevals = evstr(p.stop.maxfunevals);
  p.stop.maxiter = evstr(p.stop.maxiter);
  p.stop.tolx = evstr(p.stop.tolx);
  p.stop.tolupx = evstr(p.stop.tolupx);
  p.stop.tolfun = evstr(p.stop.tolfun);
  p.stop.tolfunhist = evstr(p.stop.tolfunhist);

  ccovfac = 1;
  p.opt.separable = evstr(p.opt.separable); // iteration(N, lambda)
  es.sp.separable = p.opt.separable; // iteration(N, lambda)
  if es.sp.separable ~= 0
    ccovfac = (N+1.5)/3;
    if es.sp.separable == 1
      es.sp.separable = %inf;
    end
  end

  // -------------------- Initialization --------------------------------
  es.counteval = 0;
  es.out.dimension = N;
  // es.out.fitnessfunction = fitfun; // makes the struct non-investigable
  // es.out.stopflagss = list(); // for restarts

  // set xmean and sigma
  // TODO: check this
  // es.sigma = in.sigma; (obsolete? in case of restart?)
  clear sigma; // to be on the save side
  clear lambda; // TODO: clear all input parameters needed for getparam up to here for savity
  es.xmean = cma_phenogenotransform(es.out.genopheno, cma_xstart(typical_x, x0, scaling_of_variables));


  // --- Strategy parameter setting: Selection
  es.sp.weights = log(es.sp.mu+0.5)-log(1:es.sp.mu)'; // muXone array for weighted recombination
  // es.sp.lambda=12; es.sp.mu=3; weights = ones(es.sp.mu,1); disp('equal w''s'); // uncomment for (3_I,12)-ES
  es.sp.weights = es.sp.weights/sum(es.sp.weights);           // normalize recombination weights array
  es.sp.mueff = sum(es.sp.weights)^2/sum(es.sp.weights.^2); // variance-effective size of mu

  // Strategy parameter setting: Adaptation
  es.sp.cc = (4 + es.sp.mueff/N) / (N+4 + 2*es.sp.mueff/N);  // time const. for cumulation for covariance matrix
  if es.sp.separable > 2*sqrt(N)  // noticable difference only for N>>100
    es.sp.cc = (1+1/N + es.sp.mueff/N) / (N^0.5 + 1/N + 2.*es.sp.mueff/N)
//    es.sp.cc = (1+1/N + es.sp.mueff/N) / (N^0.5 + 1./N + 2.*es.sp.mueff/N)
  end
  es.sp.cs = (es.sp.mueff+2)/(N+es.sp.mueff+3);  // time const. for cumulation for sigma control
//  es.sp.cs = sqrt(es.sp.mueff)/(sqrt(N-1)+sqrt(es.sp.mueff)); disp('cs=' + string(es.sp.cs));
//  es.sp.cs    = (es.sp.mueff+2)/(sqrt(N)+es.sp.mueff+3); disp('cs=' + string(es.sp.cs));

  es.sp.mucov = es.sp.mueff;            // mu used for C update, mucov=1 ==> rank-1 only
//  es.sp.mucov = 1; disp('mucov=' + string(es.sp.mucov));
  es.sp.ccov_def = (1/es.sp.mucov) * 2/(N+1.4)^2 ...     // learning rate for covariance matrix
          + (1-1/es.sp.mucov) * ((2*es.sp.mucov-1)/((N+2)^2+2*es.sp.mucov));
  es.sp.ccov = ccovfac * es.sp.ccov_def;
  if es.sp.ccov > 1
    es.sp.ccov = 1;
    if ccovfac == 1
      warning('ccov set to one, due to a BUG? ');
    end
  end
  es.sp.damps = 1 + es.sp.cs ...                      // damping for sigma, usually close to 1
          + 2*max(0, sqrt((es.sp.mueff-1)/(N+1))-1);


// qqq hack different constants here
//  es.sp.cc = 1; disp('cc=' + string(es.sp.cc));
//  es.sp.damps = 0.5*es.sp.damps; 1e97; disp('damps=' + string(es.sp.damps)) // "turn off" sigma adaptation
 //es.sp.ccov = 0.0*es.sp.ccov; disp('es.sp.ccov=' + string(es.sp.ccov)); // turn off covariance matrix adaptation

  // Initialize dynamic (internal) strategy parameters and constants
  es.pc   = zeros(N,1);                   // evolution paths for C
  es.ps = zeros(N,1);                     //   and sigma
  es.B    = eye(N,N);                     // B defines the coordinate system
  es.D    = eye(N,N);                     // diagonal matrix D defines the scaling
  es.BD   = eye(N,N); // == es.B*es.D;       selling memory for speed
  es.C    = eye(N,N); // == es.BD*(es.BD)'   covariance matrix
  es.const.chiN = N^0.5*(1-1/(4*N)+1/(21*N^2)); // expectation of ||N(0,I)|| == norm(randn(N,1))
                                       // exact: sqrt(2) * gamma((N+1)/2) / gamma(N/2)
  es.fitnesshist = %nan*ones(10+ceil(30*N/es.sp.lambda),1); // history of fitness values

  // initialize constraints
//  if ~isdef('constraints_weights', 'local')
//    constraints_weights = [];
//  end
//  constraints_obj = constraints_init(N, es.sp.lambda, constraints_weights, ceil(2+N/1));

  // initial printed message
  if (p.verb.displaymodulo)
    printf('(%d/%d_W,%d)-CMA-ES ', es.sp.mu, es.sp.mu, es.sp.lambda);
    printf('(W=[');
    for i=1:min(3,length(es.sp.weights))
      printf('%2.0f,', 100*es.sp.weights(i));
    end
    printf('...]%%, mueff=%3.1f) in %d', es.sp.mueff, N);
    if N < length(es.out.genopheno.scaling)
      printf(' of %d', length(es.out.genopheno.scaling));
    end
    printf('-D\n');
    if es.sp.separable ~= 0
      printf('Covariance matrix is diagonal');
      if es.sp.separable > 1 & es.sp.separable < %inf
        printf(' for %d iterations (%d evaluations)', ...
            es.sp.separable, es.sp.separable*es.sp.lambda);
      end
      printf('!\n');
    end

    if p.verb.plotmodulo == 0 & p.verb.logmodulo
      printf('    for plotting call CMA_PLOT using Ctrl-C or another Scilab shell\n');
    end
    es.verb.dispannotationline = ...
        'Iter, Evals: Function Value   (worst) |Axis Ratio |idx:Min SD, idx:Max SD\n';
  end


  // Write headers to output data files
  if p.verb.logmodulo  | es.logfunvals
    es.files.additions = ...
        ['axlen', 'fit', 'stddev', 'xmean', 'xrecentbest'];
    //          ['axlen', 'disp', 'fit', 'xrecentbest', 'stddev', 'xmean'];
    if p.verb.append
      es.counteval = p.verb.append;
    else
      for name = es.files.additions
        [fid, err] = mopen(p.verb.filenameprefix+name+'.xml', 'w');
        if err ~= 0
          warning('could not open '+p.verb.filenameprefix+name+'.xml' + ' ');
          es.files.additions(find(es.files.additions == name)) = [];
        else
          mfprintf(fid, '%s\n', ...
              '<CMAES-OUTPUT version=""' + string(es.out.version) + '"">');
          mfprintf(fid, '  <NAME>'+name+'</NAME>\n');
          mfprintf(fid, '  <DATE>'+date()+'</DATE>\n');
          mfprintf(fid, '  <PARAMETERS>\n');
          mfprintf(fid, '    dimension=' + string(N) + '\n');
          mfprintf(fid, '  </PARAMETERS>\n');
          // different cases for DATA columns annotations here
          mfprintf(fid, '  <DATA');
          if name == 'axlen'
             mfprintf(fid, '  columns=""iteration, evaluation, sigma, ' + ...
                 'max axis length, min axis length, ' + ...
                 'all principle axes lengths (sorted square roots of eigenvalues of C)""');
          elseif name == 'fit'
            mfprintf(fid, '  columns=""iteration, evaluation, sigma, axis ratio, bestever,' + ...
                ' best, median, worst fitness function value,' + ...
                ' further objective values of best""');
          elseif name == 'stddev'
            mfprintf(fid, '  columns=""iteration, evaluation, sigma, void, void, ' + ...
                'stds==sigma*sqrt(diag(C))""');
          elseif name == 'xmean'
            mfprintf(fid, '  columns=""iteration, evaluation, void, void, void, xmean""');
          elseif name == 'xrecentbest'
            mfprintf(fid, '  columns=""iteration, evaluation, fitness, void, void, xrecentbest""');
          end
          mfprintf(fid, '>\n'); // DATA
          if name == 'xmean'
            mfprintf(fid, '%ld %ld 0 0 0 ', 0, es.counteval);
            // mfprintf(fid, '%ld %ld 0 0 %e ', es.countiter, es.counteval, fmean);
//qqq            mfprintf(fid, msprintf('%e ', genophenotransform(es.out.genopheno, es.xmean)) + '\n');
            mfprintf(fid, msprintf('%e ', es.xmean) + '\n');
          end
          mclose(fid);
        end
      end // for files
    end
  end

  es.param = p;

  es.countiter = 0; // within the restart loop
  es.stop = %F;
  es.out.stopflags = list();

endfunction

// --------------------------------------------
function es = %cma_tell(es, X, arfitness)
  // X is a list of N-dimensional vectors
  //   (was a N x lambda array, a row vector of column vectors of size N)
  // arfitness is a row vector
  // both are transposed if the size is unique
  // this is not compliant to the OMD specification(s) (array of row vectors)

  es = cma_tell(es,X,arfitness);

endfunction

// --------------------------------------------
function es = cma_tell(es, X, arfitness)
  // X is a list of N-dimensional vectors
  //   (was a N x lambda array, a row vector of column vectors of size N)
  // arfitness is a row vector
  // both are transposed if the size is unique
  // this is not compliant to the OMD specification(s) (array of row vectors)

  N = es.N;

  // this is a bit unsave, as errors might occur later just by changing the problem dimension
  if size(arfitness, 1) == es.sp.lambda  & size(arfitness, 2) ~= es.sp.lambda
    arfitness = arfitness';
  end
  if 11 < 3 & size(X, 1) ~= N & size(X, 2) == N
    // X = X';
  end

  if length(X) ~= es.sp.lambda  // size(X, 1) ~= N | size(X, 2) ~= es.sp.lambda
    error('cma_tell: input X has wrong size, lambda vectors are expected');
  end
  if size(arfitness, 2) ~= es.sp.lambda
    error('cma_tell: arfitness must be a row vector of size lambda');
  end

  es.countiter = es.countiter + 1;
  es.counteval = es.counteval + es.sp.lambda;
  arx = [];

  for k=1:es.sp.lambda
    // arx(:,k) = cma_phenogenotransform(es.out.genopheno, X(:,k));
    if length(X(k)) ~= N
      error('vectors must have size ' + string(N))
    end
    arx(:,k) = cma_phenogenotransform(es.out.genopheno, X(k));
    if es.countiter <= es.sp.separable // diagonal CMA
      arz(:,k) = (arx(:,k)-es.xmean)./diag(es.D)/es.sigma; // == sigma^-1*D^-1*B'*(xmean-xold)
    else
      arz(:,k) = es.sigma^(-1)*diag(diag(es.D).^(-1))*es.B'*(arx(:,k)-es.xmean);   // == sigma^-1*D^-1*B'*(xmean-xold)
    end
  end

  // if or(abs(X(:,1) - es.xmean) > 10 * es.sigma * sqrt(diag(es.C)))
  if or(abs(arz) > 10)  // P(10-sigma) < 1e-20
    warning('input X is not consistent and will presumably lead to failure to converge');
  end

    // Sort by fitness and compute weighted mean into xmean
    [arfitness, arindex] = gsort(-arfitness); // minimization
    // this Matlab incompatibility sheds a dubious light on Scilab
    arfitness = - arfitness;
    es.fitnesshist = [arfitness(1); es.fitnesshist(1:$-1)];

    es.xmean = arx(:,arindex(1:es.sp.mu))*es.sp.weights;   // recombination, new mean value
    zmean = arz(:,arindex(1:es.sp.mu))*es.sp.weights;   // == sigma^-1*D^-1*B'*(xmean-xold)

    // Cumulation: Update evolution paths
    if es.countiter <= es.sp.separable // diagonal CMA
      es.ps   = (1-es.sp.cs)*es.ps + sqrt(es.sp.cs*(2-es.sp.cs)*es.sp.mueff) * (1 * zmean); // Eq. (4)
    else
      es.ps   = (1-es.sp.cs)*es.ps + sqrt(es.sp.cs*(2-es.sp.cs)*es.sp.mueff) * (es.B * zmean); // Eq. (4)
    end
    hsig = norm(es.ps) / sqrt(1-(1-es.sp.cs)^(2*es.counteval/es.sp.lambda))/es.const.chiN < 1.4 + 2/(N+1);
//qqq
//hsig = 1; disp('hsig=' + string(hsig));
    if es.countiter <= es.sp.separable // diagonal CMA
      es.pc   = (1-es.sp.cc)*es.pc + hsig * sqrt(es.sp.cc*(2-es.sp.cc)*es.sp.mueff) * (diag(es.D) .* zmean); // Eq. (2)
    else
      es.pc   = (1-es.sp.cc)*es.pc + hsig * sqrt(es.sp.cc*(2-es.sp.cc)*es.sp.mueff) * (es.BD * zmean); // Eq. (2)
    end

    // Adapt covariance matrix C
    if 1 < 3 & es.countiter <= es.sp.separable // diagonal CMA
      es.C = diag((1-es.sp.ccov) * diag(es.C) ...
          + es.sp.ccov * (1/es.sp.mucov) * (es.pc.^2 ...   // plus rank one update
                                            + (1-hsig) * es.sp.cc*(2-es.sp.cc) * diag(es.C)) ...
          + es.sp.ccov * (1-1/es.sp.mucov) ...           // plus rank mu update
                       * diag(es.C) .* (arz(:,arindex(1:es.sp.mu)).^2 * es.sp.weights));
    else
      es.C = (1-es.sp.ccov) * es.C ...                   // regard old matrix        // Eq. (3)
          + es.sp.ccov * (1/es.sp.mucov) * (es.pc*es.pc' ...   // plus rank one update
                                            + (1-hsig) * es.sp.cc*(2-es.sp.cc) * es.C) ...
          + es.sp.ccov * (1-1/es.sp.mucov) ...           // plus rank mu update
                       * (es.BD*arz(:,arindex(1:es.sp.mu))) ...
                       *  diag(es.sp.weights) * (es.BD*arz(:,arindex(1:es.sp.mu)))';
    end
    // Adapt step size sigma
    es.sigma = es.sigma * exp((es.sp.cs/es.sp.damps)*(norm(es.ps)/es.const.chiN - 1));  // Eq. (5)

    // Update B and D from C (without the modulo this is O(N^3))
    if es.countiter <= es.sp.separable // diagonal CMA
      es.C = diag(es.C); // for efficiency C is a diagonal vector for the next lines

      // Align order of magnitude of scales of sigma and C for nicer output
      if 1 < 2 & es.sigma > 1e10*sqrt(max((es.C)))
        fac = es.sigma/sqrt(median((es.C)));
        es.sigma = es.sigma/fac;
        es.pc = fac * es.pc;
        es.C = fac^2 * es.C;
      end

      es.D = sqrt(es.C); // D contains standard deviations now
      // es.D = es.D * N / sum(es.D);

      // set back to normal matrices
      es.BD = diag(es.D);
      es.D = diag(es.D);  // the output assumes matrices
      es.C = diag(es.C);

      if ceil(es.sp.separable) == es.countiter
        es.sp.ccov = es.sp.ccov_def;
      end

    // the original stuff
    elseif es.sp.ccov > 0 & modulo(es.countiter, 1/es.sp.ccov/N/5) < 1
      es.C     = triu(es.C)+triu(es.C,1)';  // enforce symmetry
      [es.B,es.D] = spec(es.C);             // eigen decomposition, B==normalized eigenvectors
      // error management
      if max(diag(es.D)) > 1e14*min(diag(es.D)) | min(diag(es.D)) <= 0
        if 11 < 3
          error('min(diag(es.D)) <= 0');
        elseif 1 < 3
          es.out.stopflags($+1) = 'conditioncov';
        else
          // limit condition of C to 1e14 + 1
          warning('condition number > 1e14 or min(diag(es.D)) = ' + ...
              string(min(diag(es.D))) + ' <= 0');
          tmp = max(diag(es.D))/1e14;
          es.C = es.C + tmp*eye(N,N); es.D = es.D + tmp*eye(N,N);
        end
      end
      es.D = diag(sqrt(diag(es.D))); // D contains standard deviations now
      es.BD = es.B*es.D;
    end

    // end implementation of core algorithm

    // Set stopflags, termination
    if (arfitness(1) <= es.param.stop.fitness)
      es.out.stopflags($+1) = 'fitness';
    end
    if es.counteval >= es.param.stop.maxfunevals
      es.out.stopflags($+1) = 'maxfunevals';
    end
    if es.countiter >= es.param.stop.maxiter
      es.out.stopflags($+1) = 'maxiter';
    end
    if and(es.sigma*(max(abs(es.pc), sqrt(diag(es.C)))) < es.param.stop.tolx)
      es.out.stopflags($+1) = 'tolx';
    end
    if or(es.sigma*sqrt(diag(es.C)) > es.param.stop.tolupx)
      es.out.stopflags($+1) = 'tolupx';
    end
    if max([es.fitnesshist(~isnan(es.fitnesshist)); arfitness']) ...
      - min([es.fitnesshist(~isnan(es.fitnesshist)); arfitness']) <= es.param.stop.tolfun
      es.out.stopflags($+1) = 'tolfun';
    end
    if es.countiter >= length(es.fitnesshist) & ...
      max(es.fitnesshist) - min(es.fitnesshist) <= es.param.stop.tolfunhist // <= makes a difference with zero
      es.out.stopflags($+1) = 'tolfunhist';
    end
    if max(diag(es.D)) / min(diag(es.D)) > 1e7 // captured above
//      for i = 1:N
//          es.C(i,i) = es.C(i,i) + max(diag(es.D))^2 / 1e14;
//      end
//      es.out.stopflags($+1) = 'conditioncov';
    end
    if or(es.xmean == es.xmean + 0.2*es.sigma*sqrt(diag(es.C)))
      es.out.stopflags($+1) = 'noeffectcoord';
    end

    if es.sp.separable <= es.countiter
      if or(es.xmean == es.xmean + 0.1*es.sigma*diag(es.D));
        es.out.stopflags($+1) = 'noeffectaxis';
      end
    else
      if and(es.xmean == es.xmean + 0.1*es.sigma*es.BD(:,1+floor(modulo(es.countiter,N))))
        es.out.stopflags($+1) = 'noeffectaxis';
      end
    end
    if max([es.fitnesshist; arfitness(1:ceil(1+es.sp.lambda/2))']) ...
        - min([es.fitnesshist; arfitness(1:ceil(1+es.sp.lambda/2))']) == 0
      es.out.stopflags($+1) = 'equalfunvals';
    end
    // read stopping message from file
    if ~isempty(es.param.verb.readfromfile)
      [fid, err] = mopen(es.param.verb.readfromfile, 'r');
      if err == 0
        [c s1 s3 idx] = mfscanf(fid, ' %s %s %s %i');
        if c > 1
          // 'stop' sets stopflags to manual
          if s1 == 'stop'
            es.out.stopflags($+1) = 'manual';
          end
        end
        mclose(fid);
      end // err == 0
    end // es.param.verb.readfromfile not empty

    es.stop = ~isempty(es.out.stopflags);

    // Keep overall best solution and some more data
    if es.countiter == 1
      es.out.solutions.bestever.x = cma_genophenotransform(es.out.genopheno, es.xmean);
      es.out.solutions.bestever.f = %inf;  // first hack
      es.out.solutions.bestever.evals = es.counteval;
    end

    fmean = %nan; // a simple hack to have fmean defined without additional feval
    es.out.evals = es.counteval;
    es.out.iterations = es.countiter;
    es.out.solutions.iterations = es.countiter;
    es.out.solutions.mean.x = cma_genophenotransform(es.out.genopheno, es.xmean);
    es.out.solutions.mean.f = fmean;
    es.out.solutions.mean.evals = es.counteval;
    es.out.solutions.recentbest.x = cma_genophenotransform(es.out.genopheno, arx(:, arindex(1)));
    es.out.solutions.recentbest.f = arfitness(1);
    es.out.solutions.recentbest.evals = es.counteval + arindex(1) - es.sp.lambda;
    es.out.solutions.recentworst.x = cma_genophenotransform(es.out.genopheno, arx(:, arindex($)));
    es.out.solutions.recentworst.f = arfitness($);
    es.out.solutions.recentworst.evals = es.counteval + arindex($) - es.sp.lambda;
    if arfitness(1) < es.out.solutions.bestever.f
      es.out.solutions.bestever.x = cma_genophenotransform(es.out.genopheno, arx(:, arindex(1)));
      es.out.solutions.bestever.f = arfitness(1);
      es.out.solutions.bestever.evals = es.counteval + arindex(1) - es.sp.lambda;
    end

    writeOutput(es, arx(:,arindex(1)), arfitness);

    // plot output data
    if es.param.verb.plotmodulo & ...
        (((es.countiter == 3 | ~modulo(1e-12*round(1e12*log10(es.countiter+1)),1))) ...
        | modulo(es.countiter-1,es.param.verb.plotmodulo) == 0 ...
        | ~isempty(es.out.stopflags))
      if es.countiter > 1
        cma_plotintern(es.param.verb.filenameprefix, 324);
      end
      if es.countiter == 2
        if isempty(strindex(get(gcf(), "figure_name"), ...
            es.param.verb.filenameprefix)) //"Plots of CMA-ES (Fig 324)"))
          // addmenu(324, 'Plot Now', 'plot now (refresh)', action=list(2, 'cma_plotintern'));
          set(gcf(), "figure_name", es.param.verb.filenameprefix);
        end
      end
    end

    // print / display

    if es.param.verb.displaymodulo & ...
        (es.countiter <= 3 | ~isempty(es.out.stopflags) | ...
         modulo(es.countiter-1,es.param.verb.displaymodulo) == 0)
      if modulo(es.countiter-1,10*es.param.verb.displaymodulo) == 0
        printf(es.verb.dispannotationline);
      end
      [Cmin, imin] = min(diag(es.C));
      [Cmax, imax] = max(diag(es.C));
      printf('%4d,%6d: %+14.7e +(%.0e)| %.2e |%3d:%.2e,%2d:%.2e\n', ...
          es.countiter, es.counteval, arfitness(1), arfitness($)-arfitness(1), ...
          max(diag(es.D))/min(diag(es.D)), ...
          imin, es.sigma*sqrt(Cmin), imax, es.sigma*sqrt(Cmax));
    end

endfunction // tell

// --------------------------------------------
function writeOutput(es, xrecent, arfitness)
// writes data to four or five output files
//

    // Write output data to file
    if es.param.verb.logmodulo & (es.countiter < 10 ...
                 | (es.countiter < es.param.verb.logmodulo & modulo(log10(es.countiter),1) == 0) ...
                 | modulo(es.countiter-1,es.param.verb.logmodulo) == 0 ...
                 | ~isempty(es.out.stopflags) ...
                 | (~isempty(es.param.verb.logfunvals) & arfitness(1) <= es.param.verb.logfunvals(1)))

      es.param.verb.logfunvals(es.param.verb.logfunvals>=arfitness(1)) = [];

      // Save output data to files
      for name = es.files.additions
        [fid, err] = mopen(es.param.verb.filenameprefix+name+'.xml', 'a');
        if err ~= 0
          warning('could not open '+es.param.verb.filenameprefix+name+'.xml' + ' ');
        else
          if name == 'axlen'
            mfprintf(fid, '%d %d %e %e %e ', es.countiter, es.counteval, es.sigma, ...
                max(diag(es.D)), min(diag(es.D)));
            mfprintf(fid, msprintf('%e ', diag(es.D)) + '\n');
          elseif name == 'disp' // TODO
          elseif name == 'fit'
            mfprintf(fid, '%ld %ld %e %e %25.18e %25.18e %25.18e %25.18e', ...
                es.countiter, es.counteval, es.sigma, max(diag(es.D))/min(diag(es.D)), ...
                es.out.solutions.bestever.f, ...
                arfitness(1), median(arfitness), arfitness($));
// TODO contraints
//            if length(vecfitness(:,1)) > 1
//              mfprintf(fid, ' %25.18e', vecfitness(1:$, arindex(1)));
//            end
            mfprintf(fid, '\n');
          elseif name == 'stddev'
            mfprintf(fid, '%ld %ld %e 0 0 ', es.countiter, es.counteval, es.sigma);
//qqq            mfprintf(fid, (msprintf('%e ', es.sigma*es.out.genopheno.scaling.*sqrt(diag(es.C)))) + '\n');
            mfprintf(fid, (msprintf('%e ', es.sigma*sqrt(diag(es.C)))) + '\n');
          elseif name == 'xmean'
            if isnan(fmean)
              mfprintf(fid, '%ld %ld 0 0 0 ', es.countiter, es.counteval);
            else
              mfprintf(fid, '%ld %ld 0 0 %e ', es.countiter, es.counteval, fmean);
            end
//qqq            mfprintf(fid, msprintf('%e ', cma_genophenotransform(es.out.genopheno, es.xmean)) + '\n');
            mfprintf(fid, msprintf('%e ', es.xmean) + '\n');
          elseif name == 'xrecentbest'
            mfprintf(fid, '%ld %ld %25.18e 0 0 ', es.countiter, es.counteval, arfitness(1));
            mfprintf(fid, msprintf('%e ', ...
                cma_genophenotransform(es.out.genopheno, xrecent)) + '\n');
          end
          mclose(fid);
        end
      end // for
    end // if logmodulo
endfunction

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

function param = check_and_merge_param(inparam, defparam)
// Checks input parameter list inparam w.r.t. valid fields:
// defparam defines all valid fields.
// Merges both parameter lists in that values from inparam
// are taken if available, otherwise from defparam.
// The function is recursively applied into lists and structs.
// Returns param == complemented parameter list

  if isempty(inparam)
    param = defparam;
    return;
  end
  names = getfield(1, inparam);
  names = names(2:$); // remove type/class name
  if typeof(inparam) == 'st' // remove also second entry of struct
    names = names(2:$)
  end
  defnames = getfield(1, defparam);
  defnames = defnames(2:$); // remove type/class name
  if typeof(defparam) == 'st' // remove also second entry of struct
    defnames = defnames(2:$)
  end

  param = defparam;
  // merge into param
  for name = names(1:$)
    // check
    if ~ or(name == defnames(1:$))
      error('unknown parameter fieldname: ' + name);
    // assign
    elseif ~isempty(inparam(name))
      if or(type(defparam(name)) == (15:17))
//pause
        param(name) = check_and_merge_param(inparam(name), defparam(name));
//pause
      else
        param(name) = inparam(name);
      end
    end
  end

endfunction

// --------------------------------------------
function res = cma_tlistgenopheno(typical_x, scales, x0)
// TODO (not urgent as it only matters for the internal
//       representation of the typical solution):
//       redesign this in terms of typical_geno and typical_pheno
//       preferably top-down from input to cmaes.
//       typical_pheno==typical_x or x0 and
//       typical_geno = typical_geno when scales==1
// set up of an affine linear transformation using typical_x (default 0) for
//   translation and scales for scaling.
// typical_x is also used for the fixed variables,
//    otherwise x0 if typical_x is empty.
//
   if isempty(typical_x)
     typical_x = zeros(scales);
     typical_x(scales == 0) = x0(scales == 0);
   end
   skip  = 0;
   if and(typical_x == 0) & and(scales == 1)
     skip = 1;
   end
   res = tlist(['genopheno';'typical_x';'scaling';'skip'], ...
       typical_x, scales, skip);
endfunction

// --------------------------------------------
function y = cma_genophenotransform(gp, x)
// input argument: tlist (an object), x vector to be transformed
// elements in gp :
//         typical_x = Nx1 vector of middle of domains, default == 0
//         scaling = Nx1 vector, default == 1
  if gp.skip
    y = x;
  else
    y = zeros(gp.typical_x); // only to set the size
    y(gp.scaling > 0) = x;
    y = gp.typical_x + y .* gp.scaling;
  end
endfunction

// --------------------------------------------
function x = cma_phenogenotransform(gp, x)
// input argument: tlist (or mlist) according to Gilles Pujol
  if ~gp.skip
    idx = gp.scaling > 0;
    x = (x(idx) - gp.typical_x(idx)) ./ gp.scaling(idx);
  end
endfunction


// --------------------------------------------
function scales = cma_check_initial_sizes(typical_x, x0, scales)
// checks sizes of x0, typical_x, scales (can be a scalar),
// checks that scales >= 0,
// returns Nx1 vector scales
// inputs must be defined but can be empty, if appropriate.

  // --- extract typical_x
  if ~isempty('typical_x') & ~sum(size(typical_x)) == 0
    if typeof(typical_x) == 'string'
      typical_x = evstr(typical_x);
    end
  else // x0 must be defined
    typical_x = zeros(x0);
  end

  // --- check typical_x
  N = size(typical_x, 1);
  if ~and(size(typical_x) == [N 1]) | N == 1
    error('CMAES: typical_x or x0 must be a column vector');
  end

  // --- extract and check x0
  if ~isempty('x0') & typeof(x0) == 'string'
    x0 = evstr(x0);
    if ~and(size(x0) == [N 1])
      error('CMAES: x0 must be a column vector and agree with typical_x');
    end
  end

  // --- extract and check scales
  if isempty(scales)
    scales = ones(N,1);
  elseif typeof(scales) == 'string'
    scales = evstr(scales);
  end
  if 11 < 3 & and(size(scales) == [1 1]) // sigma is used as initial step-size
    scales = scales * ones(N,1);
  end
  if ~and(size(scales) == [N 1])
    error('CMAES: input parameter scales must be a column vector and agree with x0 and typical_x');
  end
  if or(scales < 0)
    error('CMAES: input parameter scales must >= zero');
  end

endfunction

// --------------------------------------------
function x0 = cma_xstart(typical_x, x0, scales)
// computes xstart from the input parameters
// all in phenotypic space
  if ~isempty(x0)
    if typeof(x0) == 'string'
      x0 = evstr(x0);
    end
  else // generate x0 from typical_x
    if typeof(typical_x) == 'string'
      typical_x = evstr(typical_x);
    end
    // Rather no uniform distribution here, because of invariance
    x0 = typical_x + scales .* rand(scales,'normal');
  end
endfunction

/////////////////////////////////////////////////////////////////////////////

//_______________________________________________________
//_______________________________________________________
//
function cma_annotate(xfinal, yfinal, ann)
// not in use but potentially useful
// INPUT
//      xfinal = xdata(idx($));   // single last abscissa value or two values (first and last)
//      yfinal = ydata(idx($),:); // y-values to be annotated
//      ann = string(condef(1,:)'); // row-wise strings to be printed

if length(xfinal) == 2
  xfirst = xfinal(1);
  xfinal = xfinal(2);
end
      // annotation "function", should become a stand alone function soon!?
      yrange = get(gca(), "y_ticks");
      yrange = yrange.locations([1 $]);
      yrange(1) = yrange(1) + 1e-6*(yrange(2) - yrange(1)); // prevents widening of range
      yrange(2) = yrange(2) - 1e-6*(yrange(2) - yrange(1));
      xrange = get(gca(), "x_ticks");
      xrange = xrange.locations([1 $]);
      xlast = max([1.07 * xfinal xrange(2)]);
      [ignore, idx] = gsort(yfinal);
      idx2 = [];
      Ngeno = length(yfinal);
      idx2(idx) = (Ngeno-1:-1:0)';
      plot(xfinal*ones(1,2), yrange, 'k', 'linewidth', 1);
      plot([xfinal; xlast]', ...
          [yfinal; ...
              yrange(1) + (idx2')*(yrange(2)-yrange(1))/(Ngeno-1)]);

      set(gca(), "clip_state", "off"); // numbers visible everywhere
      str = '';
      for i = 1:Ngeno
        if i > 9
          str = '';
        end
        xstring(xlast, yrange(1) + ...
            (yrange(2)-yrange(1))*(idx2(i)/(Ngeno-1) - 0.3/max(10,Ngeno)), ...
            [str ann(i)]);
      end

endfunction

//_______________________________________________________
//_______________________________________________________
//
function cma_plot(fignb, name_prefix, name_extension, object_variables_name, plotrange)
//
// plots data from CMA-ES, so far for Scilab and Java.
// defaults: name_prefix='outcmaes', name_extension='.xml'
//           fignb=326

defaults.name_prefix = 'outcmaes';
defaults.name_extension = '.xml';
defaults.object_variables = 'xmean'; // or 'xrecentbest'

if ~isdef('name_prefix', 'local')
  name_prefix = defaults.name_prefix;
  if isdef('filenameprefix') // inherit from upper environment
    name_prefix = filenameprefix;
  end
elseif isempty(name_prefix)
  name_prefix = defaults.name_prefix;
end
if ~isdef('name_extension', 'local')
  name_extension = defaults.name_extension;
elseif isempty(name_extension)
  name_extension = defaults.name_extension;
end
if ~isdef('fignb', 'local')
  fignb = 326;
end

if ~isdef('object_variables_name', 'local')
  object_variables_name = defaults.object_variables;
end

if ~isdef('plotrange', 'local')
  plotrange = [];
elseif length(plotrange) ~= 2
  error(['plotrange must have length 2, not ' string(length(plotrange))]);
end

  names = ['axlen', 'fit', 'stddev', object_variables_name];
//  names = ['axlen', 'fit', 'stddev', 'xmean', 'xrecentbest'];
  data = tlist(['data', names]);
// columns are generation feval something y-values

  // read data from files into tlist data
  for name = names
    [fid, err] = mopen([name_prefix + name + name_extension]);
    if err ~= 0
      warning('File ' + [name_prefix + name + name_extension] + ...
          ' could not be opened');
      data(name) = [];
    else
      if 1 < 3  // reading quick and dirty
        mclose(fid);
        data(name) = fscanfMat([name_prefix + name + name_extension]);
      else // this can become a more fail save version of reading the data
        // TODO: check whether this also works if several DATA entries are present
        s = mgetl(fid); // read complete file
        mclose(fid);
        // remove headings up to < DATA ... >
        idx = grep(s, 'DATA'); // TODO: make fail save by searching for < DATA * >
        idx2 = grep(s, '>'); // lines where > character is found
        idx2 = idx2(idx2>=idx(1));
        if ~isempty(idx2)
          s = s(idx2(1)+1:$); // TODO keep eventually also the header stuff
        end

        // remove trailing non numerical stuff
        idx = grep(s, '<');
        if ~isempty(idx)
          s = s(1:idx(1)-1);
        end
        data(name) = evstr(s);
        //       data(name) = mfscanf(-1, fid);
      end
      if ~isempty(plotrange)
        idx = data(name)(:,2) >= plotrange(1) & data(name)(:,2) <= plotrange(2);
        data(name) = data(name)(idx,:);
      end
    end
  end

// TODO the legend interferes with large negativ fitnesses

// comments:
// data.*(:,1) == iteration number
// data.*(:,2) == function evaluations
// data.*(:,6:$) == the vector data like xmean
// data.fit(:,6) == recent best fitness

  flg_raw_stds = %F; // remove sigma from stds
  if flg_raw_stds & ~isempty(data.stddev)
    for i = 6:size(data.stddev,2)
      data.stddev(:,i) = data.stddev(:,i) ./ data.stddev(:,3);
    end
  end

  scf(fignb); // set current figure
  drawlater;
  clf; // is needed only because of the annotation of variables


  //////////////////
  subplot(2,2,1);
  if ~isempty(data.fit)
// TODO remove this (it is for testing purpose only)
// data.fit(:,6) = data.fit(:,6) - 2e5;
    xgrid;
    [fmin, idxmin] = min(data.fit(:,6));
    // idxmin = find(data.fit(:,6) == fmin);

    xtitle("Function Value (fval, fval minus f_min), Sigma (g), Axis Ratio (r)", ...
        'f_recent='+sprintf('%25.18e', data.fit($, 6)), ...
        'log10(abs(value))');

    // plot legend for best function value
    plot(data.fit(:,2)(idxmin), log10(abs(fmin) + 1e-49), 'r*');
    legend('f_best=' + sprintf('%25.18e', fmin));
    data.fit(idxmin,6) = %nan;

    // additional fitness data, eventually objectives and constraints values
    if size(data.fit,2) > 8
      d = log10(abs(data.fit(:,9:$) + 1e-29));
      d(data.fit(:,9:$)==0) = %nan;  // cannot be the case
      plot(data.fit(:,2), d, 'm-');  // pos. and neg. entries
      d(data.fit(:,9:$) > 0) = %nan;
      plot(data.fit(:,2), d, 'red'); // neg. entries
    end

    // plot abs function value in blue, red, and black
    idxpos = data.fit(:,6) > 0;
    idxneg = data.fit(:,6) < 0; // removes nan and zeros

    if or(idxpos) // if any entry in idxpos is true
      plot(data.fit(idxpos,2), ...
          log10(0e-49 + abs(data.fit(idxpos,6))), 'b.');
    end

    if or(idxneg)
      plot(data.fit(idxneg,2), ...
          log10(0e-49 + abs(data.fit(idxneg,6))), 'm.');
    end

    // median and worst fitness
    if size(data.fit,2) > 6
      plot(data.fit(:,2), log10(abs(data.fit(:,7:8)) + 1e-99), 'k-');
    end

    // plot function value differences disregarding all fmins
    // careful: unfortunately log10 cannot handle %nan (in some versions)
    idx = ~isnan(data.fit(:,6));
    plot(data.fit(idx,2), log10(data.fit(idx,6) - fmin + 1e-49), 'c-');

    // plot marker(s) for best function value _over_ the graph now, see legend above
    plot(data.fit(:,2)(idxmin), log10(abs(fmin) + 1e-49), 'r*');
    plot(data.fit(:,2)(idxmin), log10(cma_minnan(data.fit(:,6)-fmin) + 1e-49)*ones(idxmin), 'r*');

  end // data.fit is not empty

  // plot sigma
  if ~isempty(data.stddev)
    plot(data.stddev(:,2), log10(data.stddev(:,3)), 'g');
  else
    plot(data.fit(:,2), log10(data.fit(:,3)), 'g');
  end

  // plot axis ratio
  plot(data.fit(:,2), log10(data.fit(:,4)), 'r');

  //////////////////
  subplot(2,2,2);
  name = object_variables_name; // 'xmean' or 'xrecentbest'
  if ~isempty(data(name))
    xgrid;
    xtitle('Object Variables (' + name + ', ' + ...
        string(size(data(name),2)-5) + '-D)');

    plot(data(name)(:,2), data(name)(:,6:$)); // this should optional be the genotyp, annotation is tricky

    // annotations of variables with numbers
    Ngeno = size(data(name), 2) - 5;
    if Ngeno < 100

      yrange = get(gca(), "y_ticks");
      yrange = yrange.locations([1 $]);
      yrange(1) = yrange(1) + 1e-6*(yrange(2) - yrange(1)); // prevents widening of range
      yrange(2) = yrange(2) - 1e-6*(yrange(2) - yrange(1));
      xrange = get(gca(), "x_ticks");
      xrange = xrange.locations([1 $]);
      xlast = max([1.07 * data(name)($,2) xrange(2)]);
      [sorted_x, idx] = gsort(data(name)($,6:$));
      idx2 = [];
      idx2(idx) = (Ngeno-1:-1:0)';
      plot(data(name)($,2)*ones(1,2), yrange, 'k', 'linewidth', 1);
      plot([data(name)($,2); xlast]', ...
          [data(name)($, 6:$) ; ...
              yrange(1) + (idx2')*(yrange(2)-yrange(1))/(Ngeno-1)]);

      set(gca(), "clip_state", "off"); // numbers visible everywhere
      str = ' ';
      for i = 1:Ngeno
        if i > 9
          str = '';
        end
        xstring(xlast, yrange(1) + ...
            (yrange(2)-yrange(1))*(idx2(i)/(Ngeno-1) - 0.3/max(10,Ngeno)), ...
            [str 'x(' string((i)) ')=' string(data(name)($,i+5))]);
      end
    end // annotation

  end // isempty data. xmean

  //////////////////
  subplot(2,2,3);
  if ~isempty(data.axlen)
    xgrid;
    xtitle("Principle Axis Lengths", 'function evaluations', 'log10(value)');
    plot(data.axlen(:,2), log10(data.axlen(:,6:$)+1e-49));
  end

  //////////////////
  subplot(2,2,4);
  if ~isempty(data.stddev)
    xgrid;
    xtitle("Standard Deviations", 'function evaluations', 'log10(value)');

    plot(data.stddev(:,2), log10(data.stddev(:,6:$)));

    // annotations of variables with numbers
    Ngeno = size(data.stddev, 2) - 5;
    if Ngeno < 100
      yrange = get(gca(), "y_ticks");
      yrange = yrange.locations([1 $]);
      yrange(1) = yrange(1) + 1e-6*(yrange(2) - yrange(1));
      yrange(2) = yrange(2) - 1e-6*(yrange(2) - yrange(1));
      xrange = get(gca(), "x_ticks");
      xrange = xrange.locations([1 $]);
      xlast = max([1.07 * data.stddev($,2) xrange(2)]);
      [sorted_x, idx] = gsort(data.stddev($,6:$));
      idx2 = [];
      idx2(idx) = (Ngeno-1:-1:0);
      plot(data.stddev($,2)*ones(1,2), yrange, 'k', 'linewidth', 1);
      plot([data.stddev($,2); xlast]', ...
          [log10(data.stddev($,6:$)); ...
              yrange(1) + (yrange(2)-yrange(1))*(idx2/(Ngeno-1))]);

      set(gca(), "clip_state", "off");
      str = ' ';
      idxvars = 1:Ngeno; // find(out.genopheno.scaling > 0);
      for i = 1:Ngeno
        if i > 9
          str = '';
        end
        xstring(xlast, yrange(1) + ...
            (yrange(2)-yrange(1)) * (idx2(i)/(Ngeno-1) - 0.3/max(Ngeno,10)), ...
            [str string(idxvars(i))]);
      end
    end
  end // ~isempty(data.stddev)

  drawnow;
  return data;

endfunction

//_______________________________________________________
//_______________________________________________________
//
function res = cma_minnan(in)

  res = min(in(~isnan(in)));

endfunction

//_______________________________________________________
//_______________________________________________________
//
function cma_plotintern(name_prefix, fignb, name_extension, object_variables_name)
//
// plots data from CMA-ES, so far for Scilab and Java.
// defaults: name_prefix='outcmaes', name_extension='.xml'
//           fignb=326

// Can be called by Plot Now button, therefor the strange
// arrangement of input elements.

  defaults.name_prefix = 'outcmaes';
  defaults.name_extension = '.xml';
  defaults.object_variables = 'xmean'; // or 'xrecentbest'

  if ~isdef('name_prefix', 'local')
    name_prefix = defaults.name_prefix;
    if isdef('filenameprefix') // inherit from upper environment
      name_prefix = filenameprefix;
    end
  elseif isempty(name_prefix)
    name_prefix = defaults.name_prefix;
  elseif name_prefix == 1 // called from "Plot Now" button
    name_prefix = get(gcf(), "figure_name");
  end
  if ~isdef('name_extension', 'local')
    name_extension = defaults.name_extension;
  elseif isempty(name_extension)
    name_extension = defaults.name_extension;
  end
  if ~isdef('fignb', 'local')
    fignb = 326;
  end

  if ~isdef('object_variables_name', 'local')
    object_variables_name = defaults.object_variables;
  end

  cma_plot(fignb, name_prefix, name_extension, object_variables_name);

endfunction

// CHANGES LOG
// November 2010, v0.988: replace sort with gsort, update constant cc
// May 2009, v0.997: bug fix for idxmin in plot in case of several best f-values
// June 2008, v0.993: bug fix for runcma.sce: cma_ask does not return a
//      list for lambda==1 anymore.
// June 2008, v0.991: input parameter check in cma_optim removed (is done
//      in cma_new anyway). warning on nan in cma_optim removed. option
//      readfromfile can be empty now.
// May 2008: The interface of ask and tell changed: a list of solutions rather than
//   a row array is passed now. This makes the interface safer. The possiblibity
//   to evaluate the population as a matrix is impeded.


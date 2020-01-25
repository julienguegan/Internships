TOOLBOX
cma-toolbox  for stochastic non-linear non-convex minimization with the
  covariance matrix adapation evolution strategy (CMA-ES)

installation with exec(path_to_cma_toolbox + 'builder.sce') // nothing really weird happens here
loading with      exec(path_to_cma_toolbox + 'loader.sce')  // only load, assumes the library 
                                                            // has been build

FILES
the toolbox contains the following source files 
  cma.sci         all core functions, CAVE: visible only after cma() was
                                            called which happens in loader/builder
  cma_optim.sci   function cma_optim(...) minimizer with functional interface
  cma_optim_restarted.sci   function cma_optim_restarted(...) similar to cma_optim
  runcma.sce      small demo, a good place to go from here 
  fitfuns.sci     a volatile collection of benchmark functions, not loaded with the library  

FUNCTIONS 
external functions, interfaces
  // functional interface
  function [xopt, f, out, param] = cma_optim(costf, x0, sigma0, param) 
  function [xopt, f, out, param] = cma_optim_restarted(costf, x0, sigma0, restarts, param) 
  // object oriented ask-and-tell interface
  function [es, p] = cma_new(inparam)  // constructor
  function data    = cma_plot(fignb, name_prefix, name_extension, object_variables_name, plotrange)
  function X       = cma_ask(es, lam)  // delivers lam new points
  function es      = cma_tell(es, X, arfitness) // update, concludes iteration
  function [y, x]  = cma_best(es)      // best-ever solution, see also es.out.solutions
  function flag    = cma_stop(es)      // termination criterion met?
  function X       = %cma_ask(es)      // delivers const es.param.opt.lambda new points
  function es      = %cma_tell(es, X, arfitness) // update, concludes iteration
  function [y, x]  = %cma_best(es)     // best-ever solution, see also es.out.solutions

eventually useful functions
  function param   =  check_and_merge_param(inparam, defparam) // verifies and merges 
                                                               // two parameter lists
  function y       =  cma_genophenotransform(gp, x)   // gp=es.out.genopheno
  function x       =  cma_phenogenotransform(gp, x)   // transform into genotypic representation
  function res     =  cma_minnan(vec)                 // utility, compute min disregarding %nan

internal functions 
  function p       =  cma_getdefaults()
  function gp      =  cma_tlistgenopheno(typical_x, scales, x0)
  function scales  =  cma_check_initial_sizes(typical_x, x0, scales)
  function x0      =  cma_xstart(typical_x, x0, scales)  
  function            cma_annotate(xfinal, yfinal, ann) // for later use
  function            cma_plotintern(name_prefix, fignb, name_extension, object_variables_name)



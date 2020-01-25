function [xopt, f, out, param] = cma_optim_restarted(costf, x0, sigma0, restarts, param)
// implements restarts with increasing population size, see cma_optim()  
  if ~isdef('costf', 'local') | ~isdef('x0', 'local') ...
      | ~isdef('sigma0', 'local')
    if ~isdef('costf', 'local') 
      [xopt, f, out, param] = cma_optim();
    else
      [xopt, f, out, param] = cma_optim([]);
    end
    return
  end
  if ~isdef('param', 'local')
    param = {}; // [] works as well
  end
  fullparam = check_and_merge_param(param, cma_new([]));
  p = fullparam; 
  done = %f;
  for irun = 1:1+restarts
    [xopt, f, out, param] = cma_optim(costf, x0, sigma0, p); 
    if irun == 1
      f_prev = f; 
      xopt_prev = xopt;
    elseif f_prev < f
      f = f_prev;
      xopt = xopt_prev;
    end
    for i = 1:size(out.stopflags)
      if out.stopflags(i) == 'fitness' | out.stopflags(i) == 'maxfunevals' 
	done = %t;
      end
    end
    if done | irun > restarts
      break; 
    end 
    p.opt.lambda = 2 * param.opt.lambda; 
    p.verb.append = out.evals; 
  end // for irun
  if irun > 1
    disp('  ' + string(irun) + ' runs')
  end
endfunction


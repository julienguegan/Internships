function [xopt, f, out, param] = cma_optim(costf, x0, sigma0, param)
// functional interface to object oriented optimizer cma
// param are optional parameters

  if ~isdef('costf', 'local') | ~isdef('x0', 'local') ...
      | ~isdef('sigma0', 'local')
    o = cma_new([]); 
    xopt = [];
    names = getfield(1,o);
    for name = names(3:$)
      if name ~= 'x0' & name ~= 'sigma0'  
	xopt(name) = o(name);
      end
    end
    if ~isdef('costf', 'local') | ~isempty(costf)
      disp('cma_optim has three mandatory parameters: ' ...
	  + 'costf, x0, sigma0. Default optional param returned');
    end
    return;
  end
  if ~isdef('param', 'local')
    if isdef('param')
      warning('inherited value for param disregarded');
    end
    param = [];
  end
  param.x0 = x0;
  param.sigma0 = sigma0;
  [es, param] = cma_new(param); 
  while ~ cma_stop(es) // es.out.stopflags is empty
    X = cma_ask(es);   // returns a matrix of lambda column vectors
    y = [];            // just in case
    for i = 1:length(X)
      y(i) = costf(X(i));
      // treat %nan
      while isnan(y(i)) 
	X(i) = cma_ask(es, 1); 
	y(i) = costf(X(i)); 
        // warning('costf(x)==%nan, ' + string(i) + '. solution resampled');
      end
    end
    es = cma_tell(es, X, y');
  end
  [f, xopt] = cma_best(es);
  fmean = costf(es.xmean);
  if fmean <= f
    f = fmean;
    xopt = es.xmean;
  end
  out = es.out;
  if param.verb.displaymodulo
    for str = out.stopflags
      printf('terminated on ''' + str + '''\n');
    end
  end

endfunction 


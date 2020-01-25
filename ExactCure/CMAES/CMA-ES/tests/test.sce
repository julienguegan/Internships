// exec(get_absolute_file_path('test.sce') + 'builder.sce');
// clear; exec(get_absolute_file_path('test.sce') + 'loader.sce');

if 1 < 3 // case runcma.sce
  exec(get_absolute_file_path('test.sce') + '/runcma.sce');
  cma_plot(320);
end

if 1 < 3 // case cma_optim_restarted
  p = cma_optim([]);
  p.verb.logmodulo = 10;
  p.verb.plotmodulo = 1e9; // plot rarely
  p.stop.fitness = 1e-4;
  p.opt.separable = '10*N*10/lambda'; // adaptation time of tablet

  if ~isdef('frastrigin')
    getf(get_absolute_file_path('test.sce') + '/fitfuns.sci');
  end
  [x f out param] = cma_optim_restarted(frastrigin10, '10*ones(5,1)-2', 3, 5, p);
  // cma_plot(321);
end


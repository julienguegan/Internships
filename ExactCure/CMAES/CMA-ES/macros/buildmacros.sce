// old code:
// printf('Building macros in %s\n', get_absolute_file_path('buildmacros.sce'));
// genlib('cmalib', get_absolute_file_path('buildmacros.sce'), %t);
// cma(); // only the first reference makes all functions from cma.sci available (stupid Scilab)
// // cma_optim([]); // cma_optim_restarted

macros_path = get_absolute_file_path("buildmacros.sce");
  tbx_build_macros(TOOLBOX_NAME, macros_path);
clear macros_path;
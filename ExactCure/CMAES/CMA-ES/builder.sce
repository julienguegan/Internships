// Copyright (C) 2008 - INRIA
// Copyright (C) 2009-2011 - DIGITEO

// This file is released under the 3-clause BSD license. See COPYING-BSD.

mode(-1);
lines(0);

function main_builder()

  TOOLBOX_NAME  = "cma_es";
  TOOLBOX_TITLE = "CMA-ES Optimization Toolbox";
  toolbox_dir   = get_absolute_file_path("builder.sce");

// Check Scilab's version
// =============================================================================

  try
    v = getversion("scilab");
  catch
    error(gettext("Scilab 5.5.1 or more is required."));
  end

  if v(2) < 3 then
    // new API in scilab 5.3
    error(gettext('Scilab 5.3 or more is required.'));
  end

// Check modules_manager module availability
// =============================================================================

if ~isdef('tbx_build_loader') then
  error(msprintf(gettext('%s module not installed."), 'modules_manager'));
end

// Action
// =============================================================================

tbx_builder_macros(toolbox_dir); // not tested
tbx_builder_help(toolbox_dir);
tbx_build_loader(TOOLBOX_NAME, toolbox_dir); // generate loader.sce
tbx_build_cleaner(TOOLBOX_NAME, toolbox_dir); // generate cleaner.sce script

// Old code - it works, part way:
// pathB = get_absolute_file_path('builder.sce');
// if isdir(pathB+'macros') then
//   exec(pathB+'macros/buildmacros.sce');
// end
// if isdir(pathB+'help') then
//   exec(pathB+'help/builder_help.sce');
// end
// clear pathB;

endfunction
// =============================================================================
main_builder();
clear main_builder; // remove main_builder on stack
// =============================================================================

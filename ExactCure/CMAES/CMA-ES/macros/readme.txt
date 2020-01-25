In this directory cma.sci has been modified to prevent complaint from Scilab when compiling to binary file

Further more loadmacros.sce has been removed. It should now build when building the toolbox ... but original
loadmacros.sce file contains this comment:

    cma(); // only the first reference makes all functions from cma.sci available (stupid Scilab)

Which is maybe still the fact since cma.sci contains a cma() function, which only purpose is to:

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

... and cma_getdefaults() is the next function defined in cma.sci. It's
a bit longer, but it defines "param" and if these are not called, then
the library does not work properly.

Question?

How should this be "initialized" when started ... in cma-es.start script, maybe?

cma-es.start calls several atoms module functions, like lib(pathmacros), but there is
no "sign" that this will by any standard call cma() (maybe toolbox_name cma-es ?? how
would I know??).

Calling cma() in the start script probably shouldn't be done inside a function call,
because then the globally defined parameters will disappear with the scope of the
function.

... Or maybe it doesn't work (in this case maybe best to include loadmacros.sce directly, like before)
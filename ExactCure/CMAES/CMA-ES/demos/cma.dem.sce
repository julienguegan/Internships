// cma.dem.sce
// THIS SCRIPT DEMONSTRATES HOW TO RUN CMA-ES, USING THE ROSENBROCK FUNCTION.

// ORGANIZE THIS AS A FUNCTION + A FUNCTION CALL, SO THAT WE DO NOT CLUTTER
// OUR VARIABLES BROWSER WITH GLOBAL PARAMETERS. ALL VARIABLES STAY LOCAL.
function cma_demo()

// INITIALIZE (REMOVE IN ATOMS PACKAGE - ASSUMING ATOMS ALREADY INSTALLED)
scriptpath = get_absolute_file_path("cma.dem.sce");
demo_tlbx  = getshortpathname(scriptpath);
// cma.dem.sce must be located in the /demos subdirectory
root_tlbx_dir = strncpy( test_tlbx, length(test_tlbx)-length("\demos\") );

// fitfuns.sci must be located in the /tests subdirectory
fitfunctions = root_tlbx_dir + filesep() + "tests" + filesep() + "fitfuns.sci";
exec(fitfunctions,-1); // load library of test functions

// INTRODUCTION
printf("\n"); // Print some empty lines - instead of clearing console with clc;
printf("\n");
printf("\n");
printf(" CMA-ES Optimization Toolbox - demo\n");
printf(" ----------------------------------\n");
printf("\n");
printf(" The ATOMS package for CMA-ES comes with several test\n");
printf(" functions. You can use these functions to test any\n");
printf(" optimizer. Here we will choose the Rosenbrock function,\n");
printf(" frosen, which is a non-convex, non-separable banana shaped\n");
printf(" test function introduced by Howard H. Rosenbrock, often\n");
printf(" used in performance testing of optimization algorithms.");
response = input(" Press enter to continue:");
printf("\n");

// INITIALIZE
printf(" First initialize our input parameters in a column vector,\n");
printf(" starting in this case in origo (i.e. 0,0) in 10 dimensions:\n");
printf("\n");
x0 = zeros(10,1);
printf(" x0:\n");
disp(x0);
printf("\n");
response = input(" Press enter to continue:");
printf("\n");

// PREPARE CMA-ES EXECUTION
printf(" Then setup CMA-ES parameters - this is a minimum setup.\n");
printf(" Define a stochastic search range with the gaussian\n");
printf(" standard deviation:\n");
printf("\n");
sigma0 = 0.1;
printf(" sigma0:\n");
disp(sigma0);
printf("\n");

// EXECUTE
printf(" Now execute the optimization. You may have to move\n");
printf(" the graph window to see what goes on in the console.\n");
response = input(" Press enter to continue:");
printf("\n");
printf(" xopt = cma_optim(frosen, x0, sigma0)\n");
printf("\n");
// [xopt, f, out, param] = cma_optim(frosen, x0, sigma0);
xopt = cma_optim(frosen, x0, sigma0);
printf("\n");

// SHOW RESULT
printf(" The xopt vector contains the optimum (at 1,1 in 10-D):\n")
printf("\n");
printf(" xopt:\n");
disp(xopt);
printf("\n");
printf(" Done!\n");
endfunction

cma_demo();
clear cma_demo();

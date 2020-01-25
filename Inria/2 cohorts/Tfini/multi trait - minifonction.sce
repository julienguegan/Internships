
M = 0
J1 = 0
J2 = 0
n11 = 0
n21 = 0

M = ode("stiff",[M01;M02;800;800],0,t,tol, syst_dyn) 
J1 = M(3,length(t))
J2 = M(4,length(t))
n11 = (n1*J1*n)/(n1*J1+n2*J2) 
n21 = (n2*J2*n)/(n1*J1+n2*J2)

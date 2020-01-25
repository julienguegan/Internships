// CMA-ES: Evolution Strategy with Covariance Matrix Adaptation for
// nonlinear function minimization. To be used under the terms of the
// GNU General Public License (http://www.gnu.org/copyleft/gpl.html).
//
//
// Author: Nikolaus Hansen, 2003.
// e-mail: hansen[at]bionik.tu-berlin.de
// URL: http://www.bionik.tu-berlin.de/user/niko
// References: See end of file. Last change: October, 27, 2004

// Modification for scilab :
// Author: Yann COLLETTE, 2006.
// e-mail: yann[dot]colletet[at]renault[dot]com
// URL: http://ycollette.free.fr

// Debugged: indexing with floor($/2) is not working in Scilab (NH)

////////////////////
// Test functions //
////////////////////

function f=fvardim(x)
// Buckley 1989 page 98 according to Powell 2004
// Powells improved DFO is about 2 times faster than CMA-ES
// should work parallelized
  N = size(x,1);
  tmp = (1:N) * (x-1); // is a vector if x is a population
  f = sum((x-1).^2, 1) + tmp.^2 + tmp.^4;
endfunction

// infinite condition number function
function f=olivier(x)
  N = size(x,1);
  r = sqrt(sum(x.^2,1));
  //  ang = 90*(1 - abs(x(1))/sqrt(sum(x.^2,1)));
  // cos angle(x,e_1) = <x,e_1> / ||x|| / ||e_1||
  ang = acos(abs(x(1)/r)); // x(1) == <x,e_1>
  f = sum(x.^2)^0.5 + ang^2;
endfunction

// a parabolic ridge in each diagonal (2^N ridges)
// proposed by Olivier, realized by Niko
function f=oli2(x)
  N = size(x,1);
  r = sqrt(sum(x.^2,1));
  v = sign(x) / sqrt(N); // closest diagonal
  vdifference = x - (v'*x) * v; // difference vector to diagonal
  dist = sqrt(sum(vdifference.^2,1)); // distance to diagonal
  f = 1*r^2 + 1e4*(dist)^2 / 1;
  // f = sum(x.^40,1);
  // f = sum(x.^2,1) + (N - (sum(abs(x))/sqrt(sum(x.^2)))^2);
endfunction

function f=fsphere(x)
   // parallelized
   f=sum(x.^2, 1);
   if 1 < 3 & rand(1,1,'uniform') < 0.001;
//     f(1) = %nan;
   end
endfunction

function f=fspherecon(x)
   // parallelized
   N = size(x,1); if (N < 2) then error('dimension must be greater than one'); end
   f=1e6.^((0:N-1)/(N-1)) * x.^2;
   // f=sum(x(1:$,:).^2, 1);
   // f = [f ; x(1,:)+1];
   f = [f ; (x([2:2:$],:) - 20+1e-5)];
//   f = [f ; ([-x(1,:) x([3:2:$],:)] + 20+1e-5)];
endfunction

function f=fspherecorr(x)
   f=1e4*(x(1) - x(2))^2 + 1e0*sum(x(3:$).^2);
endfunction

// sphere with flat epsilon ball in the middle
function f=fsphereconst(x)
   f=sum(x.^2);
   if f < 0.1
      f = 0.01;
   end
endfunction

function f=fspherenoise(x)
  N = size(x,1);
  f=sum(x.^2) * (1 + (1/N/2) * rand(1,1,'normal')/rand(1,1,'normal') + 0/N*rand(1,1,'normal')) + 0. * (rand(1,1,'normal'));
endfunction

function f=fsqrtnoise(x)
  f=sum(x.^2)^0.25 + 0.01 * rand(1,1,'normal');
endfunction

function f=fsphereminusone(x)
  f=sum((x-1).^2);
endfunction

function f=fsphereodd(x)
  f = sum(x(1:2:$,:).^2);
endfunction

function f=fschwefel(x)
f = 0;
for i = 1:size(x,1),
  f = f+sum(x(1:i))^2;
end
endfunction

function f=fschwefel221(x)
  f = max(abs(x));
endfunction

function f=fschwefel222(x)
  f = sum(abs(x)) + prod(abs(x));
endfunction

function f=fschwefelmult(x)
  f = sum(x .* sin(sqrt(abs(x))));
  f = f + 1e-3*(sum(x(x<-500).^2) + sum(x(x>500).^2));
endfunction

function f=fgriewank(x)
// in -600..600
  f = 1 + sum(x.^2)/4000 + prod(cos(x./sqrt(1:size(x,2))'));
endfunction

function f=fcigar(x)
f = x(1)^2 + 1e6*sum(x(2:$).^2);
endfunction

function f=fmultcigar(x)
N = size(x,1);
M = min(3,N);
f = sum(x(1:M).^2) + 1e6*sum(x(M+1:$).^2);
endfunction

function f=fcigtab(x)
f = x(1)^2 + 1e8*x($)^2 + 1e4*sum(x(2:($-1)).^2);
endfunction

function f=ftabletrotselect(x, c) // c=countiter=generation
// select random subspace generational-wise
// works surprisingly well
   N=size(x,1);
   x=coordinatesystem(N)*x;
   s = rand('seed');
   rand('seed', c);
   idx = rand(x)<0.7;
   rand('seed', s);
   f = 1e6*idx(1,:).*x(1,:)^2 + sum(x(idx).^2);
endfunction

function f=ftablet(x)
f = 1e6*x(1)^2 + sum(x(2:$).^2);
endfunction

function f=felli(x)
N = size(x,1); if (N < 2) then error('dimension must be greater than one'); end
f=1e6.^((0:N-1)/(N-1)) * x.^2;
endfunction

function [f, g] = fconelli(x)
  N = size(x,1); if (N < 2) then error('dimension must be greater than one'); end
  f=1e6.^((0:N-1)/(N-1)) * x.^2;

  con.M = min(3,N); // number of constraints
  inc = floor((N-1)/(con.M-1));

  con.idx = 1:inc:1+(con.M-1)*inc;
  con.arweight = 1; // edit this
  con.arvalue = 0.1;
  g = con.arweight .* (x(con.idx) + con.arvalue);

endfunction

function f=felliodd(x) // only odd variables are used
N = size(x,1); if (N < 2) then error('dimension must be greater than one'); end
f=1e6.^((0:2:N-1)/(N-1)) * x(1:2:$,:).^2;
endfunction

function f=felli100(x)
N = size(x,1); if (N < 2) then error('dimension must be greater than one'); end
f=1e4.^((0:N-1)/(N-1)) * x.^2;
endfunction

function f=fplane(x)
f=x(1);
endfunction

function f=ftwoaxes(x)
f = sum(x(1:$/2).^2) + 1e6*sum(x(1+$/2:$).^2);
endfunction

function f=ftwopow(x)
f = sum(abs(x(1:$/2)).^2) + 1e4 * sum(abs(x(1+$/2:$)).^4);
endfunction

function f=fparabR(x)
f = -x(1) + 100*sum(x(2:$).^2)^(1.5/2);
endfunction

function f=fsharpR(x)
f = -x(1) + 100*norm(x(2:$));
endfunction

function f=frosen(x)
if (size(x,1) < 2) then error('dimension must be greater than one'); end
f = 1e2*sum((x(1:$-1).^2 - x(2:$)).^2) + sum((x(1:$-1)-1).^2);
// f = f + 0.01*rand(1,1,'normal'); // /rand(1,1,'normal');
if rand(1) < -0.001
  f = %nan;
end
endfunction

// lesser known rosenbrock, also less difficult
function f=flkrosen(x)
if (size(x,1) < 2) then error('dimension must be greater than one'); end
f = 1e2*sum((x(1:2:$).^2 - x(2:2:$)).^2) + sum((x(1:$-1)-1).^2);
// f = f + 0.01*rand(1,1,'normal'); // /rand(1,1,'normal');
endfunction

function f=fdiffpow(x)
// parallelized
N = size(x,1);
popsize = size(x,2);
if (N < 2) then error('dimension must be greater than one'); end
f=sum(abs(x).^((2+10*(0:N-1)'/(N-1)) * ones(1,popsize)), 1);
endfunction

function f=fdiffpowi(x) //
// parallelized
N = size(x,1);
popsize = size(x,2);
if (N < 2) then error('dimension must be greater than one'); end
f=sum(abs(x).^(((2:N+1)') * ones(1,popsize)), 1);
endfunction

function f=fdiffpow0(x) // sphere function
N = size(x,1);
if (N < 2) then error('dimension must be greater than one'); end
f=sum(abs(x).^(2+0*(0:N-1)'/(N-1)));
endfunction

function f=fdiffpow1(x)
N = size(x,1);
if (N < 2) then error('dimension must be greater than one'); end
f=sum(abs(x).^(2+1*(0:N-1)'/(N-1)));
endfunction

function f=fdiffpow2(x)
N = size(x,1);
if (N < 2) then error('dimension must be greater than one'); end
f=sum(abs(x).^(2+2*(0:N-1)'/(N-1)));
endfunction

function f=fdiffpow3(x)
N = size(x,1);
if (N < 2) then error('dimension must be greater than one'); end
f=sum(abs(x).^(2+3*(0:N-1)'/(N-1)));
endfunction

function f=fdiffpow4(x)
N = size(x,1);
if (N < 2) then error('dimension must be greater than one'); end
f=sum(abs(x).^(2+4*(0:N-1)'/(N-1)));
endfunction

function f=frastrigin10(x)
N = size(x,1); if (N < 2) then error('dimension must be greater than one'); end
scale=1.^((0:N-1)'/(N-1));
f = 10*size(x,1) + sum((scale.*x).^2 - 10*cos(2*%pi*(scale.*x)));
endfunction

function f=frastriginskew10(x)
N = size(x,1); if (N < 2) then error('dimension must be greater than one'); end
pop = size(x,2);
if pop > 1 then error('parallelized frastriginskew10 is not yet tested'); end
scale=10.^((0:N-1)'/(N-1));
// parallelize
scale = scale * ones(1,pop);
scale(x>0) = 1;
f = 10*size(x,1) + sum((scale.*x).^2 - 10*cos(2*%pi*(scale.*x)));
endfunction

function f=frand(x)
f=rand(1,1);
endfunction

function f=fcigarblock(x)
  N = size(x,1);
  f = fcigar(blockcoordinatesystem(N)*x);
endfunction

function f=ftabletblock(x)
  N = size(x,1);
  f = ftablet(blockcoordinatesystem(N)*x);
endfunction

function f=felliblock(x)
  N = size(x,1);
  f = felli(blockcoordinatesystem(N)*x);
endfunction

function f=fellirot(x)
  N = size(x,1);

  f = felli(coordinatesystem(N)*x);
endfunction

function f=frosenrot(x)
  f = frosen(coordinatesystem(size(x,1))*x);
endfunction

function M = coordinatesystem(N)
  global ORTHOCOORSYS_G
  if type(ORTHOCOORSYS_G) ~= 15
    ORTHOCOORSYS_G = list();
  end
  if ~or(N==definedfields(ORTHOCOORSYS_G))
    ORTHOCOORSYS_G(N) = gencoordinatesystem(N);
  end
  M = ORTHOCOORSYS_G(N);
endfunction

function ar = gencoordinatesystem(N)
  for N = N
    ar = rand(N,N,'normal');
    for i = 1:N
      for j = 1:i-1
	ar(:,i) = ar(:,i) - ar(:,i)'*ar(:,j) * ar(:,j);
      end
      ar(:,i) = ar(:,i) / norm(ar(:,i));
    end
  end
endfunction

function M = blockcoordinatesystem(N)
  global BLOCKCOORSYS_G
  if type(BLOCKCOORSYS_G) ~= 15
    BLOCKCOORSYS_G = list();
  end
  if ~or(N==definedfields(BLOCKCOORSYS_G))
    BLOCKCOORSYS_G(N) = genblockcoordinatesystem(N);
  end
  M = BLOCKCOORSYS_G(N);
endfunction

function ar = genblockcoordinatesystem(N)
  for N = N
    ar = eye(N,N);
    for i = 1:2:N-1
      ar(i:i+1,i:i+1) = 1./sqrt(2);
      ar(i, i) = -ar(i, i);
    end
  end
endfunction


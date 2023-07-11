function yout = newtonPolyVal(x, u, xout)
% This function builds the Newton interpolant and evaluates it at xout
%
% INPUT: 
% x: mesh points to be used to build the interpolant.
% u: divided difference needed to build the interpolant.
% xout: where we wish to evaluate the interpolant.
%
% OUTPUT:
% yout: result of evaluating the interpolant at xout.

  n = length(x);
  yout = u(n);

  for i=n-1:-1:1
    yout  = yout * (xout-x(i)) + u(i);
  end
 
end

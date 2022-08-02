function out = newtonPolyVal(x, u, xout)
%! This function builds up the newton interpolant and evaluates it at xout
%
% INPUT: 
% x: mesh points to be used to build the interpolant.
% u: divided difference need to build the interpolant.
% xout: where we wish to evaluate  the interpolant.
%
% OUTPUT:
% yout: result of evaluating the interpolant at xout.

  n = length(x);
  yout = coef(n);

  for i=n-1:-1:1
    yout  = yout * (xout-x(i)) + u(i);
  end
 
end

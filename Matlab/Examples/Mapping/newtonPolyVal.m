function out = newtonPolyVal(x, coef, x_val)
%% x_val --> is the value to be evaluated
%% coef --> coeffician of the newton polynomial
%% xin --> vect of x values


  n = length(x);
  %out = zeros(size(x_val));
  out = coef(n);

  for i=n-1:-1:1
    out  = out * (x_val-x(i)) + coef(i);
  end
 
end

function [vout] = adaptiveInterpolation2D(x, y, v, xout, yout, degree, interpolation_type, st, eps0, eps1)
%!This routine adaptively build a 2D tensor product interpolation based adaptiveInterpolation1D(...)
% INPUT: 
% nx: the number of points in the 1D vector x.
% ny: the number of points in the 1D vector y.
% mx: the number of points in the 1D vector xout.
% my: the number of points in the 1D vector xout.
% x: 1D mesh points of length nx used to build tensor product mesh. For i=1, ..., n-1 x_{i} <  x_{i+1}
% y: 1D mesh points of length ny used to build tensor product mesh. For i=1, ..., n-1 y_{i} <  y_{i+1}
% v: 2D array that has the data values associated with the tensor product mesh obtained from x and y
% xout: 1D vector of length mx used to construct the output tensor product mesh.
% yout: 1D vector of length my used to construct the output tensor product mesh.
% Interpolation_type: used to determine the type of interpolation to be used to build interpolant.
%   - interpolation_type=1: a data-bounded interpolant is built for each interpolant.
%   - interpolation_type=2: a positivity-preserving interpolant is built for each interpolant.
% degree: target polynomial degree and maximum polynomial degree used for each interval.
% st (optional): used guide point selection process in cases when adding the next point to the 
%   right or left both meet the requirements for positivity or data-boundedness.
%   - st=1 (default): the point with the smallest divided difference is added (ENO stencil).
%   - st=2 the point to the left of the current stencil is selected if the number of points to left
%     of x_{i} is smaller than the number of points to the right of x_{i} (i-si < ei-i). Similarly, 
%     the point to the right is selected if the number of points to the right of x_{i} is smaller
%     than the number of points to the left (i-si > ei-i). When both the number of points to right 
%     and left are the same, the algorithm chooses the point with the smallest lambda.  
%   - st=3 the point that is closest to the starting interval is chosen.
% eps0 (optional): positive parameter use constrains the bound of the positive interpolant in intervals with no
%   extremum detected.
% eps1 (optional): positive parameter use constrains the bound of the positive interpolant in intervals with 
%   extremum detected.
%
% OUTPUT:
% vout: results of evaluating interpolants on tensor product mesh obtained from xout and yout.
%
 
  %% Set optional parameters
  if(exist('st'))
    sten = st;
  else
    sten = 3;
  end
  if(exist('eps0'))
    eps2 = eps0;
  else
    eps2 = 0.01;
  end
  if(exist('eps1'))
    eps3 = eps1;
  else
    eps3 = 1.0;
  end

  nx = length(x);
  ny = length(y);
  mx = length(xout);
  my = length(yout);

  vout = zeros(mx, my);
  voutx = zeros(mx,ny);

  %% 1D interpolation along x
  voutx = zeros(mx, ny);
  for j=1:ny
    voutx(:,j) = adaptiveInterpolation1D(x, v(:,j), xout, degree, interpolation_type, sten, eps2, eps3);
  end

  %% 1D interpolation along y
  vout = zeros(mx, my);
  for i=1:mx
    vout(i,:) = adaptiveInterpolation1D(y, voutx(i,:), yout, degree,  interpolation_type, sten, eps2, eps3);
  end
end 

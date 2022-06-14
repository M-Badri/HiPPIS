function [vout] = adaptiveInterpolation2D(x, y, v, xout, yout, degree, limiter)
%
% This function performs adaptive polynomial inter interpolation to estimate
% The scalar values at location (xout(i), yout(i), zout(i))
%%
%% INPUT:
%% x, y are 1D vectors that denote the spacial coordinates of the input v
%% v is a 2D array, matrix, that denote the scalar values associate with the spatial coordinate. f(x(i), y(j), z(k)) = v(i,j,k)
%% xout, yout, zout are 1D vectors that denote the spacial coordinates of the output vout
%% degree is the desired degree of the output polynomial
%% limiter indicate the type of limiter to be used for the interpolation. 
%%	limiter = 0 --> adaptive interpolation without limiter (ENO-like)
%%	limiter = 1 --> data-bounded interpolation
%%	limiter = 2 --> positivity-preserving interpolation

%% OUTPUT
%% vout is a 2D array, matrix, that denote the scalar values associate with the spatial coordinate. f(xout(i), yout(j), zout(k)) = v(i,j,k)

  nx = length(x);
  ny = length(y);
  mx = length(xout);
  my = length(yout);


  %% 1D interpolation along x
  voutx = zeros(mx, ny);
  for j=1:ny
    voutx(:,j) = adaptiveInterpolation1D(x, v(:,j), xout, degree, limiter);
  end

  %% 1D interpolation along y
  vout = zeros(mx, my);
  for i=1:mx
    vout(i,:) = adaptiveInterpolation1D(y, voutx(i,:), yout, degree, limiter);
  end
end % end of function

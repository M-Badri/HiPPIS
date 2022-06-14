function [vout] = adaptiveInterpolation3D(x, y, z, v, xout, yout, zout, degree, limiter)
%
% This function performs adaptive polynomial inter interpolation to estimate
% The scalar values at location (xout(i), yout(i), zout(i))
%%
%% INPUT:
%% x, y, z are 1D vectors that denote the spacial coordinates of the input v
%% v is a 3D array, tensor, that denote the scalar values associate with the spatial coordinate. f(x(i), y(j), z(k)) = v(i,j,k)
%% xout, yout, zout are 1D vectors that denote the spacial coordinates of the output vout
%% degree is the desired degree of the output polynomial
%% limiter indicate the type of limiter to be used for the interpolation. 
%%	limiter = 0 --> adaptive interpolation without limiter (ENO-like)
%%	limiter = 1 --> data-bounded interpolation
%%	limiter = 2 --> positivity-preserving interpolation

%% OUTPUT
%% vout is a 3D array, tensor, that denote the scalar values associate with the spatial coordinate. f(xout(i), yout(j), zout(k)) = v(i,j,k)

  nx = length(x);
  ny = length(y);
  nz = length(z);
  mx = length(xout);
  my = length(yout);
  mz = length(zout);


  %% 1D interpolation along x
  voutx = zeros(mx, ny, nz);
  for j=1:ny
    for k=1:nz
      voutx(:,j,k) = adaptiveInterpolation1D(x, v(:,j,k), xout, degree, limiter);
    end
  end

  %% 1D interpolation along y
  vouty = zeros(mx, my, nz);
  for i=1:mx
    for k=1:nz
      vouty(i,:,k) = adaptiveInterpolation1D(y, voutx(i,:,k), yout, degree, limiter);
    end
  end

  %% 1D interpolation along z
  for i=1:mx
    for j=1:my
      vout(i,j,:) = adaptiveInterpolation1D(z, vouty(i,:,k), zout, degree, limiter);
    end
  end

end % end of function

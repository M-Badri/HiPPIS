function [vout] = adaptiveInterpolation3D(x, y, z, v, xout, yout, zout, degree, interpolation_type, st, eps0, eps1)
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

  if(exist('st'))
    sten = st;
  else
    sten = 1;
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
  nz = length(z);
  mx = length(xout);
  my = length(yout);
  mz = length(zout);


  %% 1D interpolation along x
  voutx = zeros(mx, ny, nz);
  for j=1:ny
    for k=1:nz
      voutx(:,j,k) = adaptiveInterpolation1D(x, v(:,j,k), xout, degree, interpolation_type, sten, eps2, eps3);
    end
  end

  %% 1D interpolation along y
  vouty = zeros(mx, my, nz);
  for i=1:mx
    for k=1:nz
      vouty(i,:,k) = adaptiveInterpolation1D(y, voutx(i,:,k), yout, degree, interpolation_type, sten, eps1, eps3);
    end
  end

  %% 1D interpolation along z
  for i=1:mx
    for j=1:my
      vout(i,j,:) = adaptiveInterpolation1D(z, vouty(i,:,k), zout, degree, interpolation_type, sten, eps2, eps3);
    end
  end

end % end of function

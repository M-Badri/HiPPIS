

clear;
close all;
clc;

%---------------------------------------------------------------------------------------------%
% Tutorial showing how to use the 1D, 2D, and 3D DBI and PPI interpolation
%---------------------------------------------------------------------------------------------%

  mex FFLAGS='-fexceptions -fbackslash -fPIC -fno-omit-frame-pointer -Wno-argument-mismatch -fdefault-integer-8'...
               mod_adaptiveInterpolation.F90 -v
  mex FFLAGS='-fexceptions -fbackslash -fPIC -fno-omit-frame-pointer -Wno-argument-mismatch -fdefault-integer-8'...
              adaptiveInterpolation1D.F90 mod_adaptiveInterpolation.F90 -v

  %-- 1D Tutorial -- %
  x = linspace(-pi, pi, 17);       % input mesh points
  v = sin(x);                      % input data values
  xout = linspace(-pi, pi, 100);   % output points
  
  d = 8;                           % target and maximum polynomial degree used for each interval
  interpolation_type = 2;          % 1 for DBI and 2 for PPI
  sten = 1;                        % optional parameter to guide stencil selection 1, 2, and 3
  eps0 = 0.01;                     % optional positive parameter to bound interpolant in PPI
  eps1 = 1.0;                      % optional positive parameter to bound interpolant in PPI

  [vout_apprx, deg] = adaptiveInterpolation1D(x, v, xout, d, interpolation_type, sten, eps0, eps1 ); 
  %vout_apprx = adaptiveInterpolation1D(x, v, xout, d, interpolation_type, sten, eps0, eps1 ); 
  

%  %-- 2D Tutorial -- %
%  n = 17;			     % number of input points
%  x = linspace(-pi, pi, n);          % input mesh points
%  y = x;                             
%  v2D = zeros(n,n);
%  for j=1:n
%    for i=1:n
%      v2D(i,j) = sin(x(i))*sin(y(j));  % input data values
%    end
%  end
%  xout = linspace(-pi, pi, 100);     % output points
%  yout = xout;
%  d = 8;                             % target and maximum polynomial degree used for each interval
%  interpolation_type = 2;            % 1 for DBI and 2 for PPI
%  sten = 1;                          % optional parameter to guide stencil selection 1, 2, and 3
%  eps0 = 0.01;                       % optional positive parameter to bound interpolant in PPI
%  eps1 = 1.0;                        % optional positive parameter to bound interpolant in PPI
%
%  vout_apprx2D = adaptiveInterpolation2D(x, y, v2D, xout,yout, d, interpolation_type, sten, eps0, eps1 ); 
%  
%
%  %-- 3D Tutorial -- %
%  n = 17;			     % number of input points
%  x = linspace(-pi, pi, n);          % input mesh points
%  y = x;
%  z = x;
%  v3D = zeros(n,n,n);
%  for k=1:n
%    for j=1:n
%      for i=1:n
%        v3D(i,j,k) = sin(x(i))*sin(y(j))*sin(z(k));  % input data values
%      end
%    end
%  end
%  xout = linspace(-pi, pi, 100);     % output points
%  yout = xout;
%  zout =xout;
%  d = 8;                             % target and maximum polynomial degree used for each interval
%  interpolation_type = 2;            % 1 for DBI and 2 for PPI
%  sten = 1;                          % optional parameter to guide stencil selection 1, 2, and 3
%  eps0 = 0.01;                       % optional positive parameter to bound interpolant in PPI
%  eps1 = 1.0;                        % optional positive parameter to bound interpolant in PPI
%
%  vout_apprx3D = adaptiveInterpolation3D(x, y, z, v3D, xout, yout, zout, d, interpolation_type, sten, eps0, eps1 ); 
 



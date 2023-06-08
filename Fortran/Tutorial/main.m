

clear;
close all;
clc;

%---------------------------------------------------------------------------------------------%
% Tutorial showing how to use the 1D, 2D, and 3D DBI and PPI interpolation
%---------------------------------------------------------------------------------------------%

  mex FFLAGS='-fexceptions -fbackslash -fPIC -fno-omit-frame-pointer -Wno-argument-mismatch -fdefault-integer-8'...
              mod_adaptiveInterpolation.F90 
  mex FFLAGS='-fexceptions -fbackslash -fPIC -fno-omit-frame-pointer -Wno-argument-mismatch -fdefault-integer-8'...
              adaptiveInterpolation1D_vec.F90 mod_adaptiveInterpolation.F90
  mex FFLAGS='-fexceptions -fbackslash -fPIC -fno-omit-frame-pointer -Wno-argument-mismatch -fdefault-integer-8'...
              adaptiveInterpolation2D_vec.F90 mod_adaptiveInterpolation.F90 
  mex FFLAGS='-fexceptions -fbackslash -fPIC -fno-omit-frame-pointer -Wno-argument-mismatch -fdefault-integer-8'...
              adaptiveInterpolation3D_vec.F90 mod_adaptiveInterpolation.F90 
  mex FFLAGS='-fexceptions -fbackslash -fPIC -fno-omit-frame-pointer -Wno-argument-mismatch -fdefault-integer-8'...
              adaptiveInterpolation1D.F90 mod_adaptiveInterpolation.F90
  mex FFLAGS='-fexceptions -fbackslash -fPIC -fno-omit-frame-pointer -Wno-argument-mismatch -fdefault-integer-8'...
              adaptiveInterpolation2D.F90 mod_adaptiveInterpolation.F90 
  mex FFLAGS='-fexceptions -fbackslash -fPIC -fno-omit-frame-pointer -Wno-argument-mismatch -fdefault-integer-8'...
              adaptiveInterpolation3D.F90 mod_adaptiveInterpolation.F90 



clear;
close all;
clc;

  
  %-- 1D Tutorial -- %
  n = 17;
  m = 101;
  x = linspace(-1.0, 1.0, n);       % input mesh points
  v = 0.1./(0.1 + 25.0*x.^2);                      % input data values
  xout = linspace(-1.0, 1.0, m);   % output points
  vt = 0.1./(0.1 + 25.0*xout.^2);                      % true solution data values
  
  d = 8;                           % target and maximum polynomial degree used for each interval
  interpolation_type = 2;          % 1 for DBI and 2 for PPI
  sten = 1;                        % optional parameter to guide stencil selection 1, 2, and 3
  eps0 = 0.01;                     % optional positive parameter to bound interpolant in PPI
  eps1 = 1.0;                      % optional positive parameter to bound interpolant in PPI

  %[vout_apprx, deg] = adaptiveInterpolation1D(x, v, xout, d, interpolation_type, sten, eps0, eps1 ); 
  [vout_apprx, deg] = adaptiveInterpolation1D_vec(x, v, xout, d, interpolation_type, sten, eps0, eps1 ); 
  

  %-- Plot approximated results --%
  figure 
  plot(xout, vt, xout, vout_apprx)
  xlabel('x')
  ylabel('y')
  legend('True', 'Approx.')
  title('Approximation of $$f_{1}(x) = \frac{0.1}{0.1 + 25x^{2}}$$', 'Interpreter', 'Latex')

  %-- 2D Tutorial -- %
  y = x;                             
  v2D = zeros(n,n);
  for j=1:n
    for i=1:n
      v2D(i,j) = 0.1/(0.1 + 25.0*(x(i)^2 + y(j)^2));  % input data values
    end
  end
  m =33;
  vt2D = zeros(m,m);
  xout = linspace(-1.0, 1.0, m);   % output points
  yout = xout;
  for j=1:m
    for i=1:m 
      vt2D(i,j) = 0.1/(0.1 + 25.0*(xout(i)^2 + yout(j)^2));  % true solution data values
    end
  end
 
  d = 8;                             % target and maximum polynomial degree used for each interval
  interpolation_type = 2;            % 1 for DBI and 2 for PPI
  sten = 1;                          % optional parameter to guide stencil selection 1, 2, and 3
  eps0 = 0.01;                       % optional positive parameter to bound interpolant in PPI
  eps1 = 1.0;                        % optional positive parameter to bound interpolant in PPI

  %vout_apprx2D = adaptiveInterpolation2D(x, y, v2D, xout,yout, d, interpolation_type, sten, eps0, eps1 ); 
  vout_apprx2D = adaptiveInterpolation2D_vec(x, y, v2D, xout,yout, d, interpolation_type, sten, eps0, eps1 ); 
  
  %-- Plot approximated results --%
  [xx, yy] = meshgrid(xout, yout);
  figure 
  subplot(1,2,1)
  surf(xx, yy, vt2D)
  xlabel('x')
  ylabel('y')
  zlabel('z')
  title('$$f_{1}(x) = \frac{0.1}{0.1 + 25(x^{2} + y^{2})}$$', 'Interpreter', 'Latex')
  subplot(1,2,2)
  surf(xx, yy, vout_apprx2D)
  xlabel('x')
  ylabel('y')
  zlabel('z')
  title('Approximation of $$f_{1}(x) = \frac{0.1}{0.1 + 25(x^{2}+y^{2})}$$', 'Interpreter', 'Latex')


  %-- 3D Tutorial -- %
  y = x;
  z = x;
  v3D = zeros(n,n,n);
  for k=1:n
    for j=1:n
      for i=1:n
        v3D(i,j,k) = 0.1/(0.1 + 25.0*(x(i)^2 + y(j)^2 + z(k)^2));  % input data values
      end
    end
  end
  xout = linspace(-1.0, 1.0, m);     % output points
  yout = xout;
  zout = xout;
  d = 8;                             % target and maximum polynomial degree used for each interval
  interpolation_type = 2;            % 1 for DBI and 2 for PPI
  sten = 1;                          % optional parameter to guide stencil selection 1, 2, and 3
  eps0 = 0.01;                       % optional positive parameter to bound interpolant in PPI
  eps1 = 1.0;                        % optional positive parameter to bound interpolant in PPI

  %vout_apprx3D = adaptiveInterpolation3D(x, y, z, v3D, xout, yout, zout, d, interpolation_type, sten, eps0, eps1 ); 
  vout_apprx3D = adaptiveInterpolation3D_vec(x, y, z, v3D, xout, yout, zout, d, interpolation_type, sten, eps0, eps1 ); 
 



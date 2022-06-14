%
% This script contains exmaple demonstrating how the use the 
% adaptiveInterpolation function
%

clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This first part of the script considers a 1D example v(x) = sinx(x).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Define the functions to be used for the examples
%
sin1D = @(xx) sin(xx);


fprintf('1D Example Adp. Interp. \n')

np = 17;    % number of points used for approximations

%
% Type of limter used for the interpolation
%
limiter =2 ;

%
% maximum polynomial degree used for each interval
%
degree = 4;

%
% number of points to be interpolated to
%
mp = 100;

%
% 1D input mesh points 
%
x = linspace(-pi/2, pi/2, np);
x = x';
nx = np;
v1D = sin1D(x);  

%
% 1D output points 
%
xout = linspace(-pi/2, pi/2, mp);
xout = xout';
v1D_out_true = sin1D(xout);  

%
% Interpolation from x to xout
%
v1D_out = adaptiveInterpolation1D(x, v1D', xout, degree, limiter);

%
% Interpolation using PCHIP
%
v1D_pchip = pchip(x, v1D', xout);

figure; 
plot(xout, v1D_out_true, 'k', xout, v1D_out, xout, v1D_pchip', 'LineWidth', 4 );
legend('True', 'Interp.');
xlabel('x');
ylabel('y');
title('1D example ( sin(x) )');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This second part of the script considers a 2D example v(x,y) = sinx(x)sin(y).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define the functions to be used for the examples
%
sin2D = @(xx, yy) sin(xx).*sin(yy);

fprintf('2D Example Adp. Interp. \n')

np = 17;    % number of points used for approximations

%
% Type of limter used for the interpolation
%
limiter =2 ;

%
% maximum polynomial degree used for each interval
%
degree = 4;

%
% number of points to be interpolated to
%
mp = 100;

%
% 2D input mesh points 
%
x = linspace(-pi/2, pi/2, np);
x = x';
y = x;
nx = np;
ny = np;
v2D = sin2D(x, y);  

%
% 2D output points 
%
xout = linspace(-pi/2, pi/2, mp);
xout = xout';
yout = xout;

%
% Create meshes and data values associted with 
% meshes
%
[xmesh, ymesh] = meshgrid(x, y);
[xoutmesh, youtmesh] = meshgrid(xout, yout);
v2D = sin2D(xmesh, ymesh);  
v2D_out_true = sin2D(xoutmesh, youtmesh);  

%
% Interpolate from 2D mesh meshgrid(x,y) to the mesh meshgrid(xout,yout)
% 
v2D_out = adaptiveInterpolation2D(x, y, v2D, xout, yout, degree, limiter);

figure
subplot(1,2,1)
surf(xoutmesh, youtmesh, v2D_out_true)
xlabel('x')
ylabel('y')
zlabel('v')
title('2D Example (True solution) ')
colorbar
%
subplot(1,2,2)
surf(xoutmesh, youtmesh, v2D_out)
xlabel('x')
ylabel('y')
zlabel('v')
title('2D example (Approximated solution)')
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This third part of the script considers a 3D example v(x,y,z) = sinx(x)sin(y)sin(z).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define the functions to be used for the examples
%
sin3D = @(xx, yy, zz) sin(xx).*sin(yy).*sin(zz);


fprintf('3D Example Adp. Interp. \n')

np = 17;    % number of points used for approximations

%
% Type of limter used for the interpolation
%
limiter =2 ;

%
% maximum polynomial degree used for each interval
%
degree = 4;

%
% number of points to be interpolated to
%
mp = 30;

%
% 2D input mesh points 
%
x = linspace(-pi/2, pi/2, np);
x = x';
y = x;
z = x;
nx = np;
ny = np;
nz = np;
v2D = sin3D(x, y, z);  

%
% 2D output points 
%
xout = linspace(-pi/2, pi/2, mp);
xout = xout';
yout = xout;
zout = xout;

%
% create 3d meshes
%
[xmesh, ymesh, zmesh] = meshgrid(x,y,z);
[xoutmesh, youtmesh, zoutmesh] = meshgrid(xout,yout,zout);

%
% Evaluate function at mesh values
%
v3D = sin3D(xmesh, ymesh, zmesh);  
%v3D_out_true = sin3D(xoutmesh, youtmesh, zoutmesh);  

v3D_out = adaptiveInterpolation3D(x, y, z, v3D, xout, yout, zout, degree, limiter);


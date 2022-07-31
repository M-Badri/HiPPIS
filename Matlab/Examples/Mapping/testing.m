clear
close all
clc

%
%
%

  %test1()
  %test2()
  %test3()
  %test4()
  %test5()
  %test6()
  test7()
%---------------------------------------------------------------------------------------------%
% test1
%---------------------------------------------------------------------------------------------%

function test1()
%
%  test considering small intervals for input and output mesh
%

  n = 1e+6;
  m = 1e+7;
  a = -1.0e-10;
  b = 1.0e-10;
  x= linspace(a, b, n);
  xout= linspace(a, b, m);
  v1D = zeros(n,1);
  v1Dout = zeros(m,1);
  eps0= 1;
  eps1= 1;
  dxn = 1; % dummy variable 
  d = 8;
  for j=1:3
    sten = j
    for i=1:n
      v1D(i) = 0.1 /(0.1 + 25*(x(i)*1e+10)^2);
    end
    
    v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, sten, eps0, eps1); 

    figure
    plot(x, v1D, '*', xout, v1Dout)
    legend('data', 'apprx')
    pause
  end
end 


function test2()
%
%  test considering small example with negative values
%

  n = 20;
  m = 1000;
  a = -1.0;
  b = 1.0;
  x= linspace(a, b, n);
  xout= linspace(a, b, m);
  v1D = zeros(n,1);
  v1Dout = zeros(m,1);
  eps0= 1;
  eps1= 1;
  dxn = 1; % dummy variable 
  d = 8;
  for j=1:3
    sten = j
    for i=1:n
      v1D(i) = sin(x(i)*pi);
    end
    
    v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, sten, eps0, eps1); 

    figure
    plot(x, v1D, '*', xout, v1Dout)
    legend('data', 'apprx')
    pause
  end
end 


function test3()
%
%  test considering small constant fucntions
%

  n = 20;
  m = 1000;
  a = -1.0;
  b = 1.0;
  x= linspace(a, b, n);
  xout= linspace(a, b, m);
  v1D = zeros(n,1);
  v1Dout = zeros(m,1);
  eps0= 1;
  eps1= 1;
  dxn = 1; % dummy variable 
  d = 8;
  for j=1:3
    sten = j
    for i=1:n
      v1D(i) = 1.0;
    end
    v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, sten, eps0, eps1); 
    figure
    plot(x, v1D, '*', xout, v1Dout)
    legend('data', 'apprx')
    pause
  end
end 

function test4()
%
%  test considering a linear function constant fucntions
%

  n = 20;
  m = 1000;
  a = -1.0;
  b = 1.0;
  x= linspace(a, b, n);
  xout= linspace(a, b, m);
  v1D = zeros(n,1);
  v1Dout = zeros(m,1);
  eps0= 1;
  eps1= 1;
  dxn = 1; % dummy variable 
  d = 8;
  for j=1:3
    sten = j
    for i=1:n
      v1D(i) = (-1.0/(b-a))*(x(i)+1) + 1;
    end
    v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, sten, eps0, eps1); 
    figure
    plot(x, v1D, '*', xout, v1Dout)
    legend('data', 'apprx')
    pause
  end
end 

function test5()
%
%  test commapring L^{2}-norm error for st=1, 2 and 3 with Runge
%

  m = 10000;
  a = -1.0;
  b = 1.0;
  for n= [17, 33, 65, 129 257] ;                                              
    fprintf('---  n = %d  -- \n', n);
    x= linspace(a, b, n);
    xout= linspace(a, b, m);
    v1D = zeros(n,1);
    v1Dout = zeros(m,1);
    v1Dout_true = zeros(m,1);
    eps0= 1;  % set to default value
    eps1= 0.01; % set to default value 

    d = 8;
    fprintf('st      error \n')
    for j=1:3
      st = j;
      for i=1:n
        v1D(i) = 0.1/(0.1 + 25.0*x(i)*x(i));
      end
      for i=1:m
        v1Dout_true(i) = 0.1/(0.1 + 25.0*xout(i)*xout(i));
      end
      v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, st); 
      err = sqrt( trapz(xout, (v1Dout-v1Dout_true').^2));
      fprintf('%d \t %.8E \n', j, err)
    end
  end
  fprintf('No significant different between st=1, 2, and 3 \n')

end 

function test6()
%
%  test 2D examples with different number of points for each dimension
%

  mx = 100;
  my = 200;
  nx = 17;
  ny = 33;

  ax = -1.0;
  bx = 1.0;
  ay = -1.0;
  by =  1.0;
  x= linspace(ax, bx, nx)';
  y= linspace(ay, by, ny)';
  xout= linspace(ax, bx, mx)';
  yout= linspace(ay, by, my)';
  v2D = zeros(nx,ny);
  v2Dout = zeros(mx,my);
  d = 8;
  for j=1:ny
   for i=1:nx
     v2D(i,j) = 0.1/(0.1 + 25.0*(x(i)*x(i) + y(j)*y(j)));
   end
  end
  %
  v2Dout = adaptiveInterpolation2D(x, y,  v2D, xout, yout, d, 2); 

  [xx, yy] = meshgrid(yout, xout);
  surf(xx, yy, v2Dout)
end 



function test7()
%
%  test 2D examples with different number of points for each dimension
%

  mx = 100;
  my = 200;
  nx = 17;
  ny = 33;

  ax = -1.0e-10;
  bx = 1.0e-10;
  ay = -1.0;
  by =  1.0;
  x= linspace(ax, bx, nx)';
  y= linspace(ay, by, ny)';
  xout= linspace(ax, bx, mx)';
  yout= linspace(ay, by, my)';
  v2D = zeros(nx,ny);
  v2Dout = zeros(mx,my);
  d = 8;
  for j=1:ny
   for i=1:nx
     v2D(i,j) = 0.1/(0.1 + 25.0*(x(i)*1.0e+10*x(i)*1.0e+10 + y(j)*y(j)));
   end
  end
  %
  v2Dout = adaptiveInterpolation2D(x, y,  v2D, xout, yout, d, 2); 
  [xx, yy] = meshgrid(yout, xout);
  surf(xx, yy, v2Dout)

end 



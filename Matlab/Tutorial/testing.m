clear
close all
clc

%
%

  test1()
  test2()
  test3()
  test4()
  test5()
  test6()
  test7()

function test1()
%
%  Test example considering small intervals for input and output mesh
%  This test verifies that small values used for the input 
%  mesh points do not cause large oscillations or failure
%  of the DBI and PPI methods.
%

  n = 1e+6;
  m = 1e+7;
  a = -1.0e-10;
  b = 1.0e-10;
  check = 1;
  x= linspace(a, b, n);
  xout= linspace(a, b, m);
  v1D = zeros(n,1);
  v1Dout = zeros(m,1);
  v1Dout_true = zeros(m,1);
  d = 8;
  for i=1:n
    v1D(i) = 0.1 /(0.1 + 25*(x(i)*1e+10)^2);
  end
  for i=1:m
    v1Dout_true(i) = 0.1 /(0.1 + 25*(xout(i)*1e+10)^2);
  end
  for j=1:3
    sten = j;
    v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2); 

    for i=1:m
      if(abs(v1Dout(i)-v1Dout_true(i)) > 1.0e-1)
         fprintf('test1() FAILED: the difference betewn approximated and true solution at %d i is %f \n', ...
                  i, abs(v1Dout(i)-v1Dout_true(i)))
         check = 0;
      end
    end
    %figure
    %plot(x, v1D, '*', xout, v1Dout)
    %legend('data', 'apprx')
  end
  if(check ==1)
    fprintf('test1() --- Passed --- \n')
  end
end 


function test2()
%
%  Test considering examples with negative values.
%  This verifies that the DBI and PPI can approximate
%  functions that are not positive.
%
  n = 20;
  m = 1000;
  a = -1.0;
  b = 1.0;
  check = 1;
  x= linspace(a, b, n);
  xout= linspace(a, b, m);
  v1D = zeros(n,1);
  v1Dout = zeros(m,1);
  v1Dout_true = zeros(m,1);
  eps0= 1;
  eps1= 1;
  d = 8;
  for i=1:n
    v1D(i) = sin(x(i)*pi);
  end
  for i=1:m
    v1Dout_true(i) = sin(xout(i)*pi);
  end
  for j=1:3
    sten = j;
    v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, sten, eps0, eps1); 
    for i=1:m
      if(abs(v1Dout(i)-v1Dout_true(i)) > 1.0e-1)
         fprintf('test2() FAILED: the difference betewn approximated and true soluiton at %d i is %f \n', ...
                  i, abs(v1Dout(i)-v1Dout_true(i)))
         check = 0;
      end
    end
 
    %figure
    %plot(x, v1D, '*', xout, v1Dout)
    %legend('data', 'apprx')
  end
  if(check ==1)
    fprintf('test2() --- Passed --- \n')
  end
end 


function test3()
%
%  Test considering a constant function.
%  This test verifies that the DBI and PPI methods
%  are able to approximate linear functions
%


  n = 20;
  m = 1000;
  a = -1.0;
  b = 1.0;
  check = 1;
  x= linspace(a, b, n);
  xout= linspace(a, b, m);
  v1D = ones(n,1);
  v1Dout = zeros(m,1);
  v1Dout_true = ones(m,1);
  d = 8;
  for j=1:3
    sten = j;
    v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, sten); 
    for i=1:m
      if(abs(v1Dout(i)-v1Dout_true(i)) > 1.0e-10)
         fprintf('test2() FAILED: the difference betewn approximated and true soluiton at %d i is %f \n', ...
                  i, abs(v1Dout(i)-v1Dout_true(i)))
         check = 0;
      end
    end

    %figure
    %plot(x, v1D, '*', xout, v1Dout)
    %legend('data', 'apprx')
    %pause
  end
  if(check ==1)
    fprintf('test3() --- Passed --- \n')
  end

end 

function test4()
%
%  test considering a linear function.
%  This test verifies that the DBI and PPI methods
%  are able to approximate linear functions
%

  n = 20;
  m = 1000;
  a = -1.0;
  b = 1.0;
  check = 1;
  x= linspace(a, b, n);
  xout= linspace(a, b, m);
  v1D = zeros(n,1);
  v1Dout = zeros(m,1);
  v1Dout_true = zeros(m,1);
  d = 8;
  for i=1:n
    v1D(i) = (-1.0/(b-a))*(x(i)+1.0) + 1.0;
  end
  for i=1:m
    v1Dout_true(i) = (-1.0/(b-a))*(xout(i)+1.0) + 1.0;
  end
  for j=1:3
    sten = j;
    v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, sten); 
    for i=1:m
      if(abs(v1Dout(i)-v1Dout_true(i)) > 1.0e-10)
         fprintf('test4() FAILED: the difference betewn approximated and true soluiton at %d i is %f \n', ...
                  i, abs(v1Dout(i)-v1Dout_true(i)))
         check = 0;
      end
    end

    %figure
    %plot(x, v1D, '*', xout, v1Dout_true, xout, v1Dout)
    %legend('data', 'true', 'apprx')
    %pause
  end
  if(check ==1)
    fprintf('test4() --- Passed --- \n')
  end
end 

function test5()
%
%  test commapring L^{2}-norm error for st=1, 2 and 3 with Runge
%

  m = 10000;
  a = -1.0;
  b = 1.0;
  check =1; 
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
    for i=1:n
      v1D(i) = 0.1/(0.1 + 25.0*x(i)*x(i));
    end
    for i=1:m
      v1Dout_true(i) = 0.1/(0.1 + 25.0*xout(i)*xout(i));
    end
 
    fprintf('st      error \n')
    tmp = zeros(3,1);
    for j=1:3
      st = j;
      v1Dout = adaptiveInterpolation1D(x, v1D, xout, d, 2, st); 
      err = sqrt( trapz(xout, (v1Dout-v1Dout_true').^2));
      tmp(j) = err;
      fprintf('%d \t %.8E \n', j, err)
    end
    if(abs(tmp(1)-tmp(2)) >10 || abs(tmp(3)-tmp(2)) >10 || abs(tmp(1)-tmp(3)) >10)
      fprintf('test5() FAILED')
      check = 0;
    end
  end
  if(check ==1)
    fprintf('No significant different between st=1, 2, and 3 \n')
    fprintf('test5() --- Passed --- \n')
  end

end 

function test6()
%
%  test 2D examples with different number of points for each dimension.
%  This test verifies that the DBI and PPI work when for 2D examples 
%  with different number points in each direction i.e. nx not equal to ny 
% and mx not equal to my 
%  

  mx = 100;
  my = 200;
  nx = 33;
  ny = 65;
  ax = -1.0;
  bx = 1.0;
  ay = -1.0;
  by =  1.0;
  check = 1;
  x= linspace(ax, bx, nx)';
  y= linspace(ay, by, ny)';
  xout= linspace(ax, bx, mx)';
  yout= linspace(ay, by, my)';
  v2D = zeros(nx,ny);
  v2Dout = zeros(mx,my);
  v2Dout_true = zeros(mx,my);

  d = 8;
  for j=1:ny
   for i=1:nx
     v2D(i,j) = 0.1/(0.1 + 25.0*(x(i)*x(i) + y(j)*y(j)));
   end
  end
  %
  for j=1:my
   for i=1:mx
     v2Dout_true(i,j) = 0.1/(0.1 + 25.0*(xout(i)*xout(i) + yout(j)*yout(j)));
   end
  end
 
  v2Dout = adaptiveInterpolation2D(x, y,  v2D, xout, yout, d, 2); 

  for j=1:my
   for i=1:mx
     if(abs(v2Dout(i,j)-v2Dout_true(i,j)) > 1.0e-1)
         fprintf('test6() FAILED: the difference betewn approximated and true soluiton at %d and %d is %f \n', ...
                  i, j, abs(v2Dout(i,j)-v2Dout_true(i,j)))
         check = 0;
      end
    
   end
  end
  %[xx, yy] = meshgrid(yout, xout);
  %surf(xx, yy, v2Dout)
  if(check ==1)
    fprintf('test4() --- Passed --- \n')
  end
end 

function test7()
%
%  test 2D examples with different number of points for each dimension.
%  with small input mesh values. This test verifies that the DBI and PPI 
%  works when for 2D examples with different number points in each direction 
%  i.e. nx not equal to ny and mx not equal to my 
%  

  mx = 100;
  my = 200;
  nx = 33;
  ny = 65;
  ax = -1.0e-10;
  bx = 1.0e-10;
  ay = -1.0;
  by =  1.0;
  check = 1;
  x= linspace(ax, bx, nx)';
  y= linspace(ay, by, ny)';
  xout= linspace(ax, bx, mx)';
  yout= linspace(ay, by, my)';
  v2D = zeros(nx,ny);
  v2Dout = zeros(mx,my);
  v2Dout_true = zeros(mx,my);

  d = 8;
  for j=1:ny
   for i=1:nx
     v2D(i,j) = 0.1/(0.1 + 25.0*(x(i)*1.0e+10*x(i)*1.0e+10 + y(j)*y(j)));
   end
  end
  %
  for j=1:my
   for i=1:mx
     v2Dout_true(i,j) = 0.1/(0.1 + 25.0*(xout(i)*1.0e+10*xout(i)*1.0e+10 + yout(j)*yout(j)));
   end
  end
 
  v2Dout = adaptiveInterpolation2D(x, y,  v2D, xout, yout, d, 2); 

  for j=1:my
   for i=1:mx
     if(abs(v2Dout(i,j)-v2Dout_true(i,j)) > 1.0e-1)
         fprintf('test6() FAILED: the difference betewn approximated and true soluiton at %d and %d is %f \n', ...
                  i, j, abs(v2Dout(i,j)-v2Dout_true(i,j)))
         check = 0;
      end
    
   end
  end
  %[xx, yy] = meshgrid(yout, xout);
  %surf(xx, yy, v2Dout)
  if(check ==1)
    fprintf('test7() --- Passed --- \n')
  end
end 

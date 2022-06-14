%%
%% This is a script that test the dfferent functions in the the PPI folder.

clear
close all


n = 11;				% number of mesh points
d = 10; 			% polynomial degree
m = 100;				% number of points outside	
a = -10 + (10+10)*rand(d+1,1);  % polynomial coeficients
x = linspace(-1,1, n);		% input mesh points
xout = linspace(-1,1, m)
y = polyval(a, x);		% input data values
y_interp_true = polyval(a, xout);
y_interp = zeros(size(xout));


table = divdiff(x, y, d);
for i=1:m
  y_interp(i) = newtonPolyval(x, table(1,:), xout(i))
end
f1 = figure

figure(f1)
plot(x, y, 'x', xout, y_interp_true, 'b', xout, y_interp, 'r', 'LineWidth', 4)
legend('data', 'True', 'Interp')
xlabel('x')
ylabel('y')

err_L2 = norm(y_interp-y_interp_true)/length(xout);

  if(abs(err_L2) < 1e-10)
    fprintf('PASSED: dividff(...), and newtonPolyval(..) passed the test. The L2 norm between the true solution and the interpolated solution err = %d \n', err_L2)
  else

    fprintf('FAILED: dividff(...), and newtonPolyval(..) failed the test. The L2 norm between the true solution and the interpolated solution err = %d \n', err_L2)
  end

clear
close all
clc


%!!** experiment for math paper and report **!!
%paperExamples1D()
paperExamples2D()


%
%
%

function paperExamples1D()
%%
%%
%%

  n = [17; 33; 65; 129; 257];                                              %!! array with number of inputpoints
  d = [1; 4; 8; 16];                                                       %!! array with interpolants degrees
  a = [-1.0; -0.2; -1.0; -2.0; 0.0];
  b = [1.0;  0.2;  1.0; 0.0; 1.0];
  eps0 = [ 1.0;  0.1;  0.01; 0.001; 0.0001; 0.00 ];
  m = 10000;
  %!!** function type 1=runge funtion , 2= heaviside, 3=Gelb Tanner **!! 
  fun = [1; 2; 3; 4; 5];                       % function type

  %!!** Examples with Trunctated peak **!!
  testpeak(d(3), n(1)-1, a(1), b(1), m);

  for j=1:4 
    for i=1:5
      testpeak2(d(j), fun(1), n(i)-1, a(1), b(1), m);
    end
  end


  %!!** Comparing against PCHIP **!!
  for k=1:2
    fprintf( '*****  fun= %d ******** \n', fun(k));
    for i=1:5                                                              %!! loop over interpolant degrees
         test001(3, fun(k), n(i), a(k), b(k), m);
    end
  end

  for k=1:2%5 
    fprintf('*****  fun= %d ***** \n', fun(k) );
    for j=1:4                                                              %!! loop over interpolant degrees
      fprintf('*****  d= %d ***** \n', d(j) );
      for i=1:5                                                           %!! loop over number of input points
         test001(d(j), fun(k), n(i), a(k), b(k), m);
    end  %  !! of i loop
    end  % !! of j loop
  end %!! of k loop

end %function


function testpeak(d, n, a, b, m)

   
  %!!** Initialize variables **!!
  degOut = zeros(n-1, 2);
  x = zeros(n, 1);
  v1D = zeros(n,1);
  v1Dout = zeros(m,1);
  v1Dout_true = zeros(m,1);


  %!!** calculates intreval sizes **!!
  dxn = (b-a) / double(n-1);
  dxm = (b-a) /double(m-1);
  

  %!!** uniform mesh **!!
  for i=1:n
    x(i) = a + double(i-1)*dxn;
  end

  %!!** output mesh points **!
  for i=1:m
    xout(i) = a + double(i-1)*dxm;
  end

  %!!** Data values associated to input meshes **!!
  for i=1:n
    v1D(i) = evalFun1D(1, x(i), dxn);
  end

  %!!** True solution **!!
  for i=1:m
    v1Dout_true(i) = evalFun1D(1, xout(i), dxn);
  end

  %!!** interpolation using PCHIP **!!
  v1Dout = pchip(x, v1D, xout);

  %!!** open file **!!
  fid = fopen('peakPCHIP', 'w');
  %!!** write to file **!!
  for i=1:m
    fprintf(fid,'%.8E \t  %.8E \t  %.8E \n', xout(i), v1Dout_true(i), v1Dout(i) );
  end
  %!!** close file **!!
  fclose(fid);


  [v1Dout, degOut(1:n-1, 1)] =adaptiveInterpolation1D(x, v1D, xout, d, 1) ;

  %plot(x, v1D, xout, v1Dout_true, xout, v1Dout)
  %pause
  %!!** open file **!!
  %fid = 10
  fid = fopen('peakDBI', 'w');
  %!!** write to file **!!
  for i=1:m
    fprintf(fid,'%.8E \t %.8E \t %.8E \n', xout(i), v1Dout_true(i), v1Dout(i) );
  end
  %!!** close file **!!
  fclose(fid);
  [v1Dout, degOut(1:n-1, 2)] =adaptiveInterpolation1D(x, v1D, xout, d, 2) ;

  %!!** open file **!!
  fid = fopen('peakPPI', 'w');
  %!!** write to file **!!
  for i=1:m
    fprintf(fid,'%.8E \t %.8E \t %.8E \n', xout(i), v1Dout_true(i), v1Dout(i) );
  end
  %!!** close file **!!
  fclose(fid);

  %!!** open file **!!
  fid = fopen('peakDeg', 'w');
  %!!** write to file **!!
  for i=1:n-1
    fprintf(fid,'%d \t %d \n', degOut(i, 1), degOut(i, 2)) ;
  end
  %!!** close file **!!
  fclose(fid);


end %function

function testpeak2(d, fun, n, a, b, m)
%!! 
%!! Test interpolation with 1D runge function
%!! 

  nn = n;
  
  %!!** Initialize variables **!!
  degOut = zeros(n-1,2);
  x = zeros(n,1);
  v1D = zeros(n,1);
  v1Dout = zeros(m,1);
  v1Dout_true = zeros(m,1);


  %!!** calculates intreval sizes **!!
  dxn = (b-a) /(n-1);
  dxm = (b-a) /(m-1);
  

  %!!** uniform mesh **!!
  for i=1:n
    x(i) = a + (i-1)*dxn;
  end

  %!!** output mesh points **!
  for i=1:m
    xout(i) = a + (i-1)*dxm;
  end

  %!!** Data values associated to input meshes **!!
  for i=1:n
    v1D(i) = evalFun1D(fun, x(i), dxn);
  end

  %!!** True solution **!!
  for i=1:m
    v1Dout_true(i) = evalFun1D(fun, xout(i), dxn);
  end
 
  
  [v1Dout, degOut(1:n-1,1)] = adaptiveInterpolation1D(x, v1D, xout, d, 1, degOut(1:n-1, 1)); 

  %!!** open file **!!
  fname = openFile(fun, 1, nn, d);
  fid = fopen(char(fname), 'w');
  %!!** write to file **!!
  for i=1:m
    fprintf(fid,'%.8E \t %.8E \t %.8E \n', xout(i), v1Dout_true(i), v1Dout(i) );
  end
  %!!** close file **!!
  fclose(fid);

  [v1Dout, degOut(1:n-1,2)] = adaptiveInterpolation1D(x, v1D, xout, d, 2, degOut(1:n-1, 2)); 

  %!!** open file **!!
  fname = openFile(fun, 2, nn, d);
  fid = fopen(char(fname), 'w');
  %!!** write to file **!!
  for i=1:m
    fprintf(fid,'%.8E \t %.8E \t %.8E \n', xout(i), v1Dout_true(i), v1Dout(i) );
  end
  %!!** close file **!!
  fclose(fid);

  %!!** open file **!!
  fname = openFile(fun, 4, nn, d);
  fid = fopen(char(fname), 'w');
  %!!** write to file **!!
  for i=1:n-1
    fprintf(fid, '%d \t %d \n', degOut(i, 1), degOut(i, 2));
  end
  %!!** close file **!!
  fclose(fid);


  end %function 


function test001(d, fun, n, a, b, m)
%!! 
%!! Test interpolation with 1D runge function
%!! 
 
  %!!** Initialize variables **!!
  degOut = zeros(n-1, 2);
  degOut_lgl = zeros(n-1, 2);
  x = zeros(n,1);
  x_lgl = zeros(n, 1);
  v1D = zeros(n, 1);
  v1D_lgl = zeros(n,1);
  v1Dout = zeros(m,1);
  v1Dout_lgl = zeros(m,1);
  v1Dout_true = zeros(m,1);

  kk = 100.0;

  %!!** calculates intreval sizes **!!
  dxn = (b-a) /(n-1);
  dxm = (b-a) /(m-1);
  

  %!!** uniform mesh **!!
  for i=1:n
    x(i) = a + (i-1)*dxn;
  end

  dd = 8;
  ne = (n-1) / dd;                        %% calculates the number of elements

  %!!** LGL mesh **!!
  x_tmp = lglnodes(dd);
  x_tmp = flip(x_tmp); 
  dxn = (b-a) /ne;              %!! calculates element size
  xl = a;                                        %!! initialaziation 
  xr = a;                                       %!! initialization 
  is = 1;
  ie = 1;
  for i=1:ne
    xl = xr;                                     %!! left boundary of element i
    xr = xl + dxn;                               %!! right boun dary of element i
    is = ie;
    ie = is + dd;
    x_lgl(is:ie) = 0.50*( x_tmp* (xr-xl) + (xr+xl) );
  end 
  %!!** output mesh points **!
  for i=1:m
    xout(i) = a + (i-1)*dxm;
  end

  %!!** Data values associated to input meshes **!!
  if(fun == 5) 
    dxn = (b-a)/4.0 ;%!! indicate spacing between spike for fun=4
  end
  for i=1:n
    v1D(i) = evalFun1D(fun, x(i), dxn);
    v1D_lgl(i) =  evalFun1D(fun, x_lgl(i), dxn);
  end

  %!!** True solution **!!
  for i=1:m
    v1Dout_true(i) = evalFun1D(fun, xout(i), dxn);
  end
 

  %!!** interpolation using PCHIP **!!
  v1Dout = pchip(x, v1D, xout);
  v1Dout_lgl = pchip(x_lgl, v1D_lgl, xout);

  %!!** open file **!!
  fname=openFile(fun, 3, n, d);
  fid = fopen(char(fname), 'w');
  %!!** write to file **!!
  for i=1:m
    fprintf(fid,'%.8E \t %.8E \t %.8E \t %.8E \n', xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i) );
  end
  %!!** close file **!!
  fclose(fid);


  [v1Dout, degOut(1:n-1,1)] = adaptiveInterpolation1D(x, v1D, xout, d, 1); 
  [v1Dout_lgl, degOut_lgl(1:n-1,1)] = adaptiveInterpolation1D(x_lgl, v1D_lgl, xout, d, 1); 

  %!!** open file **!!
  fname = openFile(fun, 1, n, d);
  fid = fopen(char(fname), 'w');
  %!!** write to file **!!
  for i=1:m
    fprintf(fid,'%.8E \t %.8E \t %.8E \t %.8E \n', xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i) );
  end
  %!!** close file **!!
  fclose(fid);

  [v1Dout, degOut(1:n-1,2)] = adaptiveInterpolation1D(x, v1D, xout, d, 2); 
  [v1Dou_lgl, degOut_lgl(1:n-1,2)] = adaptiveInterpolation1D(x_lgl, v1D_lgl, xout, d, 2); 

  %!!** open file **!!
  fname= openFile(fun, 2, n, d);
  fid = fopen(char(fname), 'w');
  %!!** write to file **!!
  for i=1:m
    fprintf(fid,'%.8E \t %.8E \t %.8E \t %.8E \n', xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i) );
  end
  %!!** close file **!!
  fclose(fid);

  %!!** open file **!!
  fname = openFile(fun, 4, n, d);
  fid = fopen(char(fname), 'w');
  %!!** write to file **!!
  for i=1:n-1
    fprintf(fid,'%d \t %d \t %d \t %d \n', degOut(i, 1), degOut(i, 2), degOut_lgl(i, 1), degOut_lgl(i, 2) );
  end
  %!!** close file **!!
  fclose(fid);

end  %fucntion

function paperExamples2D()
%
%
%
  m = 1000;
  d = [1, 4, 8, 16];
  nx = [17, 33, 65, 129, 257];
  ny = [17, 33, 65, 129, 257];
  nx2 =[13, 25, 49, 97, 193];
  ny2 =[13, 25, 49, 97, 193];

  %!!** set up interval x \in [ax(i), bx(i)] and y \in [ay(i), by(i)]**!! 
  ax = [-1.0, -1.0, 0.0, -0.2 ];
  bx = [ 1.0,  1.0, 2.0, 0.2  ];
  ay = [-1.0, -1.0, 0.0, -0.2 ];
  by = [ 1.0,  1.0, 1.0, 0.2  ];
  
  %!!** function type 1=runge funtion , 2= heaviside, 3=Gelb Tanner **!! 
  fun = [1, 2, 3, 4]

  %!!** Example with truncated peak **!!
  for j=1:4                                                              
    fprintf('*****  d= %d ***** \n', d(j) );
    for i=1:5  
       %!!** Perform interpolation and calculate different L2-error
       %!    norms different methods **!!
       fprintf('nx=%d ny=%d \n', nx(i), ny(i) );
       fprintf('ax=%d bx=%d \n', ax(1), bx(1) );
       fprintf('ay=%d by=%d \n', ay(1), by(1) );
       peak2D(d(j), fun(1), nx(i)-1, ny(i)-1, ax(1), bx(1), ay(1), by(1), m);
    end
  end


  %!!** comparing against PCHIP **!!
  for k=1:4 
    if(k == 1 || k==4) 
    fprintf('*****  fun= %d *****', fun(k))
      for i=1:5 
         %** Perform interpolation and calculate different L2-error
         %   norms different methods **!!
         fprintf('nx=%d ny=%d \n', nx(i), ny(i) );
         fprintf('ax=%d bx=%d \n', ax(k), bx(k) ); 
         fprintf('ay=%d by=%d \n', ay(k), by(k) );
         test002(3, fun(k), nx(i), ny(i), ax(k), bx(k), ay(k), by(k), m)
      end
    end
  end


  for k=1:4 
    if(k==1 || k==4) 
    fprintf('*****  fun= %d***** \n',  fun(k));
    for j=1:4                
      fprintf('*****  d=%d ***** \n', d(j) );
      for i=1:5    
         %!!** Perform interpolation and calculate different L2-error
         %!    norms different methods **!!
         fprintf('nx= %d ny= %d', nx(i), ny(i) );
         fprintf('ax= %d bx= %d', ax(k), bx(k) );
         fprintf('ay= %d by= %d', ay(k), by(k) );
         test002(d(j), fun(k), nx(i), ny(i), ax(k), bx(k), ay(k), by(k), m);
      end
    end
    end
  end 

  fprintf( 'DONE \n');

end % function

function peak2D(d, fun, nx, ny, ax, bx, ay, by, m)
%!!
%!!
%!!

  %use mod_legendre
  %use mod_adaptiveInterpolation


  %implicit none


  %integer, intent(in)           :: fun                  !! function type 
  %integer, intent(in)           :: nx                    !! number of input points
  %integer, intent(in)           :: ny                    !! number of input points
  %integer, intent(in)           :: m                    !! number of output points
  %integer, intent(in)           :: d                    !! target interpolant degree
  %real(kind=8), intent(in)      :: ax, bx                 !! interval [a, b]
  %real(kind=8), intent(in)      :: ay, by                 !! interval [a, b]

  %integer                       :: limiter             !!
  %integer 		        :: i, j, k, fid, ierr, tmp_idx
  %integer 		        ::  nwk, seed
  %integer 		        :: is, ie, nnx, nny
  %real(kind=8)			:: x(nx)                 !! input mesh points  
  %real(kind=8)			:: y(ny)                 !! input mesh points  
  %integer 			:: degx(nx-1, ny)        !! input mesh points  
  %integer 			:: degx2(nx-1, ny)       !! input mesh points  
  %integer 			:: degy(ny-1, m)         !! input mesh points  
  %integer 			:: degy2(ny-1, m)        !! input mesh points  
  %real(kind=8)			:: v2D(nx, ny)           !! input data values
  %real(kind=8)			:: xout(m)               !! output points to be approximated 
  %real(kind=8)			:: yout(m)               !! output points to be approximated 
  %real(kind=8)			:: v2Dout(m, m)          !! approximated output values
  %real(kind=8)			:: v2Dout_true(m, m)       !! True values at output points
  %real(kind=8)			:: v2D_tmp(m, ny)        !! True values at output points
  %real(kind=8)			:: dxn, dxm, dyn, dym, err_L2, eps, start_t, end_t
  %real(kind=8)			:: h                    !! element spacing


  %character*12                  :: filename


  nnx = nx;
  nny = ny;
  
 

  %!!** Initialize variables **!!
  x = zeros(nx,1);
  y = zeros(ny,1);
  xout = zeros(m,1);
  yout = zeros(m,1);
  v2D = zeros(nx, ny);
  v2Dout = zeros(m,m);
  v2Dout_true = zeros(m,m);


  %!!** calculates intreval sizes **!!
  dxn = (bx-ax) /(nx-1);
  dxm = (bx-ax) /(m-1 );
  dyn = (by-ay) /(ny-1);
  dym = (by-ay) /(m-1);
  

  %!!** uniform mesh **!!
  for i=1:nx
    x(i) = ax + (i-1)*dxn;
  end
  for i=1:ny
    y(i) = ay + (i-1)*dyn;
  end

  %!!** output mesh points **!
  for i=1:m
    xout(i) = ax + (i-1)*dxm;
    yout(i) = ay + (i-1)*dym;
  end


  %!!** only used in calculation inside of evalFun2D for fun == 4
  h = dxn;

  %!!** Data values associated to input meshes **!!
  for j=1:ny
    for i=1:nx
      v2D(i,j) = evalFun2D(fun, x(i), y(j), h);
    end
  end

  %!!** True solution **!!
  for j=1:m
    for i=1:m
      v2Dout_true(i,j) = evalFun2D(fun, xout(i), yout(j), h);
    end
  end

  %!!**  Interpolation using Tensor product and DBI **!!
  for j=1:ny
    [v2D_tmp(:,j), degx(:,j)] =adaptiveInterpolation1D(x, v2D(:,j), xout, d, 1); 
  end
  for i=1:m
    [v2Dout(i,:), degy(:,i)] = adaptiveInterpolation1D(y, v2D_tmp(i,:), yout, d, 1); 
  end
 
  %!!** Open file **!! 
  fname = openFile2D(fun, 1, nnx, nny, d, 0);
  fid = fopen(char(fname), 'w');
  %!!** Write to open file **!!
  for j=1:m
    for i=1:m
      fprintf(fid,'%.8E \t %.8E \t %.8E \t %.8E \n', xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j) );
    end
  end
  %!!** close file **!!
  fclose(fid);

  %!!**  Interpolation using Tensor product and PPI **!!
  for j=1:ny
    [v2D_tmp(:,j), degx2(:,j)] = adaptiveInterpolation1D(x, v2D(:,j), xout, d, 2); 
  end
  for i=1:m
    [v2Dout(i,:), degy2(:,i)]= adaptiveInterpolation1D(y, v2D_tmp(i,:), yout, d, 2);
  end
 
  %!!** Open file **!! 
  fname = openFile2D(fun, 2, nnx, nny, d, 0);
  fid = fopen(char(fname), 'w');
  %!!** Write to open file **!!
  for j=1:m
    for i=1:m
      fprintf(fid,'%.8E \t %.8E \t %.8E \t %.8E \n', xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j));
    end
  end
  %!!** close file **!!
  fclose(fid);

end % function

function  test002(d, fun, nx, ny, ax, bx, ay, by, m)
%
%
%

  %use mod_legendre
  %use mod_adaptiveInterpolation


  %implicit none


  %integer, intent(in)           :: fun                  !! function type 
  %integer, intent(in)           :: nx                    !! number of input points
  %integer, intent(in)           :: ny                    !! number of input points
  %integer, intent(in)           :: m                    !! number of output points
  %integer, intent(in)           :: d                    !! target interpolant degree
  %real(kind=8), intent(in)      :: ax, bx                 !! interval [a, b]
  %real(kind=8), intent(in)      :: ay, by                 !! interval [a, b]
  %!!real(kind=8), intent(out)     :: errOut(3)           !! L2-errors for PCHIP, DBI and PPI
  %!real(kind=8), intent(out)     :: degOut(n-1, 2)      !! degree used for each interval

  %integer                       :: limiter             !!
  %integer 		        :: i, j, k, fid, ierr, tmp_idx
  %integer 		        ::  nwk, seed
  %integer 		        :: is, ie, nnx, nny, dd
  %integer 		        :: nex, ney              !! number of elements in x and y directions respectively
  %real(kind=8)			:: x(nx)                 !! input mesh points  
  %real(kind=8)			:: x_lgl(nx)                 !! input mesh points  
  %real(kind=8)			:: y(ny)                 !! input mesh points  
  %real(kind=8)			:: y_lgl(ny)             !! input mesh points  
  %integer 			:: degx(nx-1, ny)        !! input mesh points  
  %integer 			:: degx_lgl(nx-1, ny)    !! input mesh points  
  %integer 			:: degx2(nx-1, ny)       !! input mesh points  
  %integer 			:: degx2_lgl(nx-1, ny)   !! input mesh points  
  %integer 			:: degy(ny-1, m)         !! input mesh points  
  %integer 			:: degy_lgl(ny-1, m)     !! input mesh points  
  %integer 			:: degy2(ny-1, m)        !! input mesh points  
  %integer 			:: degy2_lgl(ny-1, m)    !! input mesh points  
  %real(kind=8)			:: v2D(nx, ny)           !! input data values
  %real(kind=8)			:: v2D_lgl(nx, ny)          !! input data values
  %real(kind=8)			:: xout(m)               !! output points to be approximated 
  %real(kind=8)			:: yout(m)               !! output points to be approximated 
  %real(kind=8)			:: v2Dout(m, m)          !! approximated output values
  %real(kind=8)			:: v2Dout_lgl(m, m)         !! approximated output values
  %real(kind=8)			:: v2Dout_true(m, m)       !! True values at output points
  %real(kind=8)			:: v2D_tmp(m, ny)        !! True values at output points
  %!!real(kind=8)			:: dd(3)                 !!
  %real(kind=8)			:: dxn, dxm, dyn, dym, err_L2, eps, start_t, end_t
  %real(kind=8)			:: x_tmp(9), w_tmp(9)
  %real(kind=8)			:: xl, xr, yl, yr
  %real(kind=8)			:: h                    !! element spacing


  %!!real(kind=8)			:: kk, tmp
  %real(kind=8)			:: wk((nx+1)*2), d_tmp(nx+1)
  %real(kind=8)			:: wk2((ny+1)*2), d_tmp2(ny+1)
  %real(kind=8)			:: fdl(m)
  %character*12                  :: filename
  %logical                       :: spline

 
 

  %!!** Initialize variables **!!
  x = zeros(nx,1);
  y = zeros(ny,1);
  xout = zeros(m,1);
  yout = zeros(m,1);
  x_lgl = zeros(nx,1);
  y_lgl = zeros(ny,1);
  v2D = zeros(nx, ny);
  v2D_lgl = zeros(nx, ny);
  v2Dout = zeros(m,m);
  v2Dout_lgl = zeros(m,m);
  v2Dout_true = zeros(m,m);
  v2D_tmp = zeros(m, ny);
  degx = zeros(nx-1, ny);
  degy = zeros(ny-1, m);
  degx_lgl = zeros(nx-1, ny);
  degy_lgl = zeros(ny-1, m);
  degx2 = zeros(nx-1, ny);
  degy2 = zeros(ny-1, m);
  degx2_lgl = zeros(nx-1, ny);
  degy2_lgl = zeros(ny-1, m);


  %!!** calculates intreval sizes **!!
  dxn = (bx-ax) /(nx-1);
  dxm = (bx-ax) /(m-1 );
  dyn = (by-ay) /(ny-1);
  dym = (by-ay) /(m-1 );
  

  %!!** uniform mesh **!!
  for i=1:nx
    x(i) = ax + (i-1)*dxn;
  end
  for i=1:ny
    y(i) = ay + (i-1)*dyn;
  end

  %!!** number of elements **!!
  dd = 8
  nex = (nx-1) / dd;
  ney = (ny-1) / dd;

  %!!** Legendre gauss lobatto points **!!
  x_tmp = lglnodes(dd);
  x_tmp = flip(x_tmp);
  dxn = (bx-ax) /nex;
  xl = ax;
  xr = ax;  
  is = 1;
  ie = 1;
  for i=1:nex
    xl = xr; 
    xr = xl + dxn; 
    is = ie;
    ie = is + dd;
    %!!** maping from [-1,1] to [xl, xr] 
    x_lgl(is:ie) = x_tmp* (xr-xl)/2.0 + (xr+xl)/2.0;
  end
  dyn = (by-ay)/ney;              %% calculates element size
  yl = ay;
  yr = ay;
  is = 1;
  ie = 1;
  for i=1:ney
    yl = yr;                                     %% left boundary of element i
    yr = yl + dyn;                                %% right boun dary of element i
    is = ie;
    ie = is + dd;
    %!!** maping from [-1,1] to [yl, yr] 
    y_lgl(is:ie) = x_tmp* (yr-yl)/2.0 + (yr+yl)/2.0;
  end

  %!!** output mesh points **!
  for i=1:m
    xout(i) = ax + (i-1)*dxm;
    yout(i) = ay + (i-1)*dym;
  end


  %!!** only used in calculation inside of evalFun2D for fun == 4
  h = (bx-ax)/((nx-1)/d);

  %!!** Data values associated to input meshes **!!
  for j=1:ny
    for i=1:nx
      v2D(i,j) = evalFun2D(fun, x(i), y(j), h);
      v2D_lgl(i,j) = evalFun2D(fun, x_lgl(i), y_lgl(j), h);
    end
  end

  %!!** True solution **!!
  for j=1:m
    for i=1:m
      v2Dout_true(i,j) = evalFun2D(fun, xout(i), yout(j), h);
    end
  end
 
  %!!**  Interpolation using Tensor product and PCHIP **!!
  for j=1:ny
    v2D_tmp(:,j) = pchip(x, v2D(:,j), xout);
  end
  for i=1:m
    v2Dout(i,:) = pchip(y, v2D_tmp(i,:), yout);
  end
  %!!** interpolation using lgl **!!
  for j=1:ny
    v2D_tmp(:,j) = pchip(x_lgl, v2D_lgl(:,j), xout);
  end
  for i=1:m
    v2Dout_lgl(i,:) = pchip(y_lgl, v2D_tmp(i,:), yout);
  end


  %!!** Open file **!! 
  fname = openFile2D(fun, 3, nx, ny, d, 0);
  fid = fopen(char(fname), 'w');
  %!!** Write to open file **!!
  for j=1:m
    for i=1:m
      fprintf(fid,'%.8E \t %.8E \t %.8E \t %.8E \t %.8E \n ', ...
          xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j), v2Dout_lgl(i, j) );
    end
  end
  %!!** close file **!!
  fclose(fid);

  %!!**  Interpolation using Tensor product and DBI **!!
  for j=1:ny
    [v2D_tmp(:,j), degx(:,j)] = adaptiveInterpolation1D(x, v2D(:,j), xout, d, 1); 
  end
  for i=1:m
    [v2Dout(i,:), degy(:,i)] = adaptiveInterpolation1D(y, v2D_tmp(i,:), yout, d, 1); 
  end

  %!!**  Interpolation using Tensor product and DBI **!!
  for j=1:ny
    [v2D_tmp(:,j), degx_lgl(:,j)] = adaptiveInterpolation1D(x_lgl, v2D_lgl(:,j), xout, d, 1); 
  end
  for i=1:m
    [v2Dout_lgl(:,i), degy_lgl(:,i)] = adaptiveInterpolation1D(y_lgl, v2D_tmp(i,:), yout, d, 1); 
  end
  
  %!!** Open file **!! 
  fname= openFile2D(fun, 1, nx, ny, d, 0);
  fid = fopen(char(fname), 'w');
  %!!** Write to open file **!!
  for j=1:m
    for i=1:m
      fprintf(fid,'%.8E \t %.8E \t %.8E \t %.8E \t %.8E \n ', ...
               xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j), v2Dout_lgl(i, j) );
    end
  end
  %!!** close file **!!
  fclose(fid)

  %!!**  Interpolation using Tensor product and PPI **!!
  for j=1:ny
    [v2D_tmp(:,j), degx2(:,j)] = adaptiveInterpolation1D(x, v2D(:,j), xout, d, 2); 
  end
  for i=1:m
    [v2Dout(i,:), degy2(:,i)] = adaptiveInterpolation1D(y, v2D_tmp(i,:), yout, d, 2); 
  end

  %!!**  Interpolation using Tensor product and DBI **!!
  for j=1:ny
    [v2D_tmp(:,j), degx2_lgl(:,j)] = adaptiveInterpolation1D(x_lgl, v2D_lgl(:,j), xout, d, 2); 
  end
  for i=1:m
    [v2Dout_lgl(:,i), degy2_lgl(:,i)] = adaptiveInterpolation1D(y_lgl, v2D_tmp(i,:), yout, d, 2); 
  end
 

  %!!** Open file **!! 
  fname = openFile2D(fun, 2, nx, ny, d, 0);
  fid = fopen(char(fname), 'w');
  %!!** Write to open file **!!
  for j=1:m
    for i=1:m
      fprintf(fid,'%.8E \t %.8E \t %.8E \t %.8E \t %.8E \n ', ...
      xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j), v2Dout_lgl(i, j) );
    end
  end
  %!!** close file **!!
  fclose(fid)

  %!!** Open file **!! 
  fname = openFile2D(fun, 4, nx, ny, d, 1);
  fid = fopen(char(fname), 'w');
  %!!** Write to open file **!!
  for j=1:ny
    for i=1: nx-1
      fprintf(fid,'%d \t %d \t %d \t %d \n ', ...
      degx(i,j), degx2(i,j), degx_lgl(i,j), degx2_lgl(i,j) );
    end
  end
  %!!** close file **!!
  fclose(fid)

  %!!** Open file **!! 
  fname  = openFile2D(fun, 4, nx, ny, d, 2)
  fid = fopen(char(fname), 'w');
  %!!** Write to open file **!!
  for i=1:m
    for j=1:ny-1
      fprintf(fid,'%d \t %d \t %d \t %d \n ', ...
      degy(j,i), degy2(j,i), degy_lgl(j,i), degy2_lgl(j,i) );
    end
  end
  %!!** close file **!!
  fclose(fid)

end % function 

function v = evalFun1D(fun, x, h)
%!!
%!! To evaluate different fuunction. fun determine
%!! the function that will be evaluated at the points x
%!!

 
  %!!** intialize variables **!!
  k = 100;
  pi = 4.0*atan(1.0); 
  t = x/h;
  delta = 0.01;
  %!!** a and b are only used for fun =5. When using fun=4 the 
  %!!   a and b must be the same as the ones Paperexample1D **!!
  a = -2.0;
  b = 0.0;
  ne = (b-a)/h;
  

  %!!** 1D runge function **!!
  if(fun == 1) 
    v = 1.0 / (1.0 + 25.0 * x * x);
  %!!** heaviside function **!!
  elseif(fun == 2)
    v = 1.0/(1.0 + exp(-2*k*x));
  %!!** Gelb and Tanner function **!!
  elseif(fun == 3)
    if(x < -0.5) 
      v = 1.0 + (2.0* exp(2.0 * pi *(x+1.0)) - 1.0 -exp(pi)) /(exp(pi)-1.0);
    else
      v = 1.0 - sin(2*pi*x / 3.0 + pi/3.0);
    end
  %!!** Modified tanh function **!!
  elseif(fun==4)
    k = 10.0;
    if(a <= x && x <= a+h)
      v = tanh(x*k);
    elseif(a+h<= x && x <= a+2*h)
      v = 2*tanh(x*k)         - tanh((a+h)*k);
    end
    for i=3:ne 
      if(a +(i-1)*h <= x && x <= a+i*h);
        v = i*tanh(x*k) - tanh((a+h)*k);
        for j=2:i-1 
         v = v - tanh( (a+(j*h))*k );
        end
      end
    end
    v = 1.0 + v;
  %!!** modified square function **!!
  elseif(fun == 5)
    v = 1.0 - abs( 2/ pi * atan( sin(pi*t) /delta));
  end


end % function evalFun1D

function v = evalFun2D(fun, x, y, h)
%
% To evaluate different fuunction. fun determine
% the function that will be evaluated at the points x
%
 
  %!!** intialize variables **!!
  k = 100;
  pi = 4.0*atan(1.0); 
  delta = 0.01;
  

  %!!** 1D runge function **!!
  if(fun==1)
    v = 1.0 / ( 1.0 + 25.0 * ( x*x + y*y) );
  elseif(fun==2)
    v = max(0.0, sin(4.0*atan(1.0)*x)*sin(4.0*atan(1.0)*y) ) 
  elseif(fun==3)then
    if( (x-1.5)*(x-1.5) + (y-0.5)*(y-0.5) <= 1.0/16.0 )
      v = 2.0*0.5*( cos(8.0*atan(1.0)*sqrt((x-1.5)*(x-1.5) + (y-0.5)*(y-0.5))));
    elseif(y-x >= 0.5)then
      v = 1.0;
    elseif(0.0 <= y-x && y-x <= 0.5)
      v = 2.0*(y-x);
    else
      v = 0.0;
    end
  elseif(fun==4)
    v = 1.0 / ( 1.0 + exp(-2*k*(x+y)*sqrt(2.0)/2.0) );
  end

end %function evalFun2D


%%
function fname = openFile2D(fun, limiter, nx, ny, d, did)
%
%
%
 
  if(d <10 )
    sfnumber =strcat("0", string(d), string(nx), string(ny));
  else
    sfnumber =strcat(string(d), string(nx), string(ny));
  end
  %!!** Runge function **!
  if(fun==1 && limiter==0) 
    fname =  strcat( "Runge2STD", sfnumber);
    %write(fname, '("Runge2STD", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif(fun==1 && limiter==1) 
    fname =  strcat( "Runge2DBI", sfnumber);
    %write(fname, '("Runge2DBI", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==1 && limiter==2)
    fname =  strcat( "Runge2DBI", sfnumber);
    %write(fname, '("Runge2PPI", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==1 && limiter==3) 
    fname =  strcat( "Runge2PCHIP", sfnumber);
    %write(fname, '("Runge2PCHIP", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==1 && limiter==4 && did==1)
    fname =  strcat( "Runge2Deg1", sfnumber);
    %write(fname, '("Runge2Deg1", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==1 && limiter== 4 && did==2)
    fname =  strcat( "Runge2Deg2", sfnumber);
    %write(fname, '("Runge2Deg2", i2.2, i3.3, "x", i3.3)')d, nx, ny
  %!!** Terrain 1 function **!
  elseif(fun==2 && limiter==0)
    fname =  strcat( "T1STD", sfnumber);
    %write(fname, '("T1STD", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif(fun==2 && limiter==1)
    fname =  strcat( "T1DBI", sfnumber);
    %write(fname, '("T1DBI", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==2 && limiter==2) 
    fname =  strcat( "T1PPI", sfnumber);
    %write(fname, '("T1PPI", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==2 && limiter==3)
    fname =  strcat("T1PPI", sfnumber);
    %write(fname, '("T1PCHIP", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==2 && limiter==4 && did==1) 
    fname =  strcat("T1Deg1", sfnumber);
    %write(fname, '("T1Deg1", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==2 && limiter==4 && did==2)
    fname =  strcat("T1Deg2", sfnumber);
    %write(fname, '("T1Deg2", i2.2, i3.3, "x", i3.3)')d, nx, ny
  %!!** Terrain 3 function **!
  elseif(fun==3 && limiter==0)
    fname =  strcat( "T2STD", sfnumber);
    %write(fname, '("T2STD", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif(fun==3 && limiter==1)
    fname =  strcat( "T2DBI", sfnumber);
    %write(fname, '("T2DBI", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==3 && limiter==2) 
    fname =  strcat( "T2PPI", sfnumber);
    %write(fname, '("T2PPI", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==3 && limiter==3) 
    fname =  strcat( "T2PCHIP", sfnumber);
    %write(fname, '("T2PCHIP", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==3 && limiter== 4 && did==1) 
    fname =  strcat( "T2Deg1", sfnumber);
    %write(fname, '("T2Deg1", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==3 && limiter==4 && did==2)
    fname =  strcat( "T2Deg2", sfnumber);
    %write(fname, '("T2Deg2", i2.2, i3.3, "x", i3.3)')d, nx, ny
  %!!** modified 2D heavisde function  function **!
  elseif(fun==4 && limiter==0)
    fname =  strcat( "Heaviside2STD", sfnumber);
    %write(fname, '("Heaviside2DSTD", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif(fun==4 && limiter==1) 
    %write(fname, '("Heaviside2DDBI", i2.2, i3.3, "x", i3.3)')d, nx, ny
    fname =  strcat( "Heaviside2DBI", sfnumber);
  elseif (fun==4 && limiter==2) 
    fname =  strcat( "Heaviside2PPI", sfnumber);
    %write(fname, '("Heaviside2DPPI", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==4 && limiter==3)
    fname =  strcat( "Heaviside2PCHIP", sfnumber);
    %write(fname, '("Heaviside2DPCHIP", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==4 && limiter==4 && did==1)
    fname =  strcat( "Heaviside2Deg1", sfnumber);
    %write(fname, '("Heaviside2DDeg1", i2.2, i3.3, "x", i3.3)')d, nx, ny
  elseif (fun==4 && limiter==4 && did==2)
    fname =  strcat( "Heaviside2Deg2", sfnumber);
    %write(fname, '("Heaviside2DDeg2", i2.2, i3.3, "x", i3.3)')d, nx, ny
  end

  %open(unit=fid,file=fname, status='unknown')
end %subroutine


function fname = openFile(fun, limiter, n, d)
%
%
%
%

  %!!** calculates file number the first 2 digit represents 
  %!    d, the interpolant degree used. The remaining digits
  %!    represents the n the number of points used **!!
  if (n >= 100)
    fnumber = d*1000 + n;
  else
    fnumber = d*1000 + n;
  end
  if(fnumber < 10000)
    sfnumber =strcat("0", string(fnumber));
  else
    sfnumber = string(fnumber);
  end
 
  %!!** Runge function **!
  if(fun==1 && limiter== 0) 
    fname =  strcat( "RungeSTD", sfnumber);
    %write(fname, '("RungeSTD", i5.5)')fnumber
  elseif(fun==1 && limiter == 1) 
    fname =  strcat( "RungeDBI", sfnumber);
    %write(fname, '("RungeDBI", i5.5)')fnumber
  elseif (fun== 1 && limiter==2) 
    fname =  strcat( "RungePPI", sfnumber);
    %write(fname, '("RungePPI", i5.5)')fnumber
  elseif (fun==1 && limiter==3)
    fname =  strcat( "RungePCHIP", sfnumber);
    %write(fname, '("RungePCHIP", i5.5)')fnumber
  elseif (fun==1 && limiter==4) 
    fname =  strcat( "RungeDeg", sfnumber);
    %write(fname, '("RungeDeg", i5.5)')fnumber
  %!!** Heaviside function !!
  elseif(fun==2 && limiter==0)
    fname =  strcat( "HeavisideSTD", sfnumber);
    %write(fname, '("HeavisideSTD", i5.5)')fnumber
  elseif(fun==2 && limiter==1) 
    fname =  strcat( "HeavisideDBI", sfnumber);
    %write(fname, '("HeavisideDBI", i5.5)')fnumber
  elseif (fun==2 && limiter==2)
    fname =  strcat( "HeavisidePPI", sfnumber);
    %write(fname, '("HeavisidePPI", i5.5)')fnumber
  elseif (fun==2 && limiter==3)
    fname =  strcat( "HeavisidePCHIP", sfnumber);
    %write(fname, '("HeavisidePCHIP", i5.5)')fnumber
  elseif (fun==2 && limiter==4) 
    fname =  strcat( "HeavisideDeg", sfnumber);
    %write(fname, '("HeavisideDeg", i5.5)')fnumber
  %!!** Heaviside function !!
  elseif(fun==3 && limiter==0) 
    fname =  strcat( "GelbTSTD", sfnumber);
    %write(fname, '("GelbTSTD", i5.5)')fnumber
  elseif(fun==3 && limiter==1) 
    fname =  strcat( "GelbTDBI", sfnumber);
    %write(fname, '("GelbTDBI", i5.5)')fnumber
  elseif (fun==3 && limiter==2) 
    fname =  strcat( "GelbTPPI", sfnumber);
    %write(fname, '("GelbTPPI", i5.5)')fnumber
  elseif (fun==3 && limiter==3)
    fname =  strcat( "GelbTPCHIP", sfnumber);
    %write(fname, '("GelbTPCHIP", i5.5)')fnumber
  elseif (fun==3 && limiter==4) then
    fname =  strcat( "GelbTDeg", sfnumber);
    %write(fname, '("GelbTDeg", i5.5)')fnumber
  %!!** Modified tanh  function !!
  elseif(fun==4 && limiter==0)
    fname =  strcat( "TanhSTD", sfnumber);
    %write(fname, '("TanhSTD", i5.5)')fnumber
  elseif(fun==4 && limiter==1)
    fname =  strcat( "TanhSTD", sfnumber);
    %write(fname, '("TanhDBI", i5.5)')fnumber
  elseif (fun==4 && limiter==2)
    fname =  strcat( "TanhPPI", sfnumber);
    %write(fname, '("TanhPPI", i5.5)')fnumber
  elseif (fun==4 && limiter== 3) 
    fname =  strcat( "TanhPCHIP", sfnumber);
    %write(fname, '("TanhPCHIP", i5.5)')fnumber
  elseif (fun==4 && limiter==4)
    fname =  strcat( "TanhDeg", sfnumber);
    %write(fname, '("TanhDeg", i5.5)')fnumber
  %!!** Modified aquare  function !!
  elseif(fun==5 && limiter==0)
    fname =  strcat( "SquareSTD", sfnumber);
    %write(fname, '("SquareSTD", i5.5)')fnumber
  elseif(fun==5 && limiter==1) 
    fname =  strcat( "SquareSTD", sfnumber);
    %write(fname, '("SquareDBI", i5.5)')fnumber
  elseif (fun==5 && limiter==2) 
    fname =  strcat( "SquarePPI", sfnumber);
    %write(fname, '("SquarePPI", i5.5)')fnumber
  elseif (fun==5 && limiter==3) 
    fname =  strcat( "SquarePCHIP", sfnumber);
    %write(fname, '("SquarePCHIP", i5.5)')fnumber
  elseif (fun==5 && limiter==4) 
    fname =  strcat( "SquareDeg", sfnumber);
    %write(fname, '("SquareDeg", i5.5)')fnumber
  end

  %fid = fopen(char(fname), 'w');
end % function



%
%---------------------------------------------------------------------------------------------%
% Driver to produce the  1D and 2D approximation presented in the manuscript
%---------------------------------------------------------------------------------------------%
%

  fileID = fopen('approximations_tables_1d_2d.txt', 'w');

  %% 1D function approximations %
  %approximations1D();
  %movefile Runge* mapping_data/data
  %movefile Heavi* mapping_data/data
  %movefile GelbT* mapping_data/data
  
  %plot_approximations ;
  %fprintf('Tables for 1D approximations saved in approximations_tables_1d_2d.txt \n');
 
  
  fprintf('Press any key to continue to the mapping examples. \n')
  pause
  %% mapping examples %
  %for k= [64, 127, 253];
  %  mapping(k)
  %end
  %movefile runge* mapping_data/data
  %movefile qc* mapping_data/data
  %
  %fprintf('The approximated solutions are save in mapping_data/data. \n')
  %fprintf('running plot_approximations.m and plot_mapping.m to  \n')
  %fprintf('produce the figures and tables in the manuscript \n')
   
  %plot_mapping ;
  
  %plot_approximations ;
  fprintf('press any key to continue to the 2D function approximation examples.\n')
  frpintf(' WARNING: the 2D examples take long time because the solution is evaluated')
  fprintf('onto a 1000 x 1000 mesh for each example and saved. \n')

  % 2D function approximations %
  approximations2D()
  movefile Runge* mapping_data/data
  movefile T1* mapping_data/data
  movefile T2* mapping_data/data
  movefile Heavi* mapping_data/data
 
  plot_approximations2D ;
  fprintf('Tables for 1D and 2D approximations saved in approximations_tables_1d_2d.txt \n');
  fclose(fileID);
  % end of script


%---------------------------------------------------------------------------------------------%
% Functions used for the  1D approximation examples in the manuscript
%---------------------------------------------------------------------------------------------%

function approximations1D()
%
% approximation1D is used to set up the 
% diffferent configuration used to produces
% the approximation results for the 1D functions
% presented in the manuscript.
%

  fprintf('Running 1D approximation examples ...  ... \n')
  %%** Initialization **%%
  m = 10000;
  n = [17, 33, 65, 129, 257];                                              
  d = [1, 4, 8, 16];                                                       
  a = [-1.0, -0.2, -1.0];
  b = [ 1.0,  0.2,  1.0];
  eps_test = [ 1.0,  0.1,  0.01, 0.001, 0.0001, 0.00 ];
  sten = [1, 2, 3];

  %%** modify eps0 and eps1 to change the bounds on the interpolant **%%
  eps0 = 0.01;
  eps1 = 1.0;

  %%** functions 1=Runge , 2= heaviside, 3=Gelb Tanner **%% 
  fun = [1, 2, 3];        

  %%** Used to evaluate different choices of eps0 **%%
  testepsilon1D(sten(2), eps_test, eps1, d(3), n(1), a, b,  m);

  for ii=1:3
    fprintf('**** parameter st= %d ****\n', sten(ii) );
    for k=1:3

      %%** Third order resulst using DBI, PPI, and PCHIP **%%
      fprintf('*****  function = %d ******** \n', fun(k) ) ;
      for i=1:5                                                             
        test001(3, eps0, eps1, sten(ii), fun(k), n(i), a(k), b(k), m, 8);
      end

      %%** higher degree interpolation methods using DBI and PPI **%%
      for j=1:4                                                             
        fprintf('*****  d= %d ***** \n', d(j) );
        for i=1:5                                                            
           test001(d(j), eps0, eps1, sten(ii), fun(k), n(i), a(k), b(k), m, 8);
        end  
      end 
    end 
  end

end 

function testepsilon1D(sten, eps0, eps1, d, n, a, b,  m)
%
% testepsilon1D aprroximates the Runge, smoothed Heaviside, and
% Gelb and Tanner functions with different values of eps0 that
% are used to bound the interpolant in the case of the PPI method. 
% This function produces the results used to build the 1D figures 
% In the manuscript.
%
% INPUT
% sten: stencil selction procedure (sten=1, sten=2, sten=3) 
% eps0(6): array of values of eps0 
% d:  traget polynomial degree for each interpolant
% n: number of points
% a(3): left boundaries
% b(3): right boundaries
% m: number of output points 
%


  %%** Initialize parameters **%%
  degOut = zeros(n-1,7);
  v1D = zeros(n,1);
  v1Dout = zeros(n,9);
  for k=1:3
    
    %%** calculates intreval sizes **%%
    dxn = (b(k)-a(k)) /double(n-1);
    dxm = (b(k)-a(k)) /double(m-1);
  

    %%** uniform mesh **%%
    for i=1:n
      x(i) = a(k) + double(i-1)*dxn;
    end
    x(n) = b(k);

    %%** output mesh points **!
    for i=1:m
      v1Dout(i, 1) = a(k) + double(i-1)*dxm;
    end
    v1Dout(m,1) = b(k);

    %%** Data values associated to input meshes **%%
    for i=1:n
      v1D(i) = evalFun1D(k, x(i), dxn);
    end

    %%** True solution **%%
    for i=1:m
      v1Dout(i,2) = evalFun1D(k, v1Dout(i, 1), dxn);
    end

    for i=1:7
      if(i==7)
        tmp = adaptiveInterpolation1D(x, v1D, v1Dout(:,1), d, 1, sten );
        v1Dout(:,2+i) = tmp; 
      else
        tmp = adaptiveInterpolation1D(x, v1D, v1Dout(:,1), d, 2, sten, eps0(i), eps1 ); 
        v1Dout(:,2+i) = tmp; 
      end
    end


    %%** open file **%%
    if( k == 1)
      fid = fopen('RungeEps', 'w');
    elseif( k == 2)
      fid = fopen('HeavisideEps', 'w');
    elseif( k == 3)
      fid = fopen('GelbTEps', 'w');
    end
    %%** write to file **%%
    for i=1:m
      for j=1:9
        fprintf(fid, '%.8E \t', v1Dout(i,j) );
      end
      fprintf(fid, '\n');
    end
    %%** close file **%%
    fclose(fid);
  end
end 

function test001(d, eps0, eps1, sten, fun, n, a, b, m, d_el)
% 
% test001 is used to approximate the Runge, smoothed Heaviside
% and Gelbd and Tanner function using different interpolation 
% methods. This function is used toproduce the 1D results presented
% in the manuscript.
%
% INPUT
% d: maximum polynomial degree for each interval
% eps0: positive user-supplied value used to bound interpolant for 
%       intervalswith no extrema.
% eps1: positive user-supplied value used to bound interpolant for 
%       intervals with extrema.
% sten: user-supplied value used to indicate stencil selection process
%       possible choices are sten=1, sten=2, sten=3.
% fun: used to indicate function used
% n: number of input points
% m: number of output points
% a: global interval left boundary
% b: global right interval boundary
% d_el: number of LGL points in each element
% 

  %%** To be used to identify file uniquely **%%
  if (d < 10)
    fnumber =  strcat("0", string(d*1000+n));
  else
    fnumber =  string(d*1000+n);
  end

  %%** get function  name **%%
  if(fun ==1)
    fun_name = "Runge";
  elseif(fun ==2)
    fun_name = "Heaviside";
  elseif(fun ==3)
    fun_name = "GelbT";
  else
    fprintf('ERROR: Invalid fun = %d \n', fun);
    fprintf( 'Invalid function value the possible values are fun =1, 2, or 3 \n');
    stop
  end

  %%** get stencil selection procedure. Needed to create file name where
  %    where the results will saved  **%%
  if(sten ==1) 
    sst = "st=1";
  elseif(sten ==2) 
    sst = "st=2";
  elseif(sten ==3) 
    sst = "st=3";
  else
    fprintf('ERROR: Invalid paparamter sten = %d \n', fun);
    fprintf('ERROR: Invalid paparamter st. The possible options are st=1, 2, or 3 \n');
    stop
  end
  
  %%** uniform mesh **%%
  dxn = (b-a) /double(n-1);
  for i=1:n
    x(i) = a + double(i-1)*dxn;
  end
  x(n)=b;

  dd = d_el;    %% number of LGL points in each element
  ne = (n-1) / dd;    %% calculates the number of elements
  
  %%** LGL mesh **%%
  x_tmp = lglnodes(dd);    %% LGL nodes in decreasing order
  x_tmp = flip(x_tmp);    %% LGL nodes in increasing order
  dxn = (b-a) / double(ne);    %% calculates element size
  xl = a;                                                   %% initialaze element left boundary 
  xr = a;                                                   %% initialize element right boundary 
  is = 1;
  ie = 1;
  for i=1:ne
    xl = xr;                                      %% update left boundary of element i
    xr = xl + dxn;                                %% update right boun dary of element i
    is = ie;
    ie = is + dd;
    x_lgl(is:ie) = 0.50*( x_tmp* (xr-xl) + (xr+xl) );
  end
  x_lgl(n)=b;

  %%** output mesh points **!
  dxm = (b-a) /double(m-1);
  for i=1:m
    xout(i) = a + double(i-1)*dxm;
  end
  xout(m) = b;

  %%** Data values associated to input meshes **%%
  dxn = (b-a)/double(ne); %% dummy variable not used for calculations
  for i=1:n
    v1D(i) = evalFun1D(fun, x(i), dxn);
    v1D_lgl(i) = evalFun1D(fun, x_lgl(i), dxn);
  end

  %%** True solution **%%
  for i=1:m
    v1Dout_true(i) = evalFun1D(fun, xout(i), dxn);
  end

  %%** interpolation using PCHIP **%%
  if(d ==3 ) 
    v1Dout = pchip(x, v1D, xout);
    v1Dout_lgl = pchip(x_lgl, v1D_lgl, xout);


    %%** open file and write to file **%%
    fname = strcat( fun_name, "PCHIP", fnumber);
    fid = fopen(char(fname), 'w');
    for i=1:m
      fprintf(fid,' %.8E \t %.8E \t %.8E \t %.8E \n ', xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i));
    end
    fclose(fid);
  end

  %%** Interpolation using DBI **%%
  tmp = adaptiveInterpolation1D(x, v1D, xout, d, 1, sten); 
  v1Dout = tmp;
  tmp = adaptiveInterpolation1D(x_lgl, v1D_lgl, xout, d, 1, sten); 
  v1Dout_lgl = tmp;

  %%** open file and write to file **%%
  fname = strcat( fun_name, "DBI", fnumber, sst);
  fid = fopen(char(fname), 'w');
  for i=1:m
    fprintf(fid,' %.8E \t %.8E \t %.8E \t %.8E \n ', xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i));
  end
  fclose(fid);

  %%** Interpolation using PPI **%%
  tmp = adaptiveInterpolation1D(x, v1D, xout, d, 2, sten, eps0, eps1); 
  v1Dout = tmp;
  tmp = adaptiveInterpolation1D(x_lgl, v1D_lgl, xout, d, 2, sten, eps0, eps1 ); 
  v1Dout_lgl = tmp;

  %%** open file and write to file **%%
  fname = strcat( fun_name, "PPI", fnumber, sst);
  fid = fopen(char(fname), 'w');
  for i=1:m
    fprintf(fid,' %.8E \t %.8E \t %.8E \t %.8E \n ', xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i));
  end
  fclose(fid);

end 


%---------------------------------------------------------------------------------------------%
% Functions used for the  2D approximation examples in the manuscript
%---------------------------------------------------------------------------------------------%
function approximations2D()
%
% approximation2D is used to set up the 
% diffferent configurations used to produces
% the approximation results for the 2D functions
% presented in the manuscript.
%

  fprintf('Running 2D approximation examples ...  ... \n')

  d = [1, 4, 8, 16];  %% array with interpolants degrees
  nx = [17, 33, 65, 129, 257];  %% array with number of inputpoints
  ny = [17, 33, 65, 129, 257];  %% array with number of inputpoints

  eps_test = [ 1.0,  0.1,  0.01, 0.001, 0.0001, 0.00 ];  % eps0 values to be used

  %%** set up interval x \in [ax(i), bx(i)] and y \in [ay(i), by(i)]**%% 
  ax = [-1.0, -1.0, 0.0, -0.2 ];  % left boundary for  x
  bx = [ 1.0,  1.0, 2.0, 0.2  ];  % right boundary for x
  ay = [-1.0, -1.0, 0.0, -0.2 ];  % left boundary ffor y
  by = [ 1.0,  1.0, 1.0, 0.2  ];  % right boundary for y
  
  m =  1000;  % number of output points in each direction

  %%** function type 1=runge funtion , 2= heaviside, 3=Gelb Tanner **%% 
  fun = [1, 2, 3, 4];                                                     %% function type

  %%**
  sten = [1, 2, 3];  % possible stencil choices for stencil construction
  eps0 = 0.01;  % user-supplied value used to bound interpolants in intervals with no extrema
  eps1 = 1.0;  % user-supplied value used to bound interpolant in tervals with extrema
  testepsilon2D(sten(2),eps0, eps1, d(3), nx(1), ny(1), ax, bx, ay, by, 100);
  for ii=1:3
    %%** comparing against PCHIP **%%
    for k=1:4 
      %%if(k == 1 .or. k == 4) 
        for i=1:5
           %%** Perform interpolation and calculate different L2-error
           %!    norms different methods **%%
           fprintf('nx=%d \t ny=%d \n', nx(i), ny(i) );
           fprintf('ax=%d \t bx=%d \n', ax(k), bx(k) );
           fprintf('ay=%d \t by=%d \n', ay(k), by(k) );
           test002(3, eps0, eps1, sten(ii), fun(k), nx(i), ny(i), ax(k), bx(k), ay(k), by(k), m, 8)
        end


        %% Higher degree interpolants %%
        fprintf('*****  fun=%d *****', fun(k) );
        for j=1:4  
          fprintf('*****  d= %d ***** \n', d(j) );
          for i=1:5
             %%** Perform interpolation and calculate different L2-error
             %%    norms different methods **%%
             fprintf('nx= %d \t ny=%d \n', nx(i), ny(i) );
             fprintf('ax(k)= %d \t bx(k) = %d \n', ax(k), bx(k) );
             fprintf('ay(k)= %d \t bx(k) = %d \n', ay(k), by(k) );
             test002(d(j), eps0, eps1, sten(ii), fun(k), nx(i), ny(i), ax(k), bx(k), ay(k), by(k), m, 8);
          end
        end
      %%end
    end 
  end

end 



function testepsilon2D(sten, eps0, eps1, d, nx, ny, ax, bx, ay, by, m)
%
% testepsilon2D aprroximates the modified Runge, smoothed Heaviside, and
% Gelb and Tanner functions with different values of eps0 that
% are used to bound the interpolant in the case of the PPI method. 
% This function produces the results used to build the 1D figures 
% In the manuscript.
%
% INPUT
% sten: stencil selction procedure (sten=1, sten=2, sten=3) 
% eps0: array of values of eps0 
% d:  traget polynomial degree for each interpolant
% n: number of points
% ax: left boundaries for x 
% ay: left boundaries for y 
% bx: right boundaries for x
% bx: right boundaries for y
% m: number of output points 
%
  
  v2D = zeros(nx,ny);  % input data values
  v2Dout_true = zeros(m,m); % to hold trus solution
  v2D_tmp = zeros(m,ny); % tmp variable only needed for PCHIP
  v2D_s = zeros(m*m, 10);% output to be written to file
  x = zeros(nx, 1);  % uniformly spaced points
  y = zeros(ny, 1);  % unformly spaced points
  x_lgl = zeros(nx, 1); % for lgl mesh points
  y_lgl = zeros(ny, 1); % for lgl mesh points
  xout = zeros(m, 1);  % output points 
  yout = zeros(m, 1);  % output points

  for k=1:4
    %%** calculates intreval sizes **%%
    dxn = (bx(k)-ax(k)) /double(nx-1);
    dxm = (bx(k)-ax(k)) /double(m-1);
    dyn = (by(k)-ay(k)) /double(ny-1);
    dym = (by(k)-ay(k)) /double(m-1);
  

     %%** uniform mesh **%%
     for i=1:nx
       x(i) = ax(k) + double(i-1)*dxn;
     end      
     x(nx) = bx(k);
     for i=1:ny  
       y(i) = ay(k) + double(i-1)*dyn;
     end
     y(ny) = by(k);

    %%** output mesh points **!
    for i=1:m
      xout(i) = ax(k) + double(i-1)*dxm;
      yout(i) = ay(k) + double(i-1)*dym;
    end
    xout(m) = bx(k);
    yout(m) = by(k);


    %%** only used in calculation inside of evalFun2D for fun == 4
    h = dxn;

    %%** Data values associated to input meshes **%%
    for j=1:ny
      for i=1:nx
        v2D(i,j) = evalFun2D(k, x(i), y(j));
      end
    end

    %%** True solution **%%
    for j=1:m
      for i=1:m
        v2Dout_true(i,j) = evalFun2D(k, xout(i), yout(j));
      end
    end

    ii = 1;
    for j=1:m
      for i=1:m
        v2D_s(ii, 1) = xout(i);
        v2D_s(ii, 2) = yout(j);
        v2D_s(ii, 3) = v2Dout_true(i, j);
        ii = ii+1;
      end
    end
 
    for kk=1:7
      %%**  Interpolation using Tensor product and PPI **%%
      if(kk == 7)
        tmp = adaptiveInterpolation2D(x, y, v2D, xout, yout, d, 1);
        v2Dout = tmp;
      else
        tmp = adaptiveInterpolation2D(x, y, v2D, xout, yout, d, 2);
        v2Dout = tmp;
      end

      ii = 1;
      for j=1:m
        for i=1: m
          v2D_s(ii, kk+3) = v2Dout(i, j);
          ii = ii+1;
        end
      end
    end

    %%** Open file **%% 
    if(k ==1 )
      fid = fopen('Runge2DEps', 'w');
    elseif(k ==2 )
      fid = fopen('Surface1Eps', 'w');
    elseif(k ==3 )
      fid = fopen('Surface2Eps', 'w');
    elseif(k ==4 )
      fid = fopen('Heaviside2DEps', 'w');
    end
    %%** Write to open file **%%
    ii =1;
    for j=1:m
      for i=1:m
        for kk=1:10
          fprintf(fid,'%.8E \t ', v2D_s(ii, kk) );
        end
        fprintf(fid, '\n');
        ii = ii+1;
      end
    end
    %%** close file **%%
    fclose(fid);
  end
end 

function test002(d, eps0, eps1, sten, fun, nx, ny, ax, bx, ay, by, m, d_el)
%
% test002 is used to approximate the Runge, smoothed Heaviside
% and 2D Terrain functions using different interpolation 
% methods. This function is used toproduce the 2D results presented
% in the manuscript.
%
% INPUT
% d: maximum polynomial degree for each interval
% eps0: positive user-supplied value used to bound interpolant for 
%       intervalswith no extrema.
% eps1: positive user-supplied value used to bound interpolant for 
%       intervals with extrema.
% sten: user-supplied value used to indicate stencil selection process
%       possible choices are sten=1, sten=2, sten=3.
% fun: used to indicate function used
% nx: number of point in x direction
% ny: number of points in y direction
% ax: global interval left boundary in the x direction
% bx: global right interval boundary in the x direction
% ay: global interval left boundary in the y direction
% by: global right interval boundary in the y direction
% d_el: number of lgl points in each element
%

  %% Initialize %%
  x = zeros(nx, 1);
  x_lgl = zeros(nx, 1);
  y = zeros(ny, 1);
  y_lgl = zeros(ny, 1);
  xout = zeros(m,1);
  yout = zeros(m,1);
  v2Dout = zeros(m,m);
  v2Dout_lgl = zeros(m,m);
  v2Dout_true = zeros(m,m);
  v2D_tmp = zeros(m,ny);

  %%** get function  name **%%
  if(fun ==1)
    fun_name = "Runge2D";
  elseif(fun ==2)
    fun_name = "T1";
  elseif(fun ==3)
    fun_name = "T2";
  elseif(fun ==4)
    fun_name = "Heaviside2D";
  else
    fprintf('ERROR: Invalid fun = %d \n', fun);
    fprintf('Invalid function value the possible values are fun =1, 2, 3, or 4 \n');
    stop
  end

  %%** get stencil selection procedure **%%
  if(sten ==1) 
    sst = "st=1";
  elseif(sten ==2) 
    sst = "st=2";
  elseif(sten ==3) 
    sst = "st=3";
  else
    fprintf('ERROR: Invalid paparamter sten = %d \n', fun);
    fprintf('ERROR: Invalid paparamter st. The possible options are st=1, 2, or 3 \n');
    stop;
  end
 
  if(d < 10)
    fnumber =strcat("0", string(d), string(nx), string(ny));
  else
    fnumber =strcat(string(d), string(nx), string(ny));
  end

  %%** calculates intreval sizes **%%
  dxn = (bx-ax) /double(nx-1);
  dxm = (bx-ax) /double(m-1);
  dyn = (by-ay) /double(ny-1);
  dym = (by-ay) /double(m-1);
  

  %%** unifnorm mesh **%%
  for i=1:nx
    x(i) = ax + double(i-1)*dxn;
  end
  x(nx) = bx;
  for i=1:ny
    y(i) = ay + double(i-1)*dyn;
  end
  y(ny) = by;

  %%** number of elements **%%
  dd = d_el;
  nex = (nx-1) / dd;
  ney = (ny-1) / dd;

  %%** Legendre gauss lobatto points **%%
  x_tmp = lglnodes(dd);
  x_tmp = flip(x_tmp); 
  dxn = (bx-ax) / double(nex);               %% calculates element size
  xl = ax;                                        %% initialaziation 
  xr = ax;                                      %% initialization 
  is = 1;
  ie = 1;
  for i=1:nex
    xl = xr;                                     %% left boundary of element i
    xr = xl + dxn;                                %% right boun dary of element i
    is = ie;
    ie = is + dd;
    %%** maping from [-1,1] to [xl, xr] 
    x_lgl(is:ie) = x_tmp* (xr-xl)/2.0 + (xr+xl)/2.0;
  end
  x_lgl(nx) = bx;
  dyn = (by-ay) / double(ney);               %% calculates element size
  yl = ay;
  yr = ay;
  is = 1;
  ie = 1;
  for i=1:ney
    yl = yr;                                     %% left boundary of element i
    yr = yl + dyn;                                %% right boun dary of element i
    is = ie;
    ie = is + dd;
    %%** maping from [-1,1] to [yl, yr] 
    y_lgl(is:ie) = x_tmp* (yr-yl)/2.0 + (yr+yl)/2.0;
  end
  y_lgl(ny) = by;

  %%** output mesh points **!
  for i=1:m
    xout(i) = ax + double(i-1)*dxm;
    yout(i) = ay + double(i-1)*dym;
  end
  xout(m) = bx;
  yout(m) = by;

  %%** only used in calculation inside of evalFun2D for fun == 4
  h = (bx-ax)/((nx-1)/d);

  %%** Data values associated to input meshes **%%
  for j=1:ny
    for i=1:nx
      v2D(i,j) = evalFun2D(fun, x(i), y(j));
      v2D_lgl(i,j) = evalFun2D(fun, x_lgl(i), y_lgl(j));
    end
  end

  %%** True solution **%%
  for j=1:m
    for i=1:m
      v2Dout_true(i,j) = evalFun2D(fun, xout(i), yout(j));
    end
  end
 
  %%**  Interpolation using Tensor product and PCHIP **%%
  if(d == 3)  
    for j=1:ny
      tmp =pchip(x, v2D(:,j), xout);
      v2D_tmp(:,j) = tmp;
    end
    for i=1:m
      tmp = pchip(y, v2D_tmp(i,:), yout);
      v2Dout(i,:) = tmp;
    end
    %%** interpolation using lgl **%%
    for j=1:ny
      tmp =pchip(x_lgl, v2D_lgl(:,j), xout);
      v2D_tmp(:,j) = tmp;
    end
    for i=1:m
      tmp = pchip(y_lgl, v2D_tmp(i,:), yout);
      v2Dout_lgl(i,:) = tmp;
    end

    %%** Open file **%% 
    fname = strcat(fun_name, "PCHIP", fnumber);
    fid = fopen(fname, 'w');
    %%** Write to open file **%%
    for j=1:m
      for i=1:m
        fprintf(fid, '%.8E \t %.8E \t %.8E \t %.8E \t %.8E \n', ...
                xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j), v2Dout_lgl(i, j) );
      end
    end
    %%** close file **%%
    fclose(fid);
  end
  %%**  Interpolation using Tensor product and DBI **%%
  v2Dout = adaptiveInterpolation2D(x, y, v2D, xout, yout, d, 1);
  v2Dout_lgl = adaptiveInterpolation2D(x_lgl, y_lgl, v2D_lgl, xout, yout, d, 1) ;

  %%** Open file **%% 
  fname = strcat(fun_name, "DBI", fnumber, sst);
  fid =fopen(fname, 'w');
  %%** Write to open file **%%
  for j=1: m
    for i=1:m
      fprintf(fid, '%.8E \t %.8E \t %.8E \t %.8E \t %.8E \n', ...
              xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j), v2Dout_lgl(i, j) );
    end
  end
  %%** close file **%%
  fclose(fid)

  %%**  Interpolation using Tensor product and PPI **%%
  v2Dout = adaptiveInterpolation2D(x, y, v2D,  xout, yout, d, 2);
  v2Dout_lgl = adaptiveInterpolation2D(x_lgl, y_lgl, v2D_lgl,  xout, yout, d, 2);

  %%** Open file **%% 
  fname = strcat(fun_name, "PPI", fnumber, sst);
  fid=fopen(fname, 'w');
  %%** Write to open file **%%
  for j=1: m
    for i=1:m
      fprintf(fid, '%.8E \t %.8E \t %.8E \t %.8E \t %.8E \n', ...
              xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j), v2Dout_lgl(i, j) );
    end
  end
  %%** close file **%%
  fclose(fid);

end 


%---------------------------------------------------------------------------------------------%
% Functions used for the mapping examples in the manuscript
%---------------------------------------------------------------------------------------------%
function  mapping(nz)
%
% This subrtouine is used to set up the mapping for the Runge and TWP-ICE examples.
% The following files below are required for the experiment.
% 'zd_qc_qv_pres_u_v_T_zp_127' and 'zd_qc_qv_pres_u_v_T_zp_253' are obtained by fitting
% 'zd_qc_qv_pres_u_v_T_zp_64' using a radial basis function interpolation and the evaluating
% the fitted function at the desired points.
%
% FILES
% 'mapping_data/zd_qc_qv_pres_u_v_T_zp_64': obtained directly from TWP-ICE simulation at
%   at t = XX s.
% 'mapping_data/zd_qc_qv_pres_u_v_T_zp_127': obtained by adding at point at the center of each interval 
% 'mapping_data/zd_qc_qv_pres_u_v_T_zp_253': obtained by adding 3 uniformly spaced points inside each 
%   interval.
%  
% INPUT
% nz: number of point to be used for the Runge and TWP-ICE examples 
%


  fprintf('Running mapping examples ...  ... \n')
  %%** Read input data from file  **%%
  if(nz == 64) 
    dd = load('mapping_data/zd_qc_zp_64');
  elseif(nz == 127) 
    dd = load('mapping_data/zd_qc_zp_127');
  elseif(nz == 253) 
    dd = load('mapping_data/zd_qc_zp_253');
  end
  zd = dd(:,1);
  qc2 = dd(:,2);
  zp = dd(:,3);
  qcp2 = dd(:,4);

  %%** Initialize polynomial degree to be used **%%
  d = [3, 5, 7];

  %%** Initialize variables names **%%
  name_runge = 'runge';
  name_qc = 'qc';

  %%%% Map zd and zp  to xd_xx an xp_xx respectively 
  a_runge = -1.0;
  b_runge = 1.0;
  zd_runge = scaleab(zd, zd(1), zd(nz), a_runge, b_runge);
  zp_runge = scaleab(zp, zd(1), zd(nz), a_runge, b_runge);

  %% Evaluate coresponding values at x 
  dx = 0.01; %% dummy variables not used for function evaluation
  for i=1:nz
    runge2(i) = evalFun1D(1, zd_runge(i), dx);
    rungep2(i) = evalFun1D(1, zp_runge(i), dx);
  end

  for j=1:3
    for i=1:3
      runge = runge2;
      rungep = rungep2;
      qc=qc2;
      qcp=qcp2;
      fprintf('********** d= %d ******** \n', d(i) );
      mapping2(zd_runge, runge, zp_runge, rungep, d(i), j, name_runge);
      mapping2(zd, qc, zp, qcp, d(i), j, name_qc);
    end
  end

end 

function  mapping2(zd, u, zp, u2, dd, st, profile_name)
%
% Subroutine for mapping data form mesh points zd to zp and back to zp
%
% INPUT
% nz: number of points
% zd: first mesh points (dynamics mesh points)
% u: data values associated with the first mesh
% zp: fecond mesh points (physics mesh points)
% u2: data values associated with the second mesh
% dd: maximum degree used fr each interpolant
% profile_name: profile name to be used to save results
%
%

  nz = length(zd);
  ud = zeros(nz,1);
  ud_pchip = zeros(nz,1);
  ud_dbi = zeros(nz,1);
  up = zeros(nz,1);
  up_pchip = zeros(nz,1);
  up_dbi = zeros(nz,1);

  ud_out = zeros(nz,3);
  up_out = zeros(nz,3);
  ud_pchip_out = zeros(nz,3);
  up_pchip_out = zeros(nz,3);
  ud_dbi_out = zeros(nz,3);
  up_dbi_out = zeros(nz,3);
  deg_ud_dbi_out = zeros(nz-1, 3);
  deg_up_dbi_out = zeros(nz-1, 3);
  deg_ud_out = zeros(nz-1, 3);
  deg_up_out = zeros(nz-1, 3);
  deg = zeros(nz-1,1);
  deg_dbi = zeros(nz-1,1);
  
  eps0 = 0.01;
  eps1 = 1.00;

  %%** set stencil type for file name **%%
  if(st==1) 
   sst = "st=1";
  elseif(st==2) 
   sst = "st=2";
  elseif(st==3) 
   sst = "st=3";
  end
 
  %%** Initialize file number **%%
  if (dd < 10) 
    fnumber = strcat( "0", string(dd*1000 + nz) );
  else
    fnumber = strcat( string(dd*1000 + nz) );
  end
 
  ud = u;
  ud_pchip = u;
  ud_dbi = u;

  %%** save Initial data **%%
  for i=1:nz
    ud_out(i,1) = zd(i);
    up_out(i,1) = zp(i);
    ud_pchip_out(i,1) = zd(i);
    up_pchip_out(i,1) = zp(i);
    ud_dbi_out(i,1) = zd(i);
    up_dbi_out(i,1) = zp(i);
  end
  %
  for i=1:nz-1
    deg_ud_out(i,1) = 0;
    deg_up_out(i,1) = 0;
    deg_ud_dbi_out(i,1) = 0;
    deg_up_dbi_out(i,1) = 0;
  end

  iter=1;
  %%** save data on physics grid **%%
  for i=1:nz
    up_out(i,iter+1) = u2(i);
    up_pchip_out(i,iter+1) = u2(i);
    up_dbi_out(i,iter+1) = u2(i);
  end
  %%
  for i=1:nz-1;
    deg_up_out(i,iter+1) = 0;
    deg_up_dbi_out(i,iter+1) = 0;
  end
   
  %%** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using DBI  **%%
  up_dbi = adaptiveInterpolation1D(zd, ud_dbi, zp, dd, 1, st, eps0, eps1); 
  
  %%** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using PPI  **%%
  up = adaptiveInterpolation1D(zd, ud, zp, dd, 2, st, eps0, eps1);

  %%** Mapping data values from zd (dynamics mesh) to zp (physics mesh) using PCHIP  **%%
  up_pchip = pchip(zd, ud_pchip, zp);

  %%** save data that is on dynamics grid**%%
  for i=1:nz
    ud_out(i,iter+1) = ud(i); 
    ud_pchip_out(i,iter+1) = ud_pchip(i);
    ud_dbi_out(i,iter+1) = ud_dbi(i);
  end
  %% save polynomial degrees used for each interval %%
  for i=1:nz-1
    deg_ud_out(i,iter+1) = deg(i); 
    deg_ud_dbi_out(i,iter+1) = deg_dbi(i);
  end

  iter = 2;
  %%** save data on physics grid **%%
  for i=1:nz
    up_out(i,iter+1) = up(i);
    up_pchip_out(i,iter+1) = up_pchip(i);
    up_dbi_out(i,iter+1) = up_dbi(i);
  end

  %%** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using DBI  **%%
  ud_dbi(2:nz-1) = adaptiveInterpolation1D(zp, up_dbi, zd(2:nz-1), dd, 1, st, eps0, eps1); 
  
  %%** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using PPI  **%%
  ud(2:nz-1) = adaptiveInterpolation1D(zp, up, zd(2:nz-1), dd, 2, st, eps0, eps1);

  %%** Mapping data values from zp (physics mesh) to zd (dynamics mesh)  using PCHIP  **%%
  ud_pchip(2:nz-1) = pchip( zp, up_pchip, zd(2:nz-1) );

  iter = 2;
  %%** save data on physics grid **%%
  for i=1:nz-1
    deg_up_out(i,iter+1) = deg(i);
    deg_up_dbi_out(i,iter+1) = deg_dbi(i);
  end
  
  %%** save data that is on dynamics grid **%%
  for i=1:nz
   ud_out(i,iter+1) = ud(i); 
   ud_pchip_out(i,iter+1) = ud_pchip(i);
   ud_dbi_out(i,iter+1) = ud_dbi(i);
  end
  %%
  for i=1:nz-1
   deg_ud_out(i,iter+1) = 0; 
   deg_ud_dbi_out(i,iter+1) = 0;
  end


  %% Save PPI results to file %%
  fname = strcat(profile_name, "dPPI", fnumber, sst); % file name 
  fid = fopen(char(fname), 'w'); % open file 
  for i=1:nz
    fprintf(fid, '%.8E \t %.8E \t %.8E \n', ud_out(i,1), ud_out(i,2), ud_out(i,3) );
  end
  fclose(fid);
  %
  fname = strcat(profile_name, "dDEGPPI", fnumber, sst);  % file name
  fid = fopen(char(fname), 'w');  % open file 
  for i=1:nz-1
    fprintf(fid, '%d \t %d \t %d \n', deg_ud_out(i,1), deg_ud_out(i,2), deg_ud_out(i,3) );
  end
  fclose(fid);  % close file 
  %
  fname = strcat(profile_name, "pPPI", fnumber, sst);
  fid = fopen(char(fname), 'w');
  for i=1:nz
    fprintf(fid, '%.8E \t %.8E \t %.8E \n', up_out(i,1), up_out(i,2), up_out(i,3) );
  end
  fclose(fid);
  %
  fname = strcat(profile_name, "pDEGPPI", fnumber, sst);
  fid = fopen(char(fname), 'w');
  for i=1:nz-1
    fprintf(fid, '%d \t %d \t %d \n', deg_up_out(i,1), deg_up_out(i,2), deg_up_out(i,3) );
  end
  fclose(fid);
  %
  %% Save DBI results to file %%
  fname = strcat(profile_name, "dDBI", fnumber, sst);
  fid = fopen(char(fname), 'w');
  for i=1:nz
    fprintf(fid, '%.8E \t %.8E \t %.8E \n', ud_dbi_out(i,1), ud_dbi_out(i,2), ud_dbi_out(i,3) );
  end
  fclose(fid);
  %
  fname = strcat(profile_name, "dDEGDBI", fnumber, sst);
  fid = fopen(char(fname), 'w');
  for i=1:nz-1
    fprintf(fid, '%d \t %d \t %d \n', deg_ud_dbi_out(i,1), deg_ud_dbi_out(i,2), deg_ud_dbi_out(i,3) );
  end
  fclose(fid);
  %
  fname = strcat(profile_name, "pDBI", fnumber, sst);
  fid = fopen(char(fname), 'w');
  for i=1:nz
    fprintf(fid, '%.8E \t %.8E \t %.8E \n', up_dbi_out(i,1), up_dbi_out(i,2), up_dbi_out(i,3) );
  end
  fclose(fid);
  %
  fname = strcat(profile_name, "pDEGDBI", fnumber, sst);
  fid = fopen(char(fname), 'w');
  for i=1:nz-1
    fprintf(fid, '%d \t %d \t %d \n', deg_up_dbi_out(i,1), deg_up_dbi_out(i,2), deg_up_dbi_out(i,3) );
  end
  fclose(fid);
  %
  %% Save PCHIP results to file %%
  fname = strcat(profile_name, "dPCHIP", fnumber);
  fid = fopen(char(fname), 'w');
  for i=1:nz
    fprintf(fid, '%.8E \t %.8E \t %.8E \n', ud_pchip_out(i,1), ud_pchip_out(i,2), ud_pchip_out(i,3) );
  end
  fclose(fid);
  %
  fname = strcat(profile_name, "pPCHIP", fnumber);
  fid = fopen(char(fname), 'w');
  for i=1:nz
    fprintf(fid, '%.8E \t %.8E \t %.8E \n', up_pchip_out(i,1), up_pchip_out(i,2), up_pchip_out(i,3) );
  end
  fclose(fid);
  %

end 


function v = evalFun1D(fun, x, h)
%
% To evaluate different function. fun determine
% the function that will be evaluated at the points x
%
% INPUT
% fun: function type
% x: function input values
% h: function parameter 
%
  
  %%** intialize variables **%%
  k = 100;
  pi = 4.0*atan(1.0); 
  t = x/h;
  delta = 0.01;

  %%** a and b are only used for fun =5. When using fun=4 the 
  %%   a and b must be the same as the ones Paperexample1D **%%
  a = -2.0;
  b = 0.0;
  ne = int16((b-a)/h);
  

  %%!** 1D runge function **%%
  if(fun == 1)
    v = 0.1 / (0.1 + 25.0 * x * x);
    %%v = 1.0 / (1.0 + 25.0 * x * x);
  %%** heaviside function **%%
  elseif(fun == 2)
    v = 1.0/(1.0 + exp(-2*k*x));
  %%** Gelb and Tanner function **%%
  elseif(fun == 3)
    if(x < -0.5) 
      v = 1.0 + (2.0* exp(2.0 * pi *(x+1.0)) - 1.0 -exp(pi)) /(exp(pi)-1.0);
    else
      v = 1.0 - sin(2*pi*x / 3.0 + pi/3.0);
    end
  %%%%** modified square function **%%
  elseif(fun == 4)
    v = 1.0 - abs( 2.0/ pi * atan( sin(pi*t) /delta));
  %%** Modified tanh function **%%
  elseif(fun == 5)
    k = 10.0;
    if(a <= x && x <= a+h)
      v = tanh(x*k);
    elseif(a+h<= x && x <= a+2*h)
      v = 2*tanh(x*k)         - tanh((a+h)*k);
    end
    for i=3:ne 
      %%** 
      if(a +(i-1)*h <= x && x <= a+i*h)
        v = i*tanh(x*k) - tanh((a+h)*k);
        for j=2:i-1 
         v = v - tanh( (a+(j*h))*k );
        end 
      end
    end
    v = 1.0 + v;
  %%%%* sin function **%%
  elseif(fun == 6)
    v = 1.0 + sin(x*pi);
  end

end 



function v =  evalFun2D(fun, x, y)
%
% To evaluate different fuunction. fun determine
% the function that will be evaluated at the points x
% fun: function type
% x: function input values
% y: function input values
%
 
  %%** intialize variables **%%
  k = 100;
  pi = 4.0*atan(1.0) ;
  delta = 0.01;

  %%** 1D runge function **%%
  if(fun == 1) 
    v = 1.0 / ( 1.0 + 25.0 * ( x*x + y*y) );
  %%** **%%
  elseif(fun == 2)
    if( (x-1.5)*(x-1.5) + (y-0.5)*(y-0.5) <= 1.0/16.0 )
      v = 2.0*0.5*( cos(8.0*atan(1.0)*sqrt((x-1.5)*(x-1.5) ...
          + (y-0.5)*(y-0.5))));%%+1)
    elseif(y-x >= 0.5)
      v = 1.0;
    elseif(0.0 <= y-x && y-x <= 0.5)
      v = 2.0*(y-x);
    else
      v = 0.0;
    end
  %%** **%%
  elseif(fun == 3)
    v = max(0.0, sin(4.0*atan(1.0)*x)*sin(4.0*atan(1.0)*y) ) ;
  %%** Smoothed Heaviside function
  elseif(fun == 4)
    v = 1.0 / ( 1.0 + exp(-2*k*(x+y)*sqrt(2.0)/2.0) );
  end

end 


function vout = scaleab(vin, v_min, v_max, a, b)
% Scale input data from [v_min v_max] to [a, b]
%
% INPUT:
% n:      number of input elements 
% vin(n): Input data of size n
% v_min:  left boundary of input interval
% v_max:  right boundary of input interval
% a:      left boundary of output interval
% b:      right boundary of output interval
%
% OUTPUT:
% vout(n): output data of size n
%

  n = length(vin);
  vout = zeros(n,1);
  for i=1:n
     %% map from [v_min, v_max] to [a, b]
     %% \forall x \in [v_min, v_max], 
     %% map(x) = a + (b-a)/(v_max -v_min)*(x-v_min) 
     vout(i) = a + (b-a)/(v_max-v_min)*(vin(i)-v_min);
  end
end 


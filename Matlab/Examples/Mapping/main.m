%
%
%
clear;
close all;
clc;

%
% Driver to produce the mapping result based on the
% Runge and the TWP-ICE.
%

approximations1D();
movefile Runge* mapping_data/data
movefile Heavi* mapping_data/data
movefile GelbT* mapping_data/data

%approximations2D()

  %%%nz = (/64, 127, 253/)
  %%%do k=1,3
  %%%  call mapping(nz(k))
  %%%end
 


function approximations1D()
%
% approximation1D is used to set up the 
% diffferent configuration used to produces
% the approximation results for the 1D functions
% presented in the manuscript.
%


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
    fprintf('st= %d \n', sten(ii) );
    for k=1:3

      %%** Third order resulst using DBI, PPI, and PCHIP **%%
      fprintf('*****  fun= %d ******** \n', fun(k) ) ;
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
%%
%% testepsilon1D aprroximates the Runge, smoothed Heaviside, and
%% Gelb and Tanner functions with different values of eps0 that
%% are used to bound the interpolant in the case of the PPI method. 
%%
%% INPUT
%% sten: stencil selction procedure (sten=1, sten=2, sten=3) 
%% eps0(6): array of values of eps0 
%% d:  traget polynomial degree for each interpolant
%% n: number of points
%% a(3): left boundaries
%% b(3): right boundaries
%% m: number of output points 
%%

%  use mod_adaptiveInterpolation
%
%  implicit none
%
%  integer, intent(in)           :: n                    %% number of input points
%  integer, intent(in)           :: m                    %% number of output points
%  integer, intent(in)           :: d                    %% target interpolant degree
%  integer, intent(in)           :: sten                 %% stencil selection procedure
%  double(kind=8), intent(in)      :: a(3), b(3)           %% interval [a, b]
%  double(kind=8), intent(in)      :: eps0(6), eps1        %% test values used for eps0
%
%  integer                       :: degOut(n-1, 7)       %% degree used for each interval
%  integer 		        :: i, j, k, fid
%  double(kind=8)			:: x(n)                 %% uniform input mesh points  
%  double(kind=8)			:: v1D(n)               %% input data values
%  double(kind=8)			:: v1Dout(m, 9)         %% output values
%  double(kind=8)			:: dxn, dxm             %% interval sizes 



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

    %%** output mesh points **!
    for i=1:m
      v1Dout(i, 1) = a(k) + double(i-1)*dxm;
    end

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
        tmp = adaptiveInterpolation1D(x, v1D, v1Dout(:,1), d, 1 );
        v1Dout(:,2+i) = tmp; 
      else
        tmp = adaptiveInterpolation1D(x, v1D, v1Dout(:,1), d, 2 ); 
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
%
%    %%** open file **%%
%    fid = 10
%    if( k == 1)
%      open(unit=fid, file='RungeEpsDeg', status='unknown')
%    elseif( k == 2)
%      open(unit=fid, file='HeavisideEpsDeg', status='unknown')
%    elseif( k == 3)
%      open(unit=fid, file='GelbTEpsDeg', status='unknown')
%    end
%
%    %%** write to file **%%
%    do i=1, n-1
%      write(fid,'(7(1x, I2))') ( degOut(i, j), j=1, 7 )
%    end
%    %%** close file **%%
%    close(fid)
  end
end 

function test001(d, eps0, eps1, sten, fun, n, a, b, m, d_el)
%% 
%% test001 is used to approximate the Runge, smoothed Heaviside
%% and Gelbd and Tanner function using different interpolation 
%% methods.
%%
%% INPUT
%% d: maximum polynomial degree for each interval
%% eps:
%% 

%  use mod_legendre
%  use mod_adaptiveInterpolation
%
%  implicit none
%
%  integer, intent(in)           :: fun                  %% function type 
%  integer, intent(in)           :: n                    %% number of input points
%  integer, intent(in)           :: m                    %% number of output points
%  integer, intent(in)           :: d                    %% target interpolant degree
%  integer, intent(in)           :: sten
%  double(kind=8), intent(in)      :: a                    %% left bounary
%  double(kind=8), intent(in)      :: b                    %% right boundary
%  double(kind=8), intent(in)      :: eps0			%% parameters used to bound interpolant in intervals with no hidden extrema  
%  double(kind=8), intent(in)      :: eps1			%% parameters used to bound interpolant in intervals with hidden extrema  
%  integer, intent(in)           :: d_el
%
%  integer                       :: degOut(n-1, 2)      %% degree used for each interval
%  integer                       :: degOut_lgl(n-1, 2)      %% degree used for each interval
%  integer                       :: limiter             %%
%  integer                       :: ne                   %% number of elments
%  integer 		        :: i, j, k, fid, ierr, tmp_idx
%  integer 		        :: is, ie, dd
%  double(kind=8)			:: x(n), x_lgl(n)                     %% uniform and  LGL input mesh points  
%  double(kind=8)			:: v1D(n), v1D_lgl(n)              %% input data values
%  double(kind=8)			:: xout(m)                                      %% output points to be approximated 
%  double(kind=8)			:: v1Dout(m), v1Dout_lgl(m)     %% approximated output values
%  double(kind=8)			:: v1Dout_true(m)                               %% True values at output points
%  double(kind=8)			:: dxn, dxm
%  double(kind=8)			:: x_tmp(d_el+1), w_tmp(d_el+1), xl, xr
%  character*16                  :: fun_name
%  character*16                  :: fnumber
%  character*16                  :: sst
%  character*64                  :: fname
%
%  %%** Local variables need for PCHIP **%%
%  integer 		        :: nwk
%  double(kind=8)			:: wk((n+1)*2), d_tmp(n+1)
%  double(kind=8)			:: fdl(m)
%  logical                       :: spline


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
    %quit(0);
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
    %quit(0);
  end
  
  %%** uniform mesh **%%
  dxn = (b-a) /double(n-1);
  for i=1:n
    x(i) = a + double(i-1)*dxn;
  end

  dd = d_el;
  ne = (n-1) / dd;                        %% calculates the number of elements
  
  %%** LGL mesh **%%
  x_tmp = lglnodes(dd);
  x_tmp = flip(x_tmp); 
  dxn = (b-a) / double(ne);                          %% calculates element size
  xl = a;                                                    %% initialaze element left boundary 
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

  %%** output mesh points **!
  dxm = (b-a) /double(m-1);
  for i=1:m
    xout(i) = a + double(i-1)*dxm;
  end

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
  tmp = adaptiveInterpolation1D(x, v1D, xout, d, 1); 
  v1Dout = tmp;
  tmp = adaptiveInterpolation1D(x_lgl, v1D_lgl, xout, d, 1); 
  v1Dout_lgl = tmp;

  %%** open file and write to file **%%
  fname = strcat( fun_name, "DBI", fnumber, sst);
  fid = fopen(char(fname), 'w');
  for i=1:m
    fprintf(fid,' %.8E \t %.8E \t %.8E \t %.8E \n ', xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i));
  end
  fclose(fid);

  %%** Interpolation using PPI **%%
  tmp = adaptiveInterpolation1D(x, v1D, xout, d, 2); 
  v1Dout = tmp;
  tmp = adaptiveInterpolation1D(x_lgl, v1D_lgl, xout, d, 2 ); 
  v1Dout_lgl = tmp;

  %%** open file and write to file **%%
  fname = strcat( fun_name, "PPI", fnumber, sst);
  fid = fopen(char(fname), 'w');
  for i=1:m
    fprintf(fid,' %.8E \t %.8E \t %.8E \t %.8E \n ', xout(i), v1Dout_true(i), v1Dout(i), v1Dout_lgl(i));
  end
  fclose(fid);

  %%%** open file write degree used for each interpolant to file**%%
  %fname = strcat( fun_name, "PPI", fnumber, sst);
  %fid = fopen(char(fname), 'w');
  %for i=1:m
  %  fprintf(fid,' %d \t %d \t %d \t %d \n ', degOut(i, 1), degOut(i, 2), degOut_lgl(i, 1), degOut_lgl(i, 2) );
  %end
  %fclose(fid);

end 


%% 2D Examples
function approximations2D()
%%
%%
%%
%  implicit none
%
%  integer                       :: nx(5), nx2(5)
%  integer                       :: ny(5), ny2(5)
%  integer                       :: d(4)
%  integer                       :: fun(4)
%  integer                       :: i, ii, j, k, kk
%  integer                       :: sten(3)
%  integer, parameter            :: m = 100%%1000
%  double(kind=8)                  :: ax(4), bx(4)
%  double(kind=8)                  :: ay(4), by(4)
%  double(kind=8)                  :: eps0, eps1, eps_test(6)


  d = [1, 4, 8, 16];                                                       %% array with interpolants degrees
  nx = [17, 33, 65, 129, 257];                                                %% array with number of inputpoints
  ny = [17, 33, 65, 129, 257];                                                %% array with number of inputpoints
  nx2 =[13, 25, 49, 97, 193];                                                %% array with number of inputpoints
  ny2 =[13, 25, 49, 97, 193];                                                %% array with number of inputpoints

  eps_test = [ 1.0,  0.1,  0.01, 0.001, 0.0001, 0.00 ];

  %%** set up interval x \in [ax(i), bx(i)] and y \in [ay(i), by(i)]**%% 
  ax = [-1.0, -1.0, 0.0, -0.2 ];
  bx = [ 1.0,  1.0, 2.0, 0.2  ];
  ay = [-1.0, -1.0, 0.0, -0.2 ];
  by = [ 1.0,  1.0, 1.0, 0.2  ];
  
  m =  200;
  %%** function type 1=runge funtion , 2= heaviside, 3=Gelb Tanner **%% 
  fun = [1, 2, 3, 4];                                                     %% function type

  %%**
  sten = [1, 2, 3];
  eps0 = 0.01;
  eps1 = 1.0;
  testepsilon2D(sten(2),eps0, eps1, d(3), nx(1), ny(1), ax, bx, ay, by, m);
  %for ii=1:3
  %  %%** comparing against PCHIP **%%
  %  for k=1:4 
  %    %%if(k == 1 .or. k == 4) 
  %      for i=1:5
  %         %%** Perform interpolation and calculate different L2-error
  %         %!    norms different methods **%%
  %         fprintf('nx=%d \t ny=%d \n', nx(i), ny(i) );
  %         fprintf('ax=%d \t bx=%d \n', ax(k), bx(k) );
  %         fprintf('ay=%d \t by=%d \n', ay(k), by(k) );
  %         test002(3, eps0, eps1, sten(ii), fun(k), nx(i), ny(i), ax(k), bx(k), ay(k), by(k), m, 8)
  %      end


  %      fprintf('*****  fun=%d *****', fun(k) );
  %      for j=1:4  
  %        fprintf('*****  d= %d ***** \n', d(j) );
  %        for i=1:5
  %           %%** Perform interpolation and calculate different L2-error
  %           %%    norms different methods **%%
  %           fprintf('nx= %d \t ny=%d \n', nx(i), ny(i) );
  %           fprintf('ax(k)= %d \t bx(k) = %d \n', ax(k), bx(k) );
  %           fprintf('ay(k)= %d \t bx(k) = %d \n', ay(k), by(k) );
  %            test002(d(j), eps0, eps1, sten(ii), fun(k), nx(i), ny(i), ax(k), bx(k), ay(k), by(k), m, 8)
  %        end
  %      end
  %    %%end
  %  end 
  %end

end 



function testepsilon2D(sten, eps0, eps1, d, nx, ny, ax, bx, ay, by, m)
%%
%%
%%

%  use mod_legendre
%  use mod_adaptiveInterpolation
%
%
%  implicit none
%
%
%  integer, intent(in)           :: nx                    %% number of input points
%  integer, intent(in)           :: ny                    %% number of input points
%  integer, intent(in)           :: m                    %% number of output points
%  integer, intent(in)           :: d                    %% target interpolant degree
%  integer, intent(in)           :: sten
%  double(kind=8), intent(in)      :: ax(4), bx(4)                 %% interval [a, b]
%  double(kind=8), intent(in)      :: ay(4), by(4)                 %% interval [a, b]
%  double(kind=8), intent(in)      :: eps0(6), eps1
%
%  integer 		        :: i, ii, kk, j, k, fid
%  double(kind=8)			:: x(nx)                 %% input mesh points  
%  double(kind=8)			:: y(ny)                 %% input mesh points  
%  integer 			:: degx2(nx-1, ny)       %% input mesh points  
%  integer 			:: degy2(ny-1, m)        %% input mesh points  
%  double(kind=8)			:: v2D(nx, ny)           %% input data values
%  double(kind=8)			:: xout(m)               %% output points to be approximated 
%  double(kind=8)			:: yout(m)               %% output points to be approximated 
%  double(kind=8)			:: v2Dout(m, m)          %% approximated output values
%  double(kind=8)			:: v2Dout_true(m, m)       %% True values at output points
%  double(kind=8)			:: v2D_tmp(m, ny)        %% True values at output points
%
%  double(kind=8)			:: v2D_s(m*m, 10)        %% True values at output points
%  double(kind=8)			:: dxn, dxm, dyn, dym
%  double(kind=8)			:: h                    %% element spacing
%
%  
%  character*36                  :: filename

  v2D = zeros(nx,nx);
  v2Dout = zeros(m,m);
  v2Dout_true = zeros(m,m);
  v2D_tmp = zeros(m,ny);
  v2D_s = zeros(m*m, 10);
  x = zeros(nx, 1);
  y = zeros(ny, 1);
  xout = zeros(m, 1);
  yout = zeros(m, 1);

  for k=1:4
    %%** calculates intreval sizes **%%
    dxn = (bx(k)-ax(k)) /double(nx-1);
    dxm = (bx(k)-ax(k)) /double(m-1);
    dyn = (by(k)-ay(k)) /double(ny-1);
    dym = (by(k)-ay(k)) /double(m-1);
  

     %%** uniform mesh **%%
     for i=1,nx
       x(i) = ax(k) + double(i-1)*dxn;
     end      
     for i=1:ny  
       y(i) = ay(k) + double(i-1)*dyn;
     end

    %%** output mesh points **!
    for i=1:m
      xout(i) = ax(k) + double(i-1)*dxm;
      yout(i) = ay(k) + double(i-1)*dym;
    end


    %%** only used in calculation inside of evalFun2D for fun == 4
    h = dxn;

    %%** Data values associated to input meshes **%%
    for j=1:ny
      for i=1:nx
        v2D(i,j) = evalFun2D(k, x(i), y(j), h);
      end
    end

    %%** True solution **%%
    for j=1:m
      for i=1:m
        v2Dout_true(i,j) = evalFun2D(k, xout(i), yout(j), h);
      end
    end

    ii = 1
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
    fclose(fid)
  end


end 

function test002(d, eps0, eps1, sten, fun, nx, ny, ax, bx, ay, by, m, d_el)
%%
%%
%%
  %%** get function  name **%%
  if(fun ==1)
    fun_name = "Runge2D"
  elseif(fun ==2)
    fun_name = "T1"
  elseif(fun ==3)
    fun_name = "T2"
  elseif(fun ==4)
    fun_name = "Heaviside2D"
  else
    fprintf('ERROR: Invalid fun = %d \n', fun);
    fprintf('Invalid function value the possible values are fun =1, 2, 3, or 4 \n');
    %quit(0)
  end

  %%** get stencil selection procedure **%%
  if(sten ==1) 
    sst = "st=1"
  elseif(sten ==2) 
    sst = "st=2"
  elseif(sten ==3) 
    sst = "st=3"
  else
    fprintf('ERROR: Invalid paparamter sten = %d \n', fun);
    fprintf('ERROR: Invalid paparamter st. The possible options are st=1, 2, or 3 \n');
    %quit(0)
  end
 

  %%** calculates intreval sizes **%%
  dxn = (bx-ax) /double(nx-1)
  dxm = (bx-ax) /double(m-1)
  dyn = (by-ay) /double(ny-1)
  dym = (by-ay) /double(m-1)
  

  %%** unifnorm mesh **%%
  for i=1,nx
    x(i) = ax + double(i-1)*dxn
  end
  for i=1:ny
    y(i) = ay + double(i-1)*dyn
  end

  %%** number of elements **%%
  dd = d_el
  nex = (nx-1) / dd
  ney = (ny-1) / dd

  %%** Legendre gauss lobatto points **%%
  x_tmp = lglnodes(dd);
  x_tmp = flip(x_tmp); 
  dxn = (bx-ax) / double(nex)               %% calculates element size
  xl = ax                                        %% initialaziation 
  xr = ax                                       %% initialization 
  is = 1
  ie = 1
  for i=1:nex
    xl = xr                                     %% left boundary of element i
    xr = xl + dxn                                %% right boun dary of element i
    is = ie
    ie = is + dd
    %%** maping from [-1,1] to [xl, xr] 
    x_lgl(is:ie) = x_tmp* (xr-xl)/2.0 + (xr+xl)/2.0
  end
  dyn = (by-ay) / double(ney)               %% calculates element size
  yl = ay
  yr = ay
  is = 1
  ie = 1
  for i=1:ney
    yl = yr                                     %% left boundary of element i
    yr = yl + dyn                                %% right boun dary of element i
    is = ie
    ie = is + dd
    %%** maping from [-1,1] to [yl, yr] 
    y_lgl(is:ie) = x_tmp* (yr-yl)/2.0 + (yr+yl)/2.0
  end

  %%** output mesh points **!
  for i=1:m
    xout(i) = ax + double(i-1)*dxm
    yout(i) = ay + double(i-1)*dym
  end


  %%** only used in calculation inside of evalFun2D for fun == 4
  h = (bx-ax)/((nx-1)/d);

  %%** Data values associated to input meshes **%%
  for j=1:ny
    for i=1:nx
      v2D(i,j) = evalFun2D(fun, x(i), y(j), h);
      v2D_lgl(i,j) = evalFun2D(fun, x_lgl(i), y_lgl(j), h);
    end
  end

  %%** True solution **%%
  for j=1:m
    for i=1:m
      v2Dout_true(i,j) = evalFun2D(fun, xout(i), h);
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
      v2Dout = tmp;
    end
    %%** interpolation using lgl **%%
    for j=1:ny
      tmp =pchip(x_lgl, v2D_lgl(:,j), xout);
      v2D_tmp(:,j) = tmp;
    end
    for i=1:m
      tmp = pchip(y_lgl, v2D_tmp(i,:), yout);
      v2Dout_lgl = tmp;
    end

    %%** Open file **%% 
    fname = strcat(fun_name, "PCHIP", fnumber);
    fopen(fname, 'w');
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
  v2Dout =0.0;
  v2Dout_lgl =0.0;
  v2Dout = adaptiveInterpolation2D(x, y, v2D, xout, yout, d, 1);
  v2Dout_tmp = adaptiveInterpolation2D(x_lgl, y_lgl, v2D_lgl, xout, yout, d, 1) ;

  %%** Open file **%% 
  fname = strcat(fun_name, "DBI", fnumber, sst);
  fopen(fname, 'w');
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
  v2Dout = 0.0
  v2Dout_lgl = 0.0
  v2Dout = adaptiveInterpolation2D(x, y, v2D,  xout, yout, d, 2);
  v2Dout = adaptiveInterpolation2D(x_lgl, y_lgl, v2D_lgl,  xout, yout, d, 2);

  %%** Open file **%% 
  fname = strcat(fun_name, "PPI", fnumber, sst);
  fopen(fname, 'w');
  %%** Write to open file **%%
  for j=1: m
    for i=1:m
      fprintf(fid, '%.8E \t %.8E \t %.8E \t %.8E \t %.8E \n', ...
              xout(i), yout(j), v2Dout_true(i, j), v2Dout(i, j), v2Dout_lgl(i, j) );
    end
  end
  %%** close file **%%
  fclose(fid)

end 



function v = evalFun1D(fun, x, h)
%%
%% To evaluate different function. fun determine
%% the function that will be evaluated at the points x
%%
  
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
    v = 1.0 + sin(x*pi) ;
  end

end 


function v =  evalFun2D(fun, x, y, h)
%%
%% To evaluate different fuunction. fun determine
%% the function that will be evaluated at the points x
%%
%  implicit none 
%
%  integer, intent(in)           :: fun                  %% function type
%  double(kind=8), intent(in)      :: x                    %% point
%  double(kind=8), intent(in)      :: y                    %% point
%  double(kind=8), intent(in)      :: h                    %% element spacing
%  double(kind=8), intent(out)     :: v                    %% point
%  double(kind=8)                  :: pi, k, delta                %% temporary variables
  
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



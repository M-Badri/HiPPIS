clearvars -except fileID
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs= 30;
for k=1:4  %% loop of functions
    figure
    if(k ==1)
     dd = load('mapping_data/data/Runge2DEps');
    elseif(k ==2)
     dd = load('mapping_data/data/Surface1Eps');
    elseif(k ==3)
     dd = load('mapping_data/data/Surface2Eps');
    elseif(k ==4)
     dd = load('mapping_data/data/Heaviside2DEps');
    end 

    npts= sqrt(length(dd(:,1)));
    xx= zeros(npts); 
    yy= xx;
    vv0 = xx;
    vv = xx;
    idx = 0;
    for jj=1:npts
      for ii=1:npts
        idx = idx + 1;
        xx(ii, jj) =dd(idx, 1); 
        yy(ii, jj) =dd(idx, 2); 
        vv0(ii, jj) = dd(idx, 3);
      end
    end 
 

    for i=1:6
      idx = 0;
      for jj=1:npts
        for ii=1:npts
          idx = idx + 1;
          vv(ii, jj) = dd(idx, i+3);
        end
      end 
      %if(k==2)
      %figure
      %surf(xx, yy, vv0)
      %pause
      %end
      if(i==1)
        subplot(1,2,1)
        surf(xx, yy, vv)%, 'FaceAlpha', 0.5)
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('$$\epsilon_{0}=1, \epsilon_{1}=1$$', 'Interpreter', 'latex', 'Fontsize', fs)
        if(k==2)
          zlim([0 1.1])
        elseif(k==4)
          zlim([0 1.25])
        end

      %elseif(i==3) 
      %  subplot(2,2,3)
      %  surf(xx, yy, vv)%, 'FaceAlpha', 0.5)
      elseif((i==6 && k~=4) || (i==6 && k==4)) 
        subplot(1,2,2)
        surf(xx, yy, vv)%, 'FaceAlpha', 0.5)
        xlabel('x')
        ylabel('y')
        zlabel('z')
        if(k==4)
          title('$$\epsilon_{0}=10^{-4}, \epsilon_{1}=10^{-4}$$', 'Interpreter', 'latex', 'Fontsize', fs)
        else
          title('$$\epsilon_{0}=10^{-4}, \epsilon_{1}=1$$', 'Interpreter', 'latex', 'Fontsize', fs)
        end
        if(k==2)
          zlim([0 1.1])
        elseif(k==4)
          zlim([0 1.25])
        end
      end 
      %ylim([0.8, 1.0])
    end
end 
pause
%% 2D tables used in the manuscript

st1 = ["Runge2D", "T1", "Heaviside2D"];
st2 = ["01", "04", "08"];
st3 = ["017x017", "033x033", "065x065", "129x129", "257x257"];
st4 = ["st=1", "st=2", "st=3"];

% Setup strings for files to be read
strDBI   = strings(5, 3);
strPPI   = strings(5, 3);
strPCHIP = strings(5, 1);

strDBI2   = strings(5, 3, 3);
strPPI2   = strings(5, 3, 3);
%
nx = [17; 33; 65; 129; 257];
ny = [17; 33; 65; 129; 257];
d = [1, 4, 8];
%

%
for k=1:3
  %% get file names to be read
  for i=1:5
    strPCHIP(i) = strcat("mapping_data/data/",st1(k), "PCHIP", "03", st3(i));
  end
  for j=1:3
    for i=1:5
      strDBI(i,j)   = strcat("mapping_data/data/",st1(k), "DBI",   "03", st3(i), st4(j) );
      strPPI(i,j)   = strcat("mapping_data/data/",st1(k), "PPI",   "03", st3(i), st4(j) );
    end
  end
  %
  %% variables to hold error calculations %%
  err1_pchip = zeros(length(nx), 3);
  err2_pchip = err1_pchip;
  err1_dbi = err1_pchip;
  err2_dbi = err1_pchip;
  err1_ppi = err1_pchip;
  err2_ppi = err1_pchip;
  %
  %% variables to hold the rate of convergence %%
  rate1_pchip = zeros(length(nx), 3);
  rate2_pchip = rate1_pchip;
  rate1_dbi = rate1_pchip;
  rate2_dbi = rate1_pchip;
  rate1_ppi = rate1_pchip;
  rate2_ppi = rate1_pchip;


  for j=1:3
    for i=1:5

      %%----------------------------------------
      %% PCHIP 
      %%----------------------------------------
      %% Load standard interpoaltion data 
      dd = load(char(strPCHIP(i, 1)));

      %% assumes that x and y have the same length 
      npts= sqrt(length(dd(:,1)));
      xx= zeros(npts); 
      yy= xx;
      zzt= xx;
      zz1_pchip= xx;
      zz2_pchip= xx; 
      idx = 0;
      for jj=1:npts
        for ii=1:npts
          idx = idx + 1;
          xx(ii, jj) =dd(idx, 1); 
          yy(ii, jj) =dd(idx, 2); 
          zzt(ii, jj) = dd(idx, 3);
          zz1_pchip(ii, jj) = dd(idx, 4);
          zz2_pchip(ii, jj) = dd(idx, 5);
        end
      end 
 
      %% calcualate errors
      err1_pchip(i,j) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz1_pchip).^2, 2)) );
      err2_pchip(i,j) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz2_pchip).^2, 2)) );


      %%----------------------------------------
      %% DBI 
      %%----------------------------------------
      %% Load DBI
      dd = load(char(strDBI(i, j)));

      %% assumes that x and y have the same length 
      npts= sqrt(length(dd(:,1)));
      xx= zeros(npts); 
      yy= xx;
      zzt= xx;
      zz1_dbi= xx;
      zz2_dbi= xx; 
      idx = 0;
      for jj=1:npts
        for ii=1:npts
          idx = idx + 1;
          xx(ii, jj) =dd(idx, 1); 
          yy(ii, jj) =dd(idx, 2); 
          zzt(ii, jj) = dd(idx, 3);
          zz1_dbi(ii, jj) = dd(idx, 4);
          zz2_dbi(ii, jj) = dd(idx, 5);
        end
      end 
 
      %% calcualate errors
      err1_dbi(i,j) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz1_dbi).^2, 2)) );
      err2_dbi(i,j) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz2_dbi).^2, 2)) );

      %%----------------------------------------
      %% PPI 
      %%----------------------------------------
      %% Load standard interpoaltion data 
      dd = load(char(strPPI(i, j)));

      %% assumes that x and y have the same length 
      npts= sqrt(length(dd(:,1)));
      xx= zeros(npts); 
      yy= xx;
      zzt= xx;
      zz1_ppi= xx;
      zz2_ppi= xx; 
      idx = 0;
      for jj=1:npts
        for ii=1:npts
          idx = idx + 1;
          xx(ii, jj) =dd(idx, 1); 
          yy(ii, jj) =dd(idx, 2); 
          zzt(ii, jj) = dd(idx, 3);
          zz1_ppi(ii, jj) = dd(idx, 4);
          zz2_ppi(ii, jj) = dd(idx, 5);
        end
      end 
 
      %% calcualate errors
      err1_ppi(i,j) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz1_ppi).^2, 2)) );
      err2_ppi(i,j) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz2_ppi).^2, 2)) );

      %% calculate convergence rates
      if(i > 1)
       rate1_pchip(i,j) = log( err1_pchip(i-1,j)/err1_pchip(i,j) ) / log( double(nx(i))/double(nx(i-1)) );
       rate2_pchip(i,j) = log( err2_pchip(i-1,j)/err2_pchip(i,j) ) / log( double(nx(i))/double(nx(i-1)) );
       rate1_dbi(i,j) = log( err1_dbi(i-1,j)/err1_dbi(i,j) ) / log( double(nx(i))/double(nx(i-1)) );
       rate2_dbi(i,j) = log( err2_dbi(i-1,j)/err2_dbi(i,j) ) / log( double(nx(i))/double(nx(i-1)) );
       rate1_ppi(i,j) = log( err1_ppi(i-1,j)/err1_ppi(i,j) ) / log( double(nx(i))/double(nx(i-1)) );
       rate2_ppi(i,j) = log( err2_ppi(i-1,j)/err2_ppi(i,j) ) / log( double(nx(i))/double(nx(i-1)) );
      end
 
    end
  end

  fprintf(fileID, '*****  2D Uniform mesh fun = %d ***** \n', k);
  for i=1:5
     fprintf(  fileID, '%d \t && %.2E  &&  %.2E  &  %.2E  &  %.2E  &&  %.2E  &  %.2E  &  %.2E   \\\\ \n', ...
               nx(i), err1_pchip(i,1), err1_dbi(i,1), err1_dbi(i,2), err1_dbi(i, 3), ...
                                       err1_ppi(i,1), err1_ppi(i,2), err1_ppi(i, 3) );
  end 


  %%% 
  for j=1:3
    for i=1:3
      for kk=1:5
        strDBI2(kk,i,j)   = strcat("mapping_data/data/",st1(k), "DBI", st2(j), st3(kk), st4(i) );
        strPPI2(kk,i,j)   = strcat("mapping_data/data/",st1(k), "PPI", st2(j), st3(kk), st4(i) );
      end
    end

    for i=1:3
      for kk=1:5
        %%----------------------------------------
        %% PCHIP 
        %%----------------------------------------
        %% Load DBI
        dd = load(char(strDBI2(kk,i,j)));

        %% assumes that x and y have the same length 
        npts= sqrt(length(dd(:,1)));
        xx= zeros(npts); 
        yy= xx;
        zzt= xx;
        zz1_dbi= xx;
        zz2_dbi= xx; 
        idx = 0;
        for jj=1:npts
          for ii=1:npts
            idx = idx + 1;
            xx(ii, jj) =dd(idx, 1); 
            yy(ii, jj) =dd(idx, 2); 
            zzt(ii, jj) = dd(idx, 3);
            zz1_dbi(ii, jj) = dd(idx, 4);
            zz2_dbi(ii, jj) = dd(idx, 5);
          end
        end 
 
        %% calcualate errors
        err1_dbi(kk,j) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz1_dbi).^2, 2)) );
        err2_dbi(kk,j) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz2_dbi).^2, 2)) );

        %%----------------------------------------
        %% PPI 
        %%----------------------------------------
        %% Load standard interpoaltion data 
        dd = load(char(strPPI2(kk, i, j)));

        %% assumes that x and y have the same length 
        npts= sqrt(length(dd(:,1)));
        xx= zeros(npts); 
        yy= xx;
        zzt= xx;
        zz1_ppi= xx;
        zz2_ppi= xx; 
        idx = 0;
        for jj=1:npts
          for ii=1:npts
            idx = idx + 1;
            xx(ii, jj) =dd(idx, 1); 
            yy(ii, jj) =dd(idx, 2); 
            zzt(ii, jj) = dd(idx, 3);
            zz1_ppi(ii, jj) = dd(idx, 4);
            zz2_ppi(ii, jj) = dd(idx, 5);
          end
        end 
 
        %% calcualate errors
        err1_ppi(kk, i) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz1_ppi).^2, 2)) );
        err2_ppi(kk, i) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz2_ppi).^2, 2)) );

        %% calculate convergence rates
        if(kk > 1)
         rate1_dbi(kk, i) = log( err1_dbi(kk-1,i)/err1_dbi(kk,i) ) / log( double(nx(kk))/double(nx(kk-1)) );
         rate2_dbi(kk, i) = log( err2_dbi(kk-1,i)/err2_dbi(kk,i) ) / log( double(nx(kk))/double(nx(kk-1)) );
         rate1_ppi(kk, i) = log( err1_ppi(kk-1,i)/err1_ppi(kk,i) ) / log( double(nx(kk))/double(nx(kk-1)) );
         rate2_ppi(kk, i) = log( err2_ppi(kk-1,i)/err2_ppi(kk,i) ) / log( double(nx(kk))/double(nx(kk-1)) );
        end
      end 
    end

    fprintf(fileID, '****** d= %d ****** \n', d(j) );
    for i=1:5
       fprintf(  fileID, '%d \t &&  %.2E  &  %.2E  &  %.2E  &&  %.2E  &  %.2E  &  %.2E   \\\\ \n', ...
                 nx(i),  err1_dbi(i,1), err1_dbi(i,2), err1_dbi(i, 3), ...
                         err1_ppi(i,1), err1_ppi(i,2), err1_ppi(i, 3) );
    end 
 
  end
  
end




if(exist('fileID'))
  clearvars -except fileID;
else
  clear;
  fileID = fopen('approximations_tables_1d_2d.txt', 'w');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fileID, '---------- Errors from 2D approximations examples ---------- \n');
fs= 30;
for k=1:3  %% loop of functions
    figure
    if(k ==1)
     dd = load('mapping_data/data/Runge2DEps');
    elseif(k ==2)
     dd = load('mapping_data/data/Heaviside2DEps');
    elseif(k ==3)
     dd = load('mapping_data/data/Surface1Eps');
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
          zlim([0 1.25])
        elseif(k==3)
          zlim([0 1.1])
        end

      %elseif(i==3) 
      %  subplot(2,2,3)
      %  surf(xx, yy, vv)%, 'FaceAlpha', 0.5)
      elseif((i==6 && k~=2) || (i==6 && k==2)) 
        subplot(1,2,2)
        surf(xx, yy, vv)%, 'FaceAlpha', 0.5)
        xlabel('x')
        ylabel('y')
        zlabel('z')
        if(k==2)
          title('$$\epsilon_{0}=10^{-4}, \epsilon_{1}=10^{-4}$$', 'Interpreter', 'latex', 'Fontsize', fs)
        else
          title('$$\epsilon_{0}=10^{-4}, \epsilon_{1}=1$$', 'Interpreter', 'latex', 'Fontsize', fs)
        end
        if(k==2)
          zlim([0 1.25])
        elseif(k==3)
          zlim([0 1.1])
        end
      end 
      %ylim([0.8, 1.0])
    end
end 
%pause
% 2D tables used in the manuscript

st1 = ["Runge2D", "Heaviside2D", "T1" ];
st2 = ["03", "04", "08"];
st3 = ["017x017", "033x033", "065x065", "129x129", "257x257"];
st4 = "st=3"; %["st=1", "st=2", "st=3"];

% Setup strings for files to be read
strDBI   = strings(5, 1);
strPPI   = strings(5, 1);
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
      strDBI(i,j)   = strcat("mapping_data/data/",st1(k), "DBI",   st2(j), st3(i), st4 );
      strPPI(i,j)   = strcat("mapping_data/data/",st1(k), "PPI",   st2(j), st3(i), st4 );
    end
  end
  %
  %% variables to hold error calculations %%
  err1_pchip = zeros(length(nx), 1);
  err1_dbi = zeros(length(nx), 3);
  err1_ppi = err1_dbi;
  %

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
    idx = 0;
    for jj=1:npts
      for ii=1:npts
        idx = idx + 1;
        xx(ii, jj) =dd(idx, 1); 
        yy(ii, jj) =dd(idx, 2); 
        zzt(ii, jj) = dd(idx, 3);
        zz1_pchip(ii, jj) = dd(idx, 4);
      end
    end 
 
    %% calcualate errors
    err1_pchip(i,1) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz1_pchip).^2, 2)) );
  end

  for j=1:3
    for i=1:5

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
      idx = 0;
      for jj=1:npts
        for ii=1:npts
          idx = idx + 1;
          xx(ii, jj) =dd(idx, 1); 
          yy(ii, jj) =dd(idx, 2); 
          zzt(ii, jj) = dd(idx, 3);
          zz1_dbi(ii, jj) = dd(idx, 4);
        end
      end 
 
      %% calcualate errors
      err1_dbi(i,j) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz1_dbi).^2, 2)) );

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
      idx = 0;
      for jj=1:npts
        for ii=1:npts
          idx = idx + 1;
          xx(ii, jj) =dd(idx, 1); 
          yy(ii, jj) =dd(idx, 2); 
          zzt(ii, jj) = dd(idx, 3);
          zz1_ppi(ii, jj) = dd(idx, 4);
        end
      end 
 
      %% calcualate errors
      err1_ppi(i,j) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz1_ppi).^2, 2)) );
    end
  end

  fprintf(fileID, '*****  2D Uniform mesh fun = %d ***** \n', k);
  for i=1:5
     fprintf(  fileID, '%d \t && %.2E  &&  %.2E  &  %.2E  &  %.2E  &&  %.2E  &  %.2E  &  %.2E   \\\\ \n', ...
               nx(i), err1_pchip(i,1), err1_dbi(i,1), err1_dbi(i,2), err1_dbi(i, 3), ...
                                       err1_ppi(i,1), err1_ppi(i,2), err1_ppi(i, 3) );
  end 
end



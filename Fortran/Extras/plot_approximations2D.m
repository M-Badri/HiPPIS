addpath("../../Extras");
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
fs= 18;
ss = 0.1;
for k=2:3  % k=1 2D version of Runge, k=2 2D version of Heaviside, k=3 surface functions.

    if(k ==1)
     dd4 = load('mapping_data/data/Runge2DEps_4'); % Runge example with d=4 
     dd8 = load('mapping_data/data/Runge2DEps_8'); % Runge examples with d=8
     lim_right = 1.0; 
     az = -37.50;
     el = 30.0;
    elseif(k ==2)
     dd4 = load('mapping_data/data/Heaviside2DEps_4');
     dd8 = load('mapping_data/data/Heaviside2DEps_8');
     lim_right = 1.25;
     az = 43.33;
     el = 5.60;
    elseif(k ==3)
     dd4 = load('mapping_data/data/Surface1Eps_4');
     dd8 = load('mapping_data/data/Surface1Eps_8');
     lim_right = 1.1;
     az = -25.82;
     el = 2.39;
    end
    %
    npts= sqrt(length(dd4(:,1)));
    xx= zeros(npts); 
    yy= xx;
    vv0 = xx;
    vv = xx;
    idx = 0;
    for jj=1:npts
      for ii=1:npts
        idx = idx + 1;
        xx(ii, jj) =dd4(idx, 1); 
        yy(ii, jj) =dd4(idx, 2); 
        vv0(ii, jj) = dd4(idx, 3);
      end
    end 
    idx = 0;
    for jj=1:npts
      for ii=1:npts
        idx = idx + 1;
        vv4_1(ii, jj) = dd4(idx, 4);
        vv4_2(ii, jj) = dd4(idx, 9);
        vv8_1(ii, jj) = dd8(idx, 4);
        vv8_2(ii, jj) = dd8(idx, 9);
        vvpchip(ii, jj) = dd4(idx, 11);
        vvmqsi(ii, jj) = dd4(idx, 12);
      end
    end 
    %subplot(3,2,1)
    figure
    s = surf(xx, yy, vvpchip); %, 'FaceAlpha', 0.5);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    zlim([0 lim_right])
    view(az, el)
    s.EdgeColor= 'none';
    set(gca, 'FontSize', fs)
    %title('PCHIP', 'Interpreter', 'latex', 'Fontsize', fs)
    %
    %subplot(3,2,2)
    figure
    s2=surf(xx, yy, vvmqsi);%, 'FaceAlpha', 0.5);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    zlim([0 lim_right])
    view(az, el)
    s2.EdgeColor= 'none';
    set(gca, 'FontSize', fs)
    %title('MQSI', 'Interpreter', 'latex', 'Fontsize', fs)
    %
    %subplot(3,2,3)
    figure
    s3=surf(xx, yy, vv4_1);%, 'FaceAlpha', 0.5);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    zlim([0 lim_right])
    view(az, el)
    s3.EdgeColor= 'none';
    set(gca, 'FontSize', fs)
    %title('PPI $$\mathcal{P}_{4} \epsilon_{0}=1, \epsilon_{1}=1$$', 'Interpreter', 'latex', 'Fontsize', fs)
    %
    %subplot(3,2,4)
    figure
    s4=surf(xx, yy, vv4_2);%, 'FaceAlpha', 0.5);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    zlim([0 lim_right])
    view(az, el)
    s4.EdgeColor= 'none';
    set(gca, 'FontSize', fs)
    %title('$$PPI \mathcal{P}_{4} \epsilon_{0}=10^{-4}, \epsilon_{1}=10^{-4}$$', 'Interpreter', 'latex', 'Fontsize', fs)
    %
    %subplot(3,2,5)
    figure
    s5= surf(xx, yy, vv8_1);%, 'FaceAlpha', 0.5);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    zlim([0 lim_right])
    view(az, el)
    s5.EdgeColor= 'none';
    set(gca, 'FontSize', fs)
    %title('$$PPI \mathcal{P}_{8} \epsilon_{0}=1, \epsilon_{1}=1$$', 'Interpreter', 'latex', 'Fontsize', fs)
    %
    %subplot(3,2,6)
    figure
    s6=surf(xx, yy, vv8_2); %, 'FaceAlpha', 0.5);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    zlim([0 lim_right])
    view(az, el)
    s6.EdgeColor= 'none';
    set(gca, 'FontSize', fs)
    %title('$$PPI \mathcal{P}_{8} \epsilon_{0}=10^{-4}, \epsilon_{1}=10^{-4}$$', 'Interpreter', 'latex', 'Fontsize', fs)
end

pause
% 2D tables used in the manuscript

st1 = ["Runge2D", "Heaviside2D", "T1" ];
st2 = ["03", "04", "08"];
st3 = ["017x017", "033x033", "065x065", "129x129", "257x257"];
st4 = "st=3"; %["st=1", "st=2", "st=3"];

% Setup strings for files to be read
strDBI   = strings(5, 1);
strPPI   = strings(5, 1);
strPCHIP = strings(5, 1);
strMQSI = strings(5, 1);

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
    strMQSI(i) = strcat("mapping_data/data/",st1(k), "MQSI", "05", st3(i));
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
  err1_mqsi = zeros(length(nx), 1);
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

    %%----------------------------------------
    %% MQSI 
    %%----------------------------------------
    %% Load standard interpoaltion data 
    dd = load(char(strMQSI(i, 1)));

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
        zz1_mqsi(ii, jj) = dd(idx, 4);
      end
    end 
 
    %% calcualate errors
    err1_mqsi(i,1) = sqrt( trapz( yy(1, :), trapz(xx(:,1), (zzt-zz1_mqsi).^2, 2)) );
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
     fprintf(  fileID, '%d \t && %.2E  && %.2E  &&  %.2E  &  %.2E  &  %.2E  &&  %.2E  &  %.2E  &  %.2E   \\\\ \n', ...
               nx(i), err1_pchip(i,1),err1_mqsi(i,1), err1_dbi(i,1), err1_dbi(i,2), err1_dbi(i, 3), ...
                                       err1_ppi(i,1), err1_ppi(i,2), err1_ppi(i, 3) );
  end 
end




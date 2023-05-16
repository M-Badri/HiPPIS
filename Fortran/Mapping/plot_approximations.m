if(exist('fileID'))
  clearvars -except fileID;
else
  clear;
  fileID = fopen('approximations_tables_1d_2d.txt', 'w');
end
%close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1D Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fileID, '---------- Errors from 1D approximations examples ---------- \n');

c1 = [0       0.4470  0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250];
c4 = [0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880];
c6 = [0.3010, 0.7450, 0.9330];
c7 = [0.6350, 0.0780, 0.1840];
fs = 28;
ls = 6;

% 1D plots used for figures in manuscript
for i=1:3
  if(i == 1)
    dd4 = load('mapping_data/data/RungeEps_4');
    dd8 = load('mapping_data/data/RungeEps_8');
  elseif(i == 2)
    dd4 = load('mapping_data/data/HeavisideEps_4');
    dd8 = load('mapping_data/data/HeavisideEps_8');
  elseif(i == 3)
    dd4 = load('mapping_data/data/GelbTEps_4');
    dd8 = load('mapping_data/data/GelbTEps_8');
  end
  figure
  subplot(1,2,1)
  if(i ==3)
    dd4(10000/4-1:10000/4+1, 2) = NaN;
  end
  plot( dd4(:,1), dd4(:,2), 'k', 'LineWidth', ls) 
  hold on
  plot( dd4(:,1), dd4(:,10), 'Color', c1, 'LineWidth', ls)
  plot( dd4(:,1), dd4(:,11), 'Color', c2, 'LineWidth', ls)
  plot( dd4(:,1), dd4(:,9),  'Color', c3, 'LineWidth', ls)
                     %'LineWidth', 4)
  legend('True', 'PCHIP', 'MQS', 'DBI')
  xlabel('x')
  ylabel('y')
  set(gca, 'FontSize', fs)
  subplot(1,2,2)
  plot( dd4(:,1), dd4(:,2), 'k', 'LineWidth', ls) 
  hold on
  plot( dd4(:,1), dd4(:,10), 'Color', c1, 'LineWidth', ls)
  plot( dd4(:,1), dd4(:,11), 'Color', c2, 'LineWidth', ls)
  plot( dd4(:,1), dd4(:,9),  'Color', c3, 'LineWidth', ls)
  if(i==1)
    xlim([0.2, 0.65])
    %ylim([0.95 1.05])
  elseif(i==2)
    xlim([0.02, 0.2])
    ylim([0.95 1.05])
  elseif(i==3)
    xlim([-1, -0.80])
    %ylim([0.95 1.05])
  end
  xlabel('x')
  ylabel('y')
  %ylim([0.03, 0.06])
  set(gca, 'FontSize', fs)

  figure
  subplot(1,2,1)
  if(i ==3)
    dd4(10000/4-1:10000/4+1, 2) = NaN;
  end
  plot( dd4(:,1), dd4(:,2), 'k', 'LineWidth', ls) 
  hold on
  plot( dd4(:,1), dd4(:,3), 'Color', c4, 'LineWidth', ls)
  plot( dd4(:,1), dd4(:,5), 'Color', c5, 'LineWidth', ls)
  plot( dd8(:,1), dd8(:,3), 'Color', c6, 'LineWidth', ls)
  plot( dd8(:,1), dd8(:,5), 'Color', c7, 'LineWidth', ls)
  legend('True', '$$ \mathcal{P}_{4}, \epsilon_{0}=1.0$$', '$$\mathcal{P}_{4}, \epsilon_{0}=0.01$$', '$$ \mathcal{P}_{8}, \epsilon_{0}=1.0$$', '$$\mathcal{P}_{8}, \epsilon_{0}=0.01$$', 'Interpreter', 'latex')
  xlabel('x')
  ylabel('y')
  set(gca, 'FontSize', fs)
  subplot(1,2,2)
  plot( dd4(:,1), dd4(:,2), 'k', 'LineWidth', ls)
  hold on
  plot( dd4(:,1), dd4(:,3), 'Color', c4, 'LineWidth', ls)
  plot( dd4(:,1), dd4(:,5), 'Color', c5, 'LineWidth', ls)
  plot( dd8(:,1), dd8(:,3), 'Color', c6, 'LineWidth', ls)
  plot( dd8(:,1), dd8(:,5), 'Color', c7, 'LineWidth', ls)
  if(i==1)
    xlim([0.2, 0.65])
    %ylim([0.95 1.05])
  elseif(i==2)
    xlim([0.02, 0.2])
    ylim([0.95 1.05])
  elseif(i==3)
    xlim([-1, -0.80])
    %ylim([0.95 1.05])
  end
  xlabel('x')
  ylabel('y')
  %ylim([0.03, 0.06])
  set(gca, 'FontSize', fs)
  pause

end
pause
%% 1D Tables used in manscript %%
%pause
%% string to be used for file names %%
st1 = ["Runge", "Heaviside", "GelbT"];
st2 = ["03", "04", "08"];
st3 = ["017", "033", "065", "129", "257"];
st4 = "st=3"; %["st=1", "st=2", "st=3"];

% Setup strings for files to be read
strDBI   = strings(5, 1);
strPPI   = strings(5, 1);
strPCHIP = strings(5, 1);
strMQSI = strings(5, 1);
%
strDBI2   = strings(5, 3, 3);
strPPI2   = strings(5, 3, 3);
%
n = [17, 33, 65, 129, 257];	%% number of input points
d = [1, 4, 8];
%
%xplt = zeros(10000,1);
%yplt = zeros(10000,4);

for k=1:3
  %% get file names to be read
  for i=1:5
    strPCHIP(i) = strcat("mapping_data/data/",st1(k), "PCHIP", st2(1), st3(i));
    strMQSI(i) = strcat("mapping_data/data/",st1(k), "MQSI", st2(1), st3(i));
  end
  for j=1:3 
    for i=1:5
      strDBI(i,j)   = strcat("mapping_data/data/",st1(k), "DBI",   st2(j), st3(i), st4 );
      strPPI(i,j)   = strcat("mapping_data/data/",st1(k), "PPI",   st2(j), st3(i), st4 );
    end
  end
  %
 
  %% variables to hold error calculations %%
  err1_pchip = zeros(length(n), 1);
  err1_mqsi = zeros(length(n), 1);
  err1_dbi = zeros(length(n), 3);
  err1_dbi = err1_dbi;
  %
  
  for i=1:5
    %% Load PCHIP data
    dd = load(char(strPCHIP(i, 1)));
    x = dd(:,1);
    yt = dd(:, 2);
    y1_pchip = dd(:,3);
    err1_pchip(i,1) = sqrt( trapz(x, (yt-y1_pchip).^2) );

    %% Load MQSI data
    dd = load(char(strMQSI(i, 1)));
    x = dd(:,1);
    yt = dd(:, 2);
    y1_mqsi = dd(:,3);
    err1_mqsi(i,1) = sqrt( trapz(x, (yt-y1_mqsi).^2) );
 
  end
  %
  for j=1:3
    for i=1:5
      %% Load DBI
      dd = load(char(strDBI(i, j)));
      x = dd(:,1);
      yt = dd(:, 2);
      y1_dbi = dd(:,3);

      %% Load PPI
      dd = load(char(strPPI(i, j)));
      y1_ppi = dd(:,3);
 
      err1_dbi(i,j) = sqrt( trapz(x, (yt-y1_dbi).^2) );
      err1_ppi(i,j) = sqrt( trapz(x, (yt-y1_ppi).^2) );
    end
  end

  fprintf(fileID, '*****  Uniform mesh fun = %d ***** \n', k);
  for i=1:5
     fprintf(  fileID, '%d \t && %.2E  && %.2E &&  %.2E  &  %.2E  &  %.2E  &&  %.2E  &  %.2E  &  %.2E   \\\\ \n', ...
               n(i), err1_pchip(i,1), err1_mqsi(i,1), err1_dbi(i,1), err1_dbi(i,2), err1_dbi(i, 3), ...
                                      err1_ppi(i,1), err1_ppi(i,2), err1_ppi(i, 3) );
  end 

end 


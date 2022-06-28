clearvars -except fileID
close all;
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1D Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1D plots used for figures in manuscript
for i=1:3
  figure
  if(i == 1)
    dd = load('mapping_data/data/RungeEps');
  elseif(i ==2)
    dd = load('mapping_data/data/HeavisideEps');
  elseif(i ==3)
    dd = load('mapping_data/data/GelbTEps');
  end
  subplot(1,2,1)
  plot( dd(:,1), dd(:,2), 'k', ...
        dd(:,1), dd(:,9),  ...
        dd(:,1), dd(:,3),  ...
        dd(:,1), dd(:,4),  ...
        dd(:,1), dd(:,5),  ...
        dd(:,1), dd(:,8),  ...
                     'LineWidth', 4)
  legend('True', 'DBI', 'eps=1.0', 'eps=0.1', 'eps=0.01', 'eps=0.0')
  xlabel('x')
  ylabel('y')
  set(gca, 'FontSize', 32)
  subplot(1,2,2)
  plot( dd(:,1), dd(:,2), 'k', ...
        dd(:,1), dd(:,9),  ...
        dd(:,1), dd(:,3),  ...
        dd(:,1), dd(:,4),  ...
        dd(:,1), dd(:,5),  ...
        dd(:,1), dd(:,8),  ...
                     'LineWidth', 4)
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
  set(gca, 'FontSize', 32)
end

%% 1D Tables used in manscript %%

%% string to be used for file names %%
st1 = ["Runge", "Heaviside", "GelbT"];
st2 = ["01", "04", "08"];
st3 = ["017", "033", "065", "129", "257"];
st4 = ["st=1", "st=2", "st=3"];

% Setup strings for files to be read
strDBI   = strings(5, 3);
strPPI   = strings(5, 3);
strPCHIP = strings(5, 1);
%
strDBI2   = strings(5, 3, 3);
strPPI2   = strings(5, 3, 3);
%
n = [17, 33, 65, 129, 257];	%% number of input points
d = [1, 4, 8];
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
  err1_pchip = zeros(length(n), 3);
  err2_pchip = err1_pchip;
  err1_dbi = err1_pchip;
  err2_dbi = err1_pchip;
  err1_ppi = err1_pchip;
  err2_ppi = err1_pchip;
  %
  %% variables to hold the rate of convergence %%
  rate1_pchip = zeros(length(n), 3);
  rate2_pchip = rate1_pchip;
  rate1_dbi = rate1_pchip;
  rate2_dbi = rate1_pchip;
  rate1_ppi = rate1_pchip;
  rate2_ppi = rate1_pchip;

  for j=1:3
    for i=1:5
      %% Load DBI
      dd = load(char(strDBI(i, j)));
      x = dd(:,1);
      yt = dd(:, 2);
      y1_dbi = dd(:,3);
      y2_dbi = dd(:, 4);

      %% Load PPI
      dd = load(char(strPPI(i, j)));
      y1_ppi = dd(:,3);
      y2_ppi = dd(:, 4);
 
       %% Load PPI
      dd = load(char(strPCHIP(i, 1)));
      y1_pchip = dd(:,3);
      y2_pchip = dd(:, 4);

      %% calculate errors
      err1_pchip(i,j) = sqrt( trapz(x, (yt-y1_pchip).^2) );
      err2_pchip(i,j) = sqrt( trapz(x, (yt-y2_pchip).^2) );
      err1_dbi(i,j) = sqrt( trapz(x, (yt-y1_dbi).^2) );
      err2_dbi(i,j) = sqrt( trapz(x, (yt-y2_dbi).^2) );
      err1_ppi(i,j) = sqrt( trapz(x, (yt-y1_ppi).^2) );
      err2_ppi(i,j) = sqrt( trapz(x, (yt-y2_ppi).^2) );
   
      %% calculate convergence rates
      if(i > 1)
       rate1_pchip(i,j) = log( err1_pchip(i-1)/err1_pchip(i) ) / log( double(n(i))/double(n(i-1)) );
       rate2_pchip(i,j) = log( err2_pchip(i-1)/err2_pchip(i) ) / log( double(n(i))/double(n(i-1)) );
       rate1_dbi(i,j) = log( err1_dbi(i-1)/err1_dbi(i) ) / log( double(n(i))/double(n(i-1)) );
       rate2_dbi(i,j) = log( err2_dbi(i-1)/err2_dbi(i) ) / log( double(n(i))/double(n(i-1)) );
       rate1_ppi(i,j) = log( err1_ppi(i-1)/err1_ppi(i) ) / log( double(n(i))/double(n(i-1)) );
       rate2_ppi(i,j) = log( err2_ppi(i-1)/err2_ppi(i) ) / log( double(n(i))/double(n(i-1)) );
      end
    end
  end

  fprintf(fileID, '*****  Uniform mesh fun = %d ***** \n', k);
  for i=1:5
     fprintf(  fileID, '%d \t && %.2E  &&  %.2E  &  %.2E  &  %.2E  &&  %.2E  &  %.2E  &  %.2E   \\\\ \n', ...
               n(i), err1_pchip(i,1), err1_dbi(i,1), err1_dbi(i,2), err1_dbi(i, 3), ...
                                      err1_ppi(i,1), err1_ppi(i,2), err1_ppi(i, 3) );
  end 

  %%%
  %%fprintf(fileID, '*****  LGL mesh fun = %d ***** \n', k);
  %%for i=1:5
  %%   fprintf(  fileID, '%d \t && %.2E  &&  %.2E  &  %.2E  &  %.2E  &&  %.2E  &  %.2E  &  %.2E   \\\\ \n', ...
  %%             n(i), err1_pchip(i,1), err1_dbi(i,1), err1_dbi(i,2), err1_dbi(i, 3), ...
  %%                                    err1_dbi(i,1), err1_ppi(i,2), err1_ppi(i, 3) );
  %%end 
  
  for j=1:3
    for i=1:3
      for ii=1:5
        strDBI2(ii,i,j)   = strcat("mapping_data/data/",st1(k), "DBI", st2(j), st3(ii), st4(i) );
        strPPI2(ii,i,j)   = strcat("mapping_data/data/",st1(k), "PPI", st2(j), st3(ii), st4(i) );
      end
    end

    for i=1:3
      for ii=1:5
        %% Load DBI
        dd = load(char(strDBI2(ii,i,j)));
        x = dd(:,1);
        yt = dd(:, 2);
        y1_dbi = dd(:,3);
        y2_dbi = dd(:, 4);

        %% Load PPI
        dd = load(char(strPPI2(ii,i,j)));
        y1_ppi = dd(:,3);
        y2_ppi = dd(:, 4);

        %% calculate errors
        err1_dbi(ii,i) = sqrt( trapz(x, (yt-y1_dbi).^2) );
        err2_dbi(ii,i) = sqrt( trapz(x, (yt-y2_dbi).^2) );
        err1_ppi(ii,i) = sqrt( trapz(x, (yt-y1_ppi).^2) );
        err2_ppi(ii,i) = sqrt( trapz(x, (yt-y2_ppi).^2) );
   
        %% calculate convergence rates
        if(i > 1)
         rate1_dbi(ii,i) = log( err1_dbi(i-1)/err1_dbi(i) ) / log( double(n(i))/double(n(i-1)) );
         rate2_dbi(ii,i) = log( err2_dbi(i-1)/err2_dbi(i) ) / log( double(n(i))/double(n(i-1)) );
         rate1_ppi(ii,i) = log( err1_ppi(i-1)/err1_ppi(i) ) / log( double(n(i))/double(n(i-1)) );
         rate2_ppi(ii,i) = log( err2_ppi(i-1)/err2_ppi(i) ) / log( double(n(i))/double(n(i-1)) );
        end
      end
    end    
    fprintf(fileID, '****** d= %d ****** \n', d(j) );
    for i=1:5
       fprintf(  fileID, '%d \t &&  %.2E  &  %.2E  &  %.2E  &&  %.2E  &  %.2E  &  %.2E   \\\\ \n', ...
                 n(i),  err1_dbi(i,1), err1_dbi(i,2), err1_dbi(i, 3), ...
                        err1_ppi(i,1), err1_ppi(i,2), err1_ppi(i, 3) );
    end 
  end
  
  
end 


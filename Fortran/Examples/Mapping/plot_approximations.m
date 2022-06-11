clear;
close all;
clc

fileID = fopen('approximations_tables_1d_2d.txt', 'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D third order results

st1 = ["Runge", "Heaviside", "GelbT"];
st2 = ["01", "04", "08"];
st3 = ["017", "033", "065", "129", "257"];
st4 = ["st=1", "st=2", "st=3"];

% Setup strings for files to be read
strDBI   = strings(5, 3);
strPPI   = strings(5, 3);
strPCHIP = strings(5, 1);

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
        char(strDBI2(ii,i,j))
pause
        x = dd(:,1);
        yt = dd(:, 2);
        y1_dbi = dd(:,3);
        y2_dbi = dd(:, 4);

        %% Load PPI
        dd = load(char(strPPI2(ii,i,j)));
        y1_ppi = dd(:,3);
        y2_ppi = dd(:, 4);

        %% calculate errors
        err1_dbi(i,j) = sqrt( trapz(x, (yt-y1_dbi).^2) );
        err2_dbi(i,j) = sqrt( trapz(x, (yt-y2_dbi).^2) );
        err1_ppi(i,j) = sqrt( trapz(x, (yt-y1_ppi).^2) );
        err2_ppi(i,j) = sqrt( trapz(x, (yt-y2_ppi).^2) );
   
        %% calculate convergence rates
        if(i > 1)
         rate1_dbi(i,j) = log( err1_dbi(i-1)/err1_dbi(i) ) / log( double(n(i))/double(n(i-1)) );
         rate2_dbi(i,j) = log( err2_dbi(i-1)/err2_dbi(i) ) / log( double(n(i))/double(n(i-1)) );
         rate1_ppi(i,j) = log( err1_ppi(i-1)/err1_ppi(i) ) / log( double(n(i))/double(n(i-1)) );
         rate2_ppi(i,j) = log( err2_ppi(i-1)/err2_ppi(i) ) / log( double(n(i))/double(n(i-1)) );
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

fclose(fileID);



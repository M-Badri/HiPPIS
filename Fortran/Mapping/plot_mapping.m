%
% Script to plot the mapping examples based on the modified Runge and 
% TWP_ICE examples.
%
if(exist('fileID'))
  clearvars -except fileID;
else
  clear;
  fileID = fopen('approximations_tables_1d_2d.txt', 'w');
end


fprintf(fileID, '---------- Errors from mapping examples ---------- \n');
%% String for data files to be read 
st5 = "st=3"; %["st=1", "st=2", "st=3"];
st4 = ["PPI"; "DBI"; "PCHIP"];
st3 = ["064";, "127"; "253"];
st2 = ["03"; "05"; "07"];
st1 = ["qc"; "runge"];

%% Set up files names to be used to create plots 
strPPI = strings(3,3);	% PPI results files
strDBI = strings(3,3);        % DBI results files
strPCHIP = strings(3,1);        % DBI results files

err_pchip= zeros(3,1,1);
err_ppi= zeros(3,3);
err_dbi= err_ppi;
%

d = [3, 5, 7];
n = [64; 127; 253];
m_size = 40;
lw = 3;
f_size = 24;
f_size2 = 24;
fig_idx = 0;
 
for ii=1:2 % loop over profiles

  for i=1:3 % loop over degree
    for j=1:3 % loop over number of points
      strPPI(i,j) = strcat( st1(ii), "d", st4(1), st2(i), st3(j), st5);
      strDBI(i,j) = strcat( st1(ii), "d", st4(2), st2(i), st3(j), st5);
    end
  end
  for j=1:3 % loop of st
    strPCHIP(j,1) = strcat(st1(ii), "dPCHIP03", st3(j));
  end

  for j=1:3
    dd_pchip = load( char(strcat("mapping_data/data/",strPCHIP(j,1,1))) );
    for i=1:3
      dd_ppi = load( char(strcat("mapping_data/data/",strPPI(i,j))) );
      dd_dbi = load( char(strcat("mapping_data/data/",strDBI(i,j))) );
      err_ppi(i,j) = max(abs(dd_ppi(:,3)-dd_ppi(:,2)));
      err_dbi(i,j) = max(abs(dd_dbi(:,3)-dd_ppi(:,2)));
    end
    err_pchip(j,1,1) = max(abs(dd_pchip(:,3)-dd_ppi(:,2)));
  end

  fprintf(fileID, ' **** fun  =  %d  **** \n', ii );
  for j=1:3
     fprintf(  fileID, '%d \t && %.2E  &&  %.2E  &  %.2E  &  %.2E  &&  %.2E  &  %.2E  &  %.2E   \\\\ \n', ...
               n(j), err_pchip(j,1), err_dbi(1,j), err_dbi(2,j), err_dbi(3,j), ...
                                        err_ppi(1,j), err_ppi(2,j), err_ppi(3,j) );
  end

end % loop over profiles 


%
% Script to plot the mapping examples based on the modifed Runge and 
% TWP_ICE examples.
%
clearvars -except fileID
close all
clc


%% String for data files to be read 
st5 = ["st=1", "st=2", "st=3"];
st4 = ["PPI"; "DBI"; "PCHIP"];
st3 = ["064";, "127"; "253"];
st2 = ["03"; "05"; "07"];
st1 = ["qc"; "runge";];

%% Set up files names to be used to create plots 
strPPI = strings(3,2,3);	% PPI results files
strDBI = strings(3,2,3);        % DBI results files
strDEGDBI = strings(3,2,3);     % polynomial degree used
strDEGPPI = strings(3,2,3);     % polynomial degree used
%

for ii=3:3

  fprintf('The parameter st =  %d  \n', ii );
  for i=1:3
    for j=1:2
      for k=1:3
        strPPI(i,j,k) = strcat( st1(j), "d", st4(1), st2(i), st3(k), st5(ii) );
        strDBI(i,j,k) = strcat( st1(j), "d", st4(2), st2(i), st3(k), st5(ii) );
        strDEGPPI(i,j,k) = strcat( st1(j), "dDEG", st4(1), st2(i), st3(k), st5(ii) );
        strDEGDBI(i,j,k) = strcat( st1(j), "dDEG", st4(2), st2(i), st3(k), st5(ii) );
      end
    end
  end
  
  strPCHIP = strings(1, 2, 3);  % PCHIP files
  %
  for j=1:2
    for k=1:3
      strPCHIP(1, j, k) = strcat(st1(j), "dPCHIP03", st3(k) );
    end
  end
  
  d = [3, 5, 7];
  m_size = 40;
  lw = 3;
  f_size = 24;
  f_size2 = 24;
  fig_idx = 0;
  
  for k=1:2
   for j=3:3
      figure
      for i=1:3
        dd_ppi = load( char(strcat("mapping_data/data/",strPPI(i,k, j))) );
        deg_dd_ppi = load( char(strcat("mapping_data/data/",strDEGPPI(i,k,j))) );
        dd_dbi = load( char(strcat("mapping_data/data/",strDBI(i,k,j))) );
        deg_dd_dbi = load( char(strcat("mapping_data/data/",strDEGDBI(i,k,j))) );
        dd_pchip = load( char(strcat("mapping_data/data/",strPCHIP(1,k,j))) );
%  char(strcat("mapping_data/data/",strPPI(i,k, j)))
%  pause
        subplot(3, 1, i);
        plot( dd_ppi(:,1), abs( (dd_pchip(:,3)-dd_ppi(:,2))), ':.', ...
              dd_ppi(:,1), abs( (dd_dbi(:,3)  -dd_ppi(:,2))),   ':.', ...
              dd_ppi(:,1), abs( (dd_ppi(:,3)  -dd_ppi(:,2))),   ':.', ...
              'MarkerSize', m_size, 'LineWidth', lw);
        legend('|PCHIP-Target|', '|DBI-Target|', '|PPI-Target|');
        if(k == 1)
          xlim([-0.9 -0.2])
        elseif(k ==2)
          xlim([-0.4, 0.4]);
        end
        %xlim([-0.90,-0.25])
        xlabel('x')
        ylabel('error')
        set(gca, 'FontSize', f_size);
        
      end
    end
  end 
end

%%
%% Script to plot the mapping examples based on the modifed Runge and 
%% TWP_ICE examples.

clear
close all
clc


%% String for data files to be read 
st5 = ["st=1", "st=2", "st=3"]
st4 = ["PPI"; "DBI"; "PCHIP"];
st3 = ["64";, "127"; "253"];
st2 = ["03"; "05"; "07"];
st1 = ["qc"; "runge";];

%% Set up files names to be used to create plots 
strPPI = strings(3,2,3);	% PPI results files
strDBI = strings(3,2,3);        % DBI results files
strDEGDBI = strings(3,2,3);     % polynomial degree used
strDEGPPI = strings(3,2,3);     % polynomial degree used
%
for i=1:3
  for j=1:2
    for k=1:3
      strPPI(i,j,k) = strcat( st1(j), "d", st4(1), st2(i), st3(k) );
      strDBI(i,j,k) = strcat( st1(j), "d", st4(2), st2(i), st3(k) );
      strDEGPPI(i,j,k) = strcat( st1(j), "dDEG", st4(1), st2(i), st3(k) );
      strDEGDBI(i,j,k) = strcat( st1(j), "dDEG", st4(2), st2(i), st3(k) );
    end
  end
end

strPCHIP = strings(1, 2, 3);  % PCHIP files
%
for j=1:2
  for k=1:3
    strPCHIP(1, i, k) = strcat(st1(i), "dPCHIP03", st3(k) );
  end
end

d = [3, 5, 6, 7];
m_size = 40;
lw = 3;
f_size = 24;
f_size2 = 24;
fig_idx = 0;

for k=1:2
 for j=1:3
    figure
    for i=1:3
      dd_ppi = load( char(strcat("data/",strPPI(i,k, j))) );
      deg_dd_ppi = load( char(strcat("data/",strDEGPPI(j,k,1))) );
      dd_dbi = load( char(strcat("data/",strDBI(j,k,1))) );
      deg_dd_dbi = load( char(strcat("data/",strDEGDBI(j,k,1))) );
      dd_pchip = load( char(strcat("data/",strPCHIP(1,k,1))) );

      subplot(3, 1, i);
      plot( zd, abs( (dd_pchip(:,3)-dd_ppi(:,1))), ':.', ...
            zd, abs( (dd_dbi(:,3)-dd_ppi(:,1))),   ':.', ...
            zd, abs( (dd_ppi(:,3)-dd_ppi(:,1))),   ':.', ...
            'MarkerSize', m_size, 'LineWidth', lw);
      legend('|PCHIP-Target|', '|PPI-Target|', '|PPI-Target|');
      if(k == 1)
        xlim([-0.4, 0.4]);
      elseif(k ==2)
        xlim([-0.9 -0.2])
      end
      %xlim([-0.90,-0.25])
      xlabel('x')
      ylabel('error')
      set(gca, 'FontSize', f_size);
      
    end
  end
end 


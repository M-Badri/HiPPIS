%
%---------------------------------------------------------------------------------------------%
% Driver to produce the  1D and 2D approximation presented in the manuscript
%---------------------------------------------------------------------------------------------%
%
  clear;
  close all;
  clc;
  fileID = fopen('approximations_tables_1d_2d.txt', 'w');

  % 1D function approximations %
  plot_approximations ;
  fprintf('Tables for 1D approximations saved in approximations_tables_1d_2d.txt \n');
 
  
  fprintf('Press any key to continue to the mapping examples. \n')
  fprintf('The approximated solutions are save in mapping_data/data. \n')
  fprintf('running plot_approximations.m and plot_mapping.m to  \n')
  fprintf('produce the figures and tables in the manuscript \n')
  
  plot_mapping ;
  
  fprintf('press any key to continue to the 2D function approximation examples.\n')
  fprintf(' WARNING: the 2D examples take long time because the solution is evaluated \n')
  fprintf('onto a 1000 x 1000 mesh for each example and saved. \n')

  % 2D function approximations %
  plot_approximations2D ;
  fprintf('Tables for 1D and 2D approximations saved in approximations_tables_1d_2d.txt \n');
  fclose(fileID);
  % end of script




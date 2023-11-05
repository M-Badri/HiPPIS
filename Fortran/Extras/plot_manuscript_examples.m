fprintf('Running the script in "plot_tables_figures.m" to produce plots and tables \n')
cd ../Drivers/Mapping/
addpath("../../Extras");
plot_tables_figures
%
cd ../BOMEX
fprintf('Running the script file "plot_bomex_paper.m" to produce BOMEX figures \n')
%
plot_bomex_paper
cd ../../Extras
%
fprintf('The tables in the BOMEX experiment are saved in "Drivers/BOMEX/bomex_tables.txt" ')
fprintf('The remaining tables are saved in tab "Drivers/Mapping/approximations_tables_1d_2d.txt" ')


fprintf('Moving into the "Mapping" folder and running the script in "/Mapping/plot_tables_figures.m" to produce plots and tables \n')
cd Mapping
plot_tables_figures

fprintf('Leaving the "Mapping" folder and moving into the "BOMEX" folder and running the script file "/BOMEX/plot_bomex_paper.m" to produce BOMEX figures \n')
cd ../BOMEX/
plot_bomex_paper
cd ../
fprintf('The tables in the BOMEX experiment are saved in "/BOMEX/bomex_tables.txt" ')
fprintf('The remaining tables are saved in tab "Mapping/approximations_tables_1d_2d.txt" ')


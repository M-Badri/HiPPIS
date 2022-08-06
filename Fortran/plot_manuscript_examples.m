fprintf('Moving into the "Mapping" folder and running script to produce plots and tables \n')
cd Mapping
plot_tables_figures
plot_mapping

fprintf('Leaving the "Mapping" folder and moving into the "BOMEX" folder and running script to produce BOMEX figures \n')
cd ../BOMEX/
plot_bomex_paper
cd ../


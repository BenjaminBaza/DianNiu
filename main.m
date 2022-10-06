%% Reset matlab
clc;
clear;
warning('off');

cd('\\fil031.uis.no\emp05\2925477\Documents\Battery research\DFN_DianNiu')

addpath("Parameters")
addpath("Coupled_solver")
addpath("Segregated_solver")
addpath("Coupled_solver")

%%%%%%%%%%%%%%%%% Set up global functions classes
global sol_fun
sol_fun = 		newton_solver_functions;
global vis_fun
vis_fun =		visual_functions;
global BC_fun
BC_fun =		boundary_condition_functions;
global eq_build_fun
eq_build_fun =	equation_building_functions;
global BV_fun
BV_fun =  BV_functions;

%%%%%%%%%%%%%%%%% Start timer
main_timer=tic;

%%%%%%%%%%%%%%%%% Read control parameters
read_timer=tic;

run read_ctrl_development
%run read_ctrl_template

read_chrono=toc(read_timer);
%%%%%%%%%%%%%%%%% Run solver
solver_timer=tic;

run DFN_solver;

solver_chrono=toc(solver_timer);
%%%%%%%%%%%%%%%%% Write data


%%%%%%%%%%%%%%%%% Plot results
plot_timer=tic;

if deb.plot_data==1
    disp('Ploting results')
    vis_fun.plot_data(sol.cell_center_coord,fv.pe,'x','Potential','Potential in electrolyte',10,'pe.pdf')
    vis_fun.plot_data(sol.cell_center_coord,fv.ce,'x','Concentration','Li concentration in electrolyte',10,'ce.pdf')
    vis_fun.plot_data_solid(sol.cell_center_coord,fv.ps,'x','Potential','Potential in solid',10,'ps.pdf',sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p)
    vis_fun.plot_data_solid(sol.cell_center_coord,fv.cse,'x','Concentration','Li concentration in solid',10,'cse.pdf',sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p)
    vis_fun.plot_data(sol.cell_center_coord,fv.j,'x','j','Rate of charge flow throught solid particle surface',10,'j.pdf')
    vis_fun.plot_data(sol.time_array,hist.V,'time (s)','Voltage (V)','Cell voltage over time',10,'V.pdf')
    vis_fun.plot_complete_data('complete_field_values.png',sol.cell_center_coord,2,sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p,sol.time_array);
    vis_fun.plot_solid_concentration_data('solid_c_values.png',sol.cell_center_coord,2,sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p,sol.time_array);
    vis_fun.plot_solid_concentration_singlePart_data('solid_c_singlePart_values.png',sol.cell_center_coord,2,sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p,sol.time_array);
    
    if deb.prints>=0
        vis_fun.plot_resuduals_DFN('DFN_res.png');
        vis_fun.plot_resuduals_diff('Diff_res.png');
        vis_fun.plot_resuduals_newt('Newt_res.png');
    end    
end

plot_chrono=toc(plot_timer);
%%%%%%%%%%%%%%%%% Generate videos
video_timer=tic;

if deb.animate_data==1
    disp('Animating results')
    vis_fun.animate_complete_data('complete_field_values.avi',sol.cell_center_coord,2,sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p,sol.time_array,hist.V);
end

video_chrono=toc(video_timer);
%%%%%%%%%%%%%%%%% Stop timer
chrono = toc(main_timer);
fprintf(1,'End of simulation \nThe full program takes %3.2f sec, \nread section %3.2f sec, \nsolver section %3.2f sec, \nplot section %3.2f sec, \nvideo section %3.2f sec,\n' ...
                                            ,chrono,read_chrono,solver_chrono,plot_chrono,video_chrono);

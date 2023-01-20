function main()
    
    %% Reset matlab
    clc;
    clear;
    warning('off');

    cd('\\fil031.uis.no\emp05\2925477\Documents\Battery research\DFN_DianNiu')

    addpath("Parameters")
    addpath("Coupled_solver")
    addpath("Segregated_solver")
    addpath("Coupled_solver")
    addpath("Mesh")
    addpath("External_code")


    %%%%%%%%%%%%%%%%% Set up global functions classes
    global sol_fun
    sol_fun =       newton_solver_functions;
    global vis_fun
    vis_fun =       visual_functions;
    global BC_fun
    BC_fun =        boundary_condition_functions;
    global eq_build_fun
    eq_build_fun =  equation_building_functions;
    global BV_fun
    BV_fun =        BV_functions;

    %disp(vis_fun.logarithmic_ticks_generator(0.000002,100000))

    %%%%%%%%%%%%%%%%% Start timer
    main_timer=tic;

    %%%%%%%%%%%%%%%%% Read control parameters
    read_timer=tic;

    run read_ctrl_development()
    %run read_ctrl_template

    global deb
    global hist
    global p
    global fv
    global sol

    read_chrono=toc(read_timer);
    %%%%%%%%%%%%%%%%% Run solver
    solver_timer=tic;

    DFN_solver();

    solver_chrono=toc(solver_timer);

    
    %%%%%%%%%%%%%%%%% Plot results
    plot_timer=tic;

    if deb.plot_data==1
        mkdir(deb.graph_folder_name) 

        disp('Ploting results')
        vis_fun.plot_data(sol.cell_center_coord,fv.pe,'x','Potential','Potential in electrolyte',10,'pe.pdf',"","")
        vis_fun.plot_data(sol.cell_center_coord,fv.ce,'x','Concentration','Li concentration in electrolyte',10,'ce.pdf',"","")
        vis_fun.plot_data_solid(sol.cell_center_coord,fv.ps,'x','Potential','Potential in solid',10,'ps.pdf',sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p)
        vis_fun.plot_data_solid(sol.cell_center_coord,fv.cse,'x','Concentration','Li concentration in solid',10,'cse.pdf',sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p)
        vis_fun.plot_data(sol.cell_center_coord,fv.j,'x','j','Rate of charge flow throught solid particle surface',10,'j.pdf',"","")
        % vis_fun.plot_data(sol.time_array(1:sol.time_ite),hist.V(1:sol.time_ite),'time (s)','Voltage (V)','Cell voltage over time',10,'V.pdf',"","")
        vis_fun.plot_data(sol.time_array(1:sol.time_ite),hist.V(1:sol.time_ite),'time (s)','Voltage (V)','Cell voltage over time',10,'V.pdf',"Chen_2020_model_1C.txt","Chen_2020_expe_1C.txt")
        % vis_fun.plot_complete_data('complete_field_values.png',sol.cell_center_coord,22289898,sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p,sol.time_array);
        vis_fun.plot_complete_data_spec_ite('complete_field_values_n_1.png',sol.cell_center_coord,24536789,sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p,sol.time_array,sol.time_ite_save-1);
        vis_fun.plot_complete_data_spec_ite('complete_field_values.png',sol.cell_center_coord,25679895,sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p,sol.time_array,sol.time_ite_save);
        vis_fun.plot_solid_concentration_data('solid_c_values.png',sol.cell_center_coord,2,sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p,sol.time_array);
        vis_fun.plot_solid_concentration_singlePart_data('solid_c_singlePart_values.png',sol.cell_center_coord,2,sol.nb_cell_n,sol.nb_cell_s,sol.nb_cell_p,sol.time_array);
        
        if deb.prints>=0
            vis_fun.plot_resuduals_DFN('DFN_res.png');
            vis_fun.plot_resuduals_diff('Diff_res.png');
            vis_fun.plot_resuduals_newt('Newt_res.png',sol.newt_ite);
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

    %%%%%%%%%%%%%%%%% Write data
    write_timer=tic;

    if deb.write_output_data>=1
        save_data();
    end

    write_chrono=toc(write_timer);


    %%%%%%%%%%%%%%%%% Stop timer
    chrono = toc(main_timer);
    fprintf(1,'End of simulation \n\nThe full program takes %3.2f sec \nread section %3.2f sec \nsolver section %3.2f sec \nplot section %3.2f sec \nvideo section %3.2f sec \nwrite section %3.2f sec\n\n' ...
                                                , chrono, read_chrono, solver_chrono, plot_chrono, video_chrono, write_chrono);
    fprintf(1,'Total time of solver setup %3.2f sec \nTotal time of matrix inversion %3.2f sec \nTotal time of solver update %3.2f sec \nTotal time in solver iteration %3.2f sec \n' ...
                                                , deb.chrono_newtsol_setup, deb.chrono_matrix_inversion, deb.chrono_newtsol_update...
                                                , deb.chrono_newtsol_setup+deb.chrono_matrix_inversion+deb.chrono_newtsol_update);
    deb.chrono_matrix_inversion=0;
    deb.chrono_newtsol_setup=0;
    deb.chrono_newtsol_update=0;

    clear video_chrono video_timer main_timer read_chrono read_timer solver_chrono solver_timer write_chrono write_timer chrono
    clear plot_chrono plot_timer pe_next_save ps_next_save M_electrlyte_diff Mn_part_diff Mp_part_diff csp_next csn_next ce_next_save sol.time_ite sol.time_ite_save
    clear rms folder_name

    beep

end
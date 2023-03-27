function read_ctrl_development()

    %% This list of parameters is meant to simulate a cylindrical 21700 commercial cell (LGM50).
    %Ref from Chen et al. 2020 "Development of experimental techniques for parameterization"
    % and Dandelion library (https://simulation.dandeliion.com/simulation/bf981e79-ffc8-4f57-b884-42af8a4554d9/)

    %reminders : in the case of a discharge cathode=positive electrode and anode=negative.

    global BV_fun

    %% Parameter functions class
    global param_functions
    param_functions=param_functions_LGM50;

    %% Geometric characteristics

    % Lengths
    global p
    p.Ln = 85.2e-6;     % negative electrode thickness [m]
    p.Lp = 75.6e-6;     % positive electrode thickness [m]
    p.Ls = 12e-6;      % separator thickness [m]
    p.Ltot = p.Ln+p.Ls+p.Lp; % total cell thickness
    p.R_s_n = 5.86e-6;   % negative electrode particle radius [m]
    p.R_s_p = 5.22e-6;   % positive electrode particle radius [m]

    % Volume fractions
    p.eps_e_n = 0.25;   % negative electrode volume fraction of electrolyte 
    p.eps_e_s = 0.47;   % separator volume fraction of electrolyte 
    p.eps_e_p = 0.335;   % negative electrode volume fraction of electrolyte 
    p.eps_s_n = 0.75;      % negative electrode volume fraction (in the absence of filler, this is 1-p.eps_e_n)
    p.eps_s_p = 0.665;      % positive electrode volume fraction (in the absence of filler, this is 1-p.eps_e_p)


    % Specific interfacial surface area
    p.coll_A_n = 0.1027; % negative electrode current collector area [m^2]
    p.coll_A_p = 0.1027; % positive electrode current collector area [m^2]
    p.R_collector_contact=0;	% Current collector total contact resistance [V.A^(-1)] [Ohm]
    p.A_s_n = 3*p.eps_s_n / p.R_s_n;  % Negative electrode specific interfacial surface area [m^2/m^3]
    p.A_s_p = 3*p.eps_s_p / p.R_s_p;  % Positive electrode specific interfacial surface area [m^2/m^3]


    %% Solver settings  
    global sol

    sol.newt_ite = 0;
    sol.time_ite = 0;


    sol.time_tot    = 2000.0;                    %Total time of the simulation [s]
    sol.time        = 0.0; %10800;                     %Total time of the simulation [s]
    sol.dt          = 0.    ;%1.                        %Time step for the time discretization [s]
    sol.max_dt      = 10.    ;%1.                        %Maximum time step for the time discretization [s]
    sol.quick_dt    = 0.5    ;%1.                        %Faster time step for the time discretization [s] for the modeling of the end of a charge / discharge
    sol.nb_steps    = 10000 ; %length(sol.time_array);    % Number of time states visited throughout the simulation
    sol.time_array  = zeros(1,sol.nb_steps) ; %sol.dt:sol.dt:sol.time_tot;  % array containing the time coordinate of each time step

    sol.max_allowed_voltage = 100.0 ; %[V] (Not used in the solver at the moment)
    sol.min_allowed_voltage = 2.5 ; %[V] (Not used in the solver at the moment)
    sol.max_allowed_voltage_time_differential = 2.;

    sol.nb_cell_n   = 40;
    sol.nb_cell_s   = 3;
    sol.nb_cell_p   = 40;
    sol.nb_cell     = sol.nb_cell_n + sol.nb_cell_s + sol.nb_cell_p ;   %
    sol.nb_cell_ps     = sol.nb_cell_n + sol.nb_cell_p ;   %

    sol.dxn         = p.Ln/sol.nb_cell_n;
    sol.dxs         = p.Ls/sol.nb_cell_s;
    sol.dxp         = p.Lp/sol.nb_cell_p;

    sol.vertex_coord_n  = 0:sol.dxn:p.Ln;
    sol.vertex_coord_s  = p.Ln+sol.dxs:sol.dxs:p.Ln+p.Ls;
    sol.vertex_coord_p  = p.Ln+p.Ls+sol.dxp:sol.dxp:p.Ln+p.Ls+p.Lp;


    sol.vertex_coord  = cat(2,sol.vertex_coord_n,sol.vertex_coord_s);
    sol.vertex_coord  = cat(2,sol.vertex_coord,sol.vertex_coord_p);
    sol.cell_center_coord  =ones(sol.nb_cell,1) ;
    for iiii = 1:1:sol.nb_cell
    	sol.cell_center_coord(iiii)=(sol.vertex_coord(iiii+1)+sol.vertex_coord(iiii))/2. ;
    end
    clear iiii

    sol.cell_dx  = cat(1,sol.dxn*ones(sol.nb_cell_n,1),sol.dxs*ones(sol.nb_cell_s,1));
    sol.cell_dx  = cat(1,sol.cell_dx,sol.dxp*ones(sol.nb_cell_p,1));
    sol.cell_center_coord_solid  = cat(1,sol.cell_center_coord(1:sol.nb_cell_n),sol.cell_center_coord(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell_n+sol.nb_cell_s+sol.nb_cell_p));


    sol.part_nb_cell= 39;%20;
    sol.particle_mesh_distribution_order= 2.0;
    particle_mesh_generator(sol.particle_mesh_distribution_order)

    sol.newton_meth_res_threshold=1e-6	;
    sol.newton_meth_max_ite=39;
    sol.newton_relax_factor = 0.8;


    %% Debug parameters

    global deb

    deb.prints=0;
    deb.run_name= "LGM50";
    deb.plot_data=1;
    deb.animate_data=0;
    deb.break_time_loop=0;
    deb.write_output_data=2; %1 to write final data, 2 to write all history data

    deb.timing_jacobian=0;
    deb.chrono_matrix_inversion=0;
    deb.chrono_newtsol_setup=0;
    deb.chrono_newtsol_update=0;
    deb.chrono_matrix_inversion_singleite=0;
    deb.chrono_newtsol_setup_singleite=0;
    deb.chrono_newtsol_update_singleite=0;

    deb.chrono_updt_dtso_singleite=0;
    deb.chrono_Calc_Jac_singleite=0;
    deb.chrono_Calc_f_singleite=0;
    deb.chrono_Calc_gJacg_singleite=0;
    deb.chrono_assemble_coupled_singleite=0;

    deb.chrono_Jacdiag_singleite = 0;
    deb.chrono_Jacce_singleite = 0;
    deb.chrono_Jacpe_singleite = 0;
    deb.chrono_Jacps_singleite = 0;
    deb.chrono_Jaccsn_singleite = 0;
    deb.chrono_Jaccsp_singleite = 0;

    deb.chrono_Jaccedce_singleite = 0;
    deb.chrono_Jacpedpe_singleite = 0;
    deb.chrono_Jacpsdps_singleite = 0;
    deb.chrono_Jaccsdcs_singleite = 0;

    deb.case_name=datestr(now,'mm_dd_yy_')+deb.run_name+'_'+num2str(sol.nb_cell_n)+'x'+num2str(sol.nb_cell_s)+'x'+num2str(sol.nb_cell_p)+'p'+num2str(sol.part_nb_cell)+'_save_';

    save_nb=1;
    while exist('../Saved_data_DFN_DianNiu/'+deb.case_name+num2str(save_nb), 'dir')
        save_nb=save_nb+1;
    end
    deb.case_name=deb.case_name+num2str(save_nb);
    deb.folder_name=append('../Saved_data_DFN_DianNiu/',deb.case_name);
    deb.read_ctrl_file_name = 'Parameters/read_ctrl_development.m' ;
    deb.graph_folder_name=deb.folder_name+'/Graphs/';

    %% Material characteristics

    % Constants
    p.k0 = 1e-10 * [0.0672,0.3545] ;	% Effective reaction rates (in the negative electrode first and in the positive second) [mol^−1/2 .m^5/2 .s−1]
    p.alpha	=0.5 ;			% Asymmetric charge-transfer coefficient
    p.brug = 1.5 ;       	% Bruggeman porosity
    p.t_plus_function_mode=0 ; %The transference number is handled as a function 1 or as a constant 0
    p.t_plus = 0.2590 ;       	% Transference number
    p.Faraday = 96487 ;    	% Faraday constant, [Coulumbs/mol]  [s.A/mol]  [kg.m^2.s^-2.V^-1/mol]
    p.Rg = 8.3145 ;			% Perfect gas constant [j.mol/K]  [kg.m^2.s^-2/mol/K]
    p.Rfilm	= 0 ;			% Film layer ionic resistance [Ohm.m^2]


    % Maximum concentrations
    p.csn_max = 33133.0;    % Max concentration in anode, [mol/m^3]
    p.csp_max = 63104.0;    % Max concentration in cathode, [mol/m^3]
    p.neg_stoichiometry_min = 0.0279;	% Minimum concentration at   0% stochiometry in anode, [mol/m^3]
    p.neg_stoichiometry_max = 0.9014;	% Maximum concentration at 100% stochiometry in anode, [mol/m^3]
    p.pos_stoichiometry_min = 0.2661;	% Minimum concentration at 100% stochiometry in cathode, [mol/m^3]
    p.pos_stoichiometry_max = 0.9084;	% Maximum concentration at   0% stochiometry in cathode, [mol/m^3]

    p.csn_stoic_min = p.csn_max*p.neg_stoichiometry_min;   % concentration at 100% stochiometry in anode, [mol/m^3]
    p.csp_stoic_min = p.csp_max*p.pos_stoichiometry_min;   % concentration at   0% stochiometry in cathode, [mol/m^3]
    p.csn_stoic_max = p.csn_max*p.neg_stoichiometry_max;   % concentration at   0% stochiometry in anode, [mol/m^3]
    p.csp_stoic_max = p.csp_max*p.pos_stoichiometry_max;   % concentration at 100% stochiometry in cathode, [mol/m^3]


    % Solid phase diffusion coefficients
    p.Dsn_array = zeros(sol.nb_cell_n*(sol.part_nb_cell+1),1);  % neg. electrode diffusion coeff [m^2/s]
    p.Dsp_array = zeros(sol.nb_cell_n*(sol.part_nb_cell+1),1);  % pos. electrode diffusion coeff [m^2/s]

    p.Ds_constant=1;

    % Electrolyte diffusion coefficient
    p.De_eff  = cat(1,ones(sol.nb_cell_n,1),ones(sol.nb_cell_s,1));
    p.De_eff  = cat(1,p.De_eff,ones(sol.nb_cell_p,1));

    p.De_eff_coeff = cat(1,p.eps_e_n^p.brug*ones(sol.nb_cell_n,1),p.eps_e_s^p.brug*ones(sol.nb_cell_s,1));
    p.De_eff_coeff  = cat(1,p.De_eff_coeff,p.eps_e_p^p.brug*ones(sol.nb_cell_p,1));

    p.De_constant=0;

    % Solid phase conductivity
    p.sig_n = 215.;    % Conductivity of solid in neg. electrode, [1/Ohm*m]=[S*m] (S=Siemens)
    p.sig_p = 0.18;    % Conductivity of solid in pos. electrode, [1/Ohm*m]=[S*m] (S=Siemens)

    p.sig_eff_n = p.sig_n * p.eps_s_n^p.brug;    % Effective conductivity in neg. electrode, [A/m/V]
    p.sig_eff_p = p.sig_p * p.eps_s_p^p.brug;    % Effective conductivity in pos. electrode, [A/m/V]
    p.sig_eff  = p.sig_eff_n*ones(sol.nb_cell_n,1);
    p.sig_eff  = cat(1,p.sig_eff,p.sig_eff_p*ones(sol.nb_cell_p,1));


    %% External input
    global ex
    ex.I_array  = 10. * ones(size(sol.time_array)) ; %0.5 * ones(size(sol.time_array));%ex.temporary_Crate*ones(size(sol.time_array));  % Array containing the input current intensity at each time step (may be redundant)
    % 1C=5A


    %% Initial conditions
    global ini
    ini.V0 = 4.0;       % Initial tension [V]
    ini.ce0 = 1000. ;%1000 ;      % Initial electrolyte concentration of Li, [mol/m^3]
    ini.SOC = 1.;
    ini.csp0 = p.csp_stoic_min+ (p.csp_stoic_max-p.csp_stoic_min)*(1-ini.SOC) ;  % Initial pos. electrode concentration of Li, [mol/m^3]
    ini.csn0 = p.csn_stoic_min+ (p.csn_stoic_max-p.csn_stoic_min)*(  ini.SOC) ;  % Initial neg. electrode concentration of Li, [mol/m^3]
    ini.T0 = 298.15;%1e3;       % Initial temperature, [K]


    % Conductivity of electrolyte
    p.kappa_eff  = cat(1,ones(sol.nb_cell_n,1),ones(sol.nb_cell_s,1));
    p.kappa_eff  = cat(1,p.kappa_eff,ones(sol.nb_cell_p,1));
    p.kappa_D_eff= p.kappa_eff;

    p.kappa_constant=0;

    %% Field variables allocation

    global fv
    fv.csn = ini.csn0 * ones(length(sol.part_coord_n),sol.nb_cell_n);   % Lithium concentration in the neg. electrode
    fv.csp = ini.csp0 * ones(length(sol.part_coord_p),sol.nb_cell_p);   % Lithium concentration in the pos. electrode
    fv.cse = cat(2,fv.csn(length(sol.part_coord_n),:),fv.csp(length(sol.part_coord_p),:));
     

    fv.ce  = ini.ce0 * ones(sol.nb_cell,1);
    fv.Ueq = 0.*ones(sol.nb_cell,1); % Open-circuit potential [V]


  	fv.Ueq(1:sol.nb_cell_n) = param_functions.neg_electrode_Ueq(ini.csn0,0)* ones(sol.nb_cell_n,1);
   	fv.Ueq(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell) = param_functions.pos_electrode_Ueq(ini.csp0,0)* ones(sol.nb_cell_p,1);

    ini.pe0 = 0.;      % Initial electrolyte electric potential, [V]
    ini.psn0 = param_functions.neg_electrode_Ueq(ini.csn0,0);      % Initial positive electrode electric potential equal to the equilibrium potential, [V]
    ini.psp0 = param_functions.pos_electrode_Ueq(ini.csp0,0);      % Initial negative electrode electric potential equal to the equilibrium potential, [V]

    if deb.prints>5
    	figure(268989898);
        fs = 16;
        set(gcf,'Position',[50 50 1800 1000]);     

        subplot(2,5,[1,2])
        plot(0:p.csn_max/50.:p.csn_max,param_functions.neg_electrode_Ueq(0:p.csn_max/50.:p.csn_max),'LineWidth',2);
        title('neg Ueq','fontsize',fs);
        grid on
        grid minor


        subplot(2,5,[4,5])
        plot(0:p.csp_max/50.:p.csp_max,param_functions.pos_electrode_Ueq(0:p.csp_max/50.:p.csp_max),'LineWidth',2);
        title('pos Ueq','fontsize',fs);
        grid on
        grid minor

        c_array=0:p.csp_max/50.:p.csp_max;
        for c_ind=1:length(c_array)
        	cond_array(c_ind)=param_functions.electrolyte_conductivity(c_array(c_ind));
        	diff_array(c_ind)=param_functions.electrolyte_diffusivity(c_array(c_ind));
        end

        subplot(2,5,[6,7])
        plot(c_array,cond_array,'LineWidth',2);
        title('electrolyteconductivity','fontsize',fs);
        grid on
        grid minor


        subplot(2,5,[9,10])
        plot(c_array,diff_array,'LineWidth',2);
        title('electrolyte diffusivity','fontsize',fs);
        grid on
        grid minor

        saveas(gcf,deb.graph_folder_name+'Parameter_functions');
    end

    fv.pe  = ini.pe0 * ones(sol.nb_cell,1);
    fv.ps  = cat(1,ini.psn0 * ones( sol.nb_cell_n , 1),ini.psp0 * ones( sol.nb_cell_p , 1));

    fv.j=BV_fun.butler_volmer_equation(fv.pe,fv.ps,fv.ce,fv.cse, ini.T0,fv.Ueq,"read_ctrl");
    fv.V=0;
    fv.SOC_neg=ini.SOC;
    fv.SOC_pos=ini.SOC;

   	p.De_eff = param_functions.electrolyte_diffusivity(ini.ce0) *ones(sol.nb_cell,1);

   	p.kappa_eff = param_functions.electrolyte_conductivity(ini.ce0) *ones(sol.nb_cell,1);
   	p.kappa_D_eff = p.kappa_eff *2* p.Rg * ini.T0/p.Faraday * (p.t_plus-1);

    p.Dsn_array = param_functions.neg_electrode_diffusivity(ini.csn0) *ones(sol.nb_cell_n*(sol.part_nb_cell+1),1);
    p.Dsp_array = param_functions.pos_electrode_diffusivity(ini.csp0) *ones(sol.nb_cell_p*(sol.part_nb_cell+1),1);

    %% Field variables history tables allocation

    global hist

    hist.csn = zeros(length(reshape(fv.csn,sol.nb_cell_n*(sol.part_nb_cell+1),1)),sol.nb_steps+1);
    hist.csp = zeros(length(reshape(fv.csp,sol.nb_cell_p*(sol.part_nb_cell+1),1)),sol.nb_steps+1);
    hist.cse = zeros(length(fv.cse),sol.nb_steps+1);
    hist.ce  = zeros(length(fv.ce),sol.nb_steps+1);
    hist.pe  = zeros(length(fv.pe),sol.nb_steps+1);
    hist.ps  = zeros(length(fv.ps),sol.nb_steps+1);
    hist.sourceps  = zeros(length(fv.ps),sol.nb_steps+1);
    hist.sourcepsBC  = zeros(sol.nb_steps+1,1);
    hist.Ueq = zeros(length(fv.Ueq),sol.nb_steps+1);
    hist.j   = zeros(length(fv.j),sol.nb_steps+1);
    hist.V   = zeros(1,sol.nb_steps);
    hist.V(1)   = ini.psp0 - ini.psn0;
    hist.SOC_neg   = zeros(1,sol.nb_steps);
    hist.SOC_pos   = zeros(1,sol.nb_steps);
    hist.charge_time=0;
    hist.charge_ite=2;


    hist.csn(:,1) = reshape(fv.csn,sol.nb_cell_n*(sol.part_nb_cell+1),1);
    hist.csp(:,1) = reshape(fv.csp,sol.nb_cell_p*(sol.part_nb_cell+1),1);
    hist.cse(:,1) = fv.cse;
    hist.ce(:,1)  = fv.ce;
    hist.pe(:,1)  = fv.pe;
    hist.ps(:,1)  = fv.ps;
    hist.Ueq(:,1) = fv.Ueq;
    hist.j(:,1)   = fv.j;
    hist.residuals   		= zeros(6,sol.newton_meth_max_ite);
    hist.residuals_time		= zeros(6,sol.nb_steps);
    hist.residuals_diff		= zeros(6,sol.nb_steps);
    hist.newt_it_number		= zeros(1,sol.nb_steps);

    hist.sum_j_pos_electrode    = zeros(1,sol.nb_steps);
    hist.sum_j_neg_electrode    = zeros(1,sol.nb_steps);

    hist.maxjn= -1000000000;
    hist.minjn=  1000000000;
    hist.maxjp= -1000000000;
    hist.minjp=  1000000000;

    hist.maxdt= -1000000000;
    hist.mindt=  1000000000;
    
    hist.maxcsn=-1000000000;
    hist.mincsn= 1000000000;
    hist.maxcsp=-1000000000;
    hist.mincsp= 1000000000;

    hist.maxdcsn=0;
    hist.mindcsn=0;
    hist.maxdcsp=0;
    hist.mindcsp=0;
    
    hist.maxdcs_0csn=0;
    hist.maxdcs_Rcsn=0;
    hist.mindcs_0csn=0;
    hist.mindcs_Rcsn=0;

    hist.maxdcs_0csp=0;
    hist.maxdcs_Rcsp=0;
    hist.mindcs_0csp=0;
    hist.mindcs_Rcsp=0;

end
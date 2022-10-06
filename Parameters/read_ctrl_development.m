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
p.L_ccn = 11.7e-6;    % negative current collector thickness [m]
p.L_ccp = 16.3e-6;    % positive current collector thickness [m]

% Volume fractions
p.eps_e_n = 0.25;   % negative electrode volume fraction of electrolyte 
p.eps_e_s = 0.47;   % separator volume fraction of electrolyte 
p.eps_e_p = 0.335;   % negative electrode volume fraction of electrolyte 
p.eps_s_n = 0.75;      % negative electrode volume fraction
p.eps_s_p = 0.665;      % positive electrode volume fraction
p.eps_f_n = 0.;  % negative electrode volume fraction of filler
p.eps_f_p = 0.;  % negative electrode volume fraction of filler


% Specific interfacial surface area
p.coll_A_n = 0.1027; % negative electrode current collector area [m^2]
p.coll_A_p = 0.1027; % positive electrode current collector area [m^2]
p.R_collector_contact=0;
p.A_s_n = 3*p.eps_s_n / p.R_s_n;  % Negative electrode specific interfacial surface area [m^2/m^3]
p.A_s_p = 3*p.eps_s_p / p.R_s_p;  % Positive electrode specific interfacial surface area [m^2/m^3]


%% Solver settings
global sol

sol.coupling_scheme=1;

sol.time_tot    = 50.;%10800;                     %Total time of the simulation [s]
sol.dt          = 10.;%1.                        %Time step for the time discretization [s]
sol.time_array  = sol.dt:sol.dt:sol.time_tot;  % array containing the time coordinate of each time step (may be redundant)
sol.nb_steps    =length(sol.time_array);    % Number of time states visited throughout the simulation

sol.max_allowed_voltage = 4.40 ; %[V] (Not used in the solver at the moment)
sol.min_allowed_voltage = 2.0 ; %[V] (Not used in the solver at the moment)

sol.nb_cell_n   = 5;%30;%50;
sol.nb_cell_s   = 3;%20;%50;
sol.nb_cell_p   = 5;%30;%50;
sol.nb_cell     = sol.nb_cell_n + sol.nb_cell_s + sol.nb_cell_p ;   %

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


sol.part_nb_cell= 10;%20;
sol.dxsn        = p.R_s_n/(sol.part_nb_cell);
sol.dxsp        = p.R_s_p/(sol.part_nb_cell);
sol.part_coord_n  = 0:sol.dxsn:p.R_s_n;
sol.part_coord_p  = 0:sol.dxsp:p.R_s_p;
sol.newton_meth_res_threshold=1e-2;
sol.newton_meth_max_ite=500;
sol.newton_update_sources= 1;
sol.newton_relax_factor = 0.8;

%% Debug parameters

global deb

deb.prints=3;
deb.videos_generation=0;
deb.plot_data=0;
deb.animate_data=0;
deb.safe_BC_mode=0;
deb.break_time_loop=0;

%% Material characteristics

% Mass densities [kg/m^3]
p.rho_sn = 0;  %NOT USED  % negative electrode solid phase density  
p.rho_sp = 0;  %NOT USED  % positive electrode solid phase density 
p.rho_e =  0;  %NOT USED  % Electrolyte density
p.rho_f = 0;   %NOT USED  % Filler density
p.rho_ccn = 0; %NOT USED  % negative electrode current collector density
p.rho_ccp = 0; %NOT USED  % positive electrode current collector density

% Masses [kg/m^2]
p.Mn = 		p.Ln * (p.rho_e*p.eps_e_n + p.rho_sn*p.eps_s_n + p.rho_f*p.eps_f_n);
p.Ms = 		p.Ls * (p.rho_e*p.eps_e_n);
p.Mp = 		p.Lp * (p.rho_e*p.eps_e_p + p.rho_sp*p.eps_s_p + p.rho_f*p.eps_f_p);
p.Mcc = 	p.rho_ccn*p.L_ccn + p.rho_ccp*p.L_ccp;
p.Mtot = 	p.Mn + p.Ms + p.Mp + p.Mcc;


% Constants
p.k0 = 1e-10 * [0.0672,0.3545] ;			% Effective reaction rates (in the negative electrode first and in the positive second) [mol^−1/2 .m^5/2 .s−1]
p.alpha	=0.5 ;			% Asymmetric charge-transfer coefficient
p.brug = 1.5 ;       	% Bruggeman porosity
p.t_plus_function_mode=0; %The transference number is handled as a function 1 or as a constant 0
p.t_plus = 0.2590 ;       	% Transference number
p.Faraday = 96487 ;    	% Faraday constant, [Coulumbs/mol]  [s.A/mol]  [kg.m^2.s^-2.V^-1/mol]
p.Rg = 8.3145 ;			% Perfect gas constant [j.mol/K]  [kg.m^2.s^-2/mol/K]
p.Rfilm	= 0 ;			% Film layer ionic resistance [Ohm.m^2]


% Maximum concentrations
temp_c_factor=1.;
p.csn_max = 33133.0* temp_c_factor; % 3.6e3 * 372 * 1800 / p.Faraday;   % Max concentration in anode, [mol/m^3]
p.csp_max = 63104.0* temp_c_factor; % 3.6e3 * 247 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]


% Solid phase diffusion coefficients
p.Dsn = 3.3e-4; %3.9e-14;  % neg. electrode diffusion coeff [m^2/s]
p.Dsp = 4.0e-5; %1e-13;  % pos. electrode diffusion coeff [m^2/s]


% Electrolyte diffusion coefficient
p.activation_energy_electrolyte =17100.0; % Activation energy for conductivity/diffusivity in electrolyte, [J·mol-1]
p.De_function_mode=1; %The diffusivity of the electrolyte is handled as a function 1 or as a constant 0
p.De = 1 ;%2.7877e-10;    % Diffusion coeff for electrolyte, [m^2/s]
p.De_eff_n = p.De * p.eps_e_n^p.brug;    % Effective diffusion coeff for the electrolyte in the neg. electrode, [m^2/s]
p.De_eff_s = p.De * p.eps_e_s^p.brug;    % Effective diffusion coeff for the electrolyte in the separator electrode, [m^2/s]
p.De_eff_p = p.De * p.eps_e_p^p.brug;    % Effective diffusion coeff for the electrolyte in the pos. electrode, [m^2/s]
p.De_eff  = cat(1,p.De_eff_n*ones(sol.nb_cell_n,1),p.De_eff_s*ones(sol.nb_cell_s,1));
p.De_eff  = cat(1,p.De_eff,p.De_eff_p*ones(sol.nb_cell_p,1));

p.De_eff_coeff = cat(1,p.eps_e_n^p.brug*ones(sol.nb_cell_n,1),p.eps_e_s^p.brug*ones(sol.nb_cell_s,1));
p.De_eff_coeff  = cat(1,p.De_eff_coeff,p.eps_e_p^p.brug*ones(sol.nb_cell_p,1));

% Solid phase conductivity
p.sig_n = 215.;    % Conductivity of solid in neg. electrode, [1/Ohm*m]=[S*m] (S=Siemens)
p.sig_p = 0.18;    % Conductivity of solid in pos. electrode, [1/Ohm*m]=[S*m] (S=Siemens)

p.sig_eff_n = p.sig_n * p.eps_s_n^p.brug;    % Effective conductivity in neg. electrode, [A/m/V]
p.sig_eff_p = p.sig_p * p.eps_s_p^p.brug;    % Effective conductivity in pos. electrode, [A/m/V]
p.sig_eff  = p.sig_eff_n*ones(sol.nb_cell_n,1);
p.sig_eff  = cat(1,p.sig_eff,p.sig_eff_p*ones(sol.nb_cell_p,1));


%% External input
ex.temporary_Crate = 100000   ;                      % 
ex.I_array  = 0.5 * ones(size(sol.time_array));%ex.temporary_Crate*ones(size(sol.time_array));  % Array containing the input current intensity at each time step (may be redundant)


%% Initial conditions
global ini
ini.V0 = 4.0;       % Initial tension [V]
%ini.ce0 = 1e3;      % Initial electrolyte concentration of Li, [mol/m^3]
%ini.csp0 = 16792.0;% 3.45e4;  % Initial pos. electrode concentration of Li, [mol/m^3]
%ini.csn0 = 29866.1;%1.29e4;  % Initial neg. electrode concentration of Li, [mol/m^3]
ini.ce0 = 1000 * temp_c_factor;      % Initial electrolyte concentration of Li, [mol/m^3]
ini.csp0 =  3.45e4 * temp_c_factor;  % Initial pos. electrode concentration of Li, [mol/m^3]
ini.csn0 =  1.29e4 * temp_c_factor;  % Initial neg. electrode concentration of Li, [mol/m^3]
ini.T0 = 298.15;%1e3;       % Initial temperature, [K]



% Conductivity of electrolyte
p.kappa_function_mode=1; %The conductivity of the electrolyte is handled as a function 1 or as a constant 0
p.kappa = 1; % Ionic conductivity of electrolyte [A/m/V]
p.kappa_eff_n = p.kappa * p.eps_e_n ^p.brug; % Ionic conductivity of electrolyte in neg. electrode [A/m/V]
p.kappa_eff_s = p.kappa * p.eps_e_s ^p.brug; % Ionic conductivity of electrolyte in separator electrode [A/m/V]
p.kappa_eff_p = p.kappa * p.eps_e_p ^p.brug; % Ionic conductivity of electrolyte in pos. electrode [A/m/V]
p.kappa_eff  = cat(1,p.kappa_eff_n*ones(sol.nb_cell_n,1),p.kappa_eff_s*ones(sol.nb_cell_s,1));
p.kappa_eff  = cat(1,p.kappa_eff,p.kappa_eff_p*ones(sol.nb_cell_p,1));
p.kappa_D_eff= p.kappa_eff*2* p.Rg * ini.T0/p.Faraday * (p.t_plus-1);

%% Field variables allocation

global fv
fv.csn = ini.csn0 * ones(length(sol.part_coord_n),sol.nb_cell_n);   % Lithium concentration in the neg. electrode
fv.csp = ini.csp0 * ones(length(sol.part_coord_p),sol.nb_cell_p);   % Lithium concentration in the pos. electrode
fv.cse = cat(2,fv.csn(length(sol.part_coord_n),:),fv.csp(length(sol.part_coord_p),:));
 

fv.ce  = ini.ce0 * ones(sol.nb_cell,1);
fv.Ueq = 0.*ones(sol.nb_cell,1); % Open-circuit potential [V]


if p.kappa_function_mode==1
	fv.Ueq(1:sol.nb_cell_n) = param_functions.neg_electrode_Ueq(ini.csn0,0)* ones(sol.nb_cell_n,1);
	fv.Ueq(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell) = param_functions.pos_electrode_Ueq(ini.csp0,0)* ones(sol.nb_cell_p,1);
	if deb.prints>2
		disp(fv.Ueq)
	end
end

ini.pe0 = 0.;      % Initial electrolyte electric potential, [V]
ini.psn0 = param_functions.neg_electrode_Ueq(ini.csn0,0);      % Initial positive electrode electric potential, [V]
ini.psp0 = param_functions.pos_electrode_Ueq(ini.csp0,0);      % Initial negative electrode electric potential, [V]

fv.pe  = ini.pe0 * ones(sol.nb_cell,1);
fv.ps  = cat(1,ini.psn0 * ones( sol.nb_cell_n , 1),ini.psp0 * ones( sol.nb_cell_p , 1));

fv.j=BV_fun.butler_volmer_equation(fv.pe,fv.ps,fv.ce,fv.cse,p.k0,p.alpha,p.Faraday,p.Rg, ini.T0,fv.Ueq,p.Rfilm,p.csn_max,p.csp_max,sol.nb_cell_n,sol.nb_cell_s);
fv.V=0;

if p.De_function_mode==1
	p.De_eff = param_functions.electrolyte_diffusivity(ini.ce0) *ones(sol.nb_cell,1);
end

if p.kappa_function_mode==1
	p.kappa_eff = param_functions.electrolyte_conductivity(ini.ce0) *ones(sol.nb_cell,1);
	p.kappa_D_eff = p.kappa_eff *2* p.Rg * ini.T0/p.Faraday * (p.t_plus-1);
end



%% Field variables history tables allocation

global hist

hist.csn = zeros(length(reshape(fv.csn,sol.nb_cell_n*(sol.part_nb_cell+1),1)),sol.nb_steps+1);
hist.csp = zeros(length(reshape(fv.csp,sol.nb_cell_p*(sol.part_nb_cell+1),1)),sol.nb_steps+1);
hist.cse = zeros(length(fv.cse),sol.nb_steps+1);
hist.ce  = zeros(length(fv.ce),sol.nb_steps+1);
hist.pe  = zeros(length(fv.pe),sol.nb_steps+1);
hist.ps  = zeros(length(fv.ps),sol.nb_steps+1);
hist.Ueq = zeros(length(fv.Ueq),sol.nb_steps+1);
hist.j   = zeros(length(fv.j),sol.nb_steps+1);
hist.V   = zeros(1,sol.nb_steps);

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




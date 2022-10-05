disp('Solver starts')

ite=0;
%%start time stepping
for time=sol.time_array
    ite=ite+1;

    disp("DFN_solver iteration "+int2str(ite)+" at time "+num2str(time)+" sec")

    if sol.coupling_scheme==0
        %% Solve concentration in the solid particles.
        run concentration_solid;
        
        %% Solve concentration in the electrolyte
        run concentration_electrolyte;

    elseif sol.coupling_scheme==1
        %% Solve the eqautions of the DFN model (potential in the solid and electrolyte (ps and pe) and concentration in the electrolyte and in the solid (ce, cs)).
        run DFN_equations_coupled_solver;
    end

    if deb.break_time_loop==1
        disp("The Newton solver has diverged.Time solver iterrupted.")
        break
    end

    %% Calculate voltage at time t
    fv.V=voltage_calc(fv.ps,p.R_collector_contact ,ex.I_array(ite));

    %% Add field variables to history records for later animation and output data
    save_fv_to_hist(ite)
    
end

clear time ite;

disp('Solver ends');
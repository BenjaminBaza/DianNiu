disp('Solver starts')

time_ite=0;

%% Start Newton solver ite timer
newt_timer=tic;
newt_chrono=0;

%%start time stepping
for time=sol.time_array
    time_ite=time_ite+1;


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
    [fv.V,fv.SOC_neg,fv.SOC_pos]=voltage_calc(fv.ps,p.R_collector_contact ,ex.I_array(time_ite));

    newt_ite_chrono= toc(newt_timer)-newt_chrono;
    disp("DFN_solver ite "+int2str(time_ite)+" at time "+num2str(time)+"s converged in "+num2str(newt_ite)+"ite "+num2str(newt_ite_chrono)+"s  voltage="+num2str(fv.V)+"V  SOCn="+num2str(fv.SOC_neg)+"  SOCp="+num2str(fv.SOC_pos))
    newt_chrono=toc(newt_timer);


    %% Add field variables to history records for later animation and output data
    save_fv_to_hist(time_ite)
    
end

clear time ;

disp('Solver ends');
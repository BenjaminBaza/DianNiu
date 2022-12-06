
function DFN_solver()

    global sol_fun
    global BV_fun
    global ini
    global deb
    global sol
    global ex
    global p
    global fv
    global hist
    global eq_build_fun

    user = memory;

    equations_nb=sol.nb_cell*3-sol.nb_cell_s + (sol.nb_cell_n+sol.nb_cell_p)*(sol.part_nb_cell+1);
    disp("Solver starts. The problem contains "+num2str(equations_nb)+" equations")


    sol.time_ite=0;

    %% Start Newton solver ite timer
    newt_timer=tic;
    newt_chrono=0;

    %%start time stepping
    while sol.time<sol.time_tot
        sol.time_ite=sol.time_ite+1;

        %sol.dt=sol.max_dt;
        if sol.time_ite<10
            sol.dt=sol.max_dt;
        else
            if not(sol.dt==sol.max_dt) && mod(sol.time_ite,10)==0
                sol.dt=min(sol.dt*8,sol.max_dt);
            end
        end

        if sol.coupling_scheme==0
            %% Solve concentration in the solid particles.
            run concentration_solid;
            
            %% Solve concentration in the electrolyte
            run concentration_electrolyte;

        elseif sol.coupling_scheme==1
            %% Solve the eqautions of the DFN model (potential in the solid and electrolyte (ps and pe) and concentration in the electrolyte and in the solid (ce, cs)).
            DFN_equations_coupled_solver();
        end

        sol.time=sol.time+sol.dt;
        sol.time_array(sol.time_ite) = sol.time;

        %% Calculate voltage at time t
        [fv.V,fv.SOC_neg,fv.SOC_pos]=voltage_calc(fv.ps,p.R_collector_contact ,ex.I_array(sol.time_ite));

        %% Add field variables to history records for later animation and output data
        save_fv_to_hist(sol.time_ite_save)

        if deb.break_time_loop==1
            disp("The Newton solver has diverged.Time solver iterrupted. "+num2str(sol.dt))
            break
        end

        dV_dt=(hist.V(sol.time_ite_save)-hist.V(max(sol.time_ite_save-1,1)))/sol.dt;

        sol.newt_ite_chrono= toc(newt_timer)-newt_chrono;
        %disp("ite "+int2str(sol.time_ite)+" at time "+num2str(sol.time)+"s I="+num2str(ex.I_array(sol.time_ite))+ ...
        %     "A CVed in "+num2str(sol.newt_ite)+"ite "+num2str(sol.newt_ite_chrono)+"s voltage="+num2str(fv.V)+ ...
        %     "V SOCn="+num2str(fv.SOC_neg)+" SOCp="+num2str(fv.SOC_pos)+" dt="+num2str(sol.dt)+ ...
        %     " dV/dt="+num2str(dV_dt))
        cprintf('*key',"ite "+int2str(sol.time_ite)+" at time "+num2str(sol.time)+"s I="+num2str(ex.I_array(sol.time_ite))+ ...
                "A CVed in "+num2str(sol.newt_ite)+"ite "+num2str(sol.newt_ite_chrono)+"s voltage="+num2str(fv.V)+ ...
                "V SOCn="+num2str(fv.SOC_neg)+" SOCp="+num2str(fv.SOC_pos)+" dt="+num2str(sol.dt)+ ...
                " dV/dt="+num2str(dV_dt)+"\n")

        disp("      avg(ce)="+num2str(mean(fv.ce))+" avg(csn)="+num2str(mean(mean(fv.csn)))+" avg(csp)="+num2str(mean(mean(fv.csp)))+" memory="+num2str(user.MemUsedMATLAB))

        %if fv.V<sol.min_allowed_voltage
        if abs(dV_dt)>sol.max_allowed_voltage_time_differential && fv.SOC_neg<0.2
            ex.I_array(sol.time_ite+1:length(sol.time_array)) = 0 ;
            sol.dt=sol.max_dt;
        end

        newt_chrono=toc(newt_timer);
        
    end
    disp('Solver ends');
end

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

        max_dt=sol.max_dt;
        if ((ex.I_array(sol.time_ite)>0 && ((min(fv.SOC_neg,fv.SOC_pos)<0.005) || ((fv.V - sol.min_allowed_voltage)<0.05)))) ...
        || ((ex.I_array(sol.time_ite)<0 && ((max(fv.SOC_neg,fv.SOC_pos)>0.995) || ((sol.max_allowed_voltage - fv.V)<0.05))))
            max_dt=sol.quick_dt;
        end
        disp((sol.time_ite))
        if sol.time_ite<=20
            sol.dt=sol.quick_dt;
        end
        if sol.time_ite==21
            sol.dt=max_dt;   
        end     
        
        if sol.time_ite>21
            if sol.dt>max_dt
                sol.dt=max_dt;
            elseif not(sol.dt==max_dt) && mod(sol.time_ite,20)==0
                if ex.I_array(sol.time_ite+1) == 0
                    sol.dt=min(max_dt);
                else
                    sol.dt=min(sol.dt*8,max_dt);
                end
            end
        end

        %% Solve the eqautions of the DFN model (potential in the solid and electrolyte (ps and pe) and concentration in the electrolyte and in the solid (ce, cs)).
        DFN_equations_coupled_solver();

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

        dV_dt=(hist.V(sol.time_ite_save+1)-hist.V(max(sol.time_ite_save,1)))/sol.dt;

        sol.newt_ite_chrono= toc(newt_timer)-newt_chrono;
        %disp("ite "+int2str(sol.time_ite)+" at time "+num2str(sol.time)+"s I="+num2str(ex.I_array(sol.time_ite))+ ...
        %     "A CVed in "+num2str(sol.newt_ite)+"ite "+num2str(sol.newt_ite_chrono)+"s voltage="+num2str(fv.V)+ ...
        %     "V SOCn="+num2str(fv.SOC_neg)+" SOCp="+num2str(fv.SOC_pos)+" dt="+num2str(sol.dt)+ ...
        %     " dV/dt="+num2str(dV_dt))
        cprintf('*key',"ite "+int2str(sol.time_ite)+" at time "+num2str(sol.time)+"s I="+num2str(ex.I_array(sol.time_ite))+ ...
                "A CVed in "+num2str(sol.newt_ite)+"ite "+num2str(sol.newt_ite_chrono)+"s voltage="+num2str(fv.V)+ ...
                "V SOCn="+num2str(fv.SOC_neg)+" SOCp="+num2str(fv.SOC_pos)+" dt="+num2str(sol.dt)+ ...
                " dV/dt="+num2str(dV_dt)+"\n")

        disp("   av(ce)="+num2str(mean(fv.ce))+" av(csn)="+num2str(mean(mean(fv.csn)))+" av(csp)="+num2str(mean(mean(fv.csp)))+ ...
            " av(pe)="+num2str(mean(fv.pe))+" av(psn)="+num2str(mean(mean(fv.ps(1:sol.nb_cell_n))))+" av(psp)="+num2str(mean(mean(fv.ps(sol.nb_cell_n+1:sol.nb_cell_ps))))+...
            " av(jn)="+num2str(mean(mean(fv.j(1:sol.nb_cell_n))))+" av(jp)="+num2str(mean(mean(fv.j(sol.nb_cell-sol.nb_cell_p+1:sol.nb_cell))))+" mi(csn)="+num2str(min(min(fv.csn))))


        if deb.timing_jacobian==1
            disp(" memory="+num2str(user.MemUsedMATLAB) ...
                +" PreNewt="+num2str(deb.chrono_newtsol_setup_singleite)+"s MatInv="+num2str(deb.chrono_matrix_inversion_singleite) ...
                +"s PostNewt="+num2str(deb.chrono_newtsol_update_singleite) )

            disp("      updt dtso="+num2str(deb.chrono_updt_dtso_singleite)+"s Calc Jac="+num2str(deb.chrono_Calc_Jac_singleite) ...
                +"s Calc f="+num2str(deb.chrono_Calc_f_singleite) +"s Calc gJacg="+num2str(deb.chrono_Calc_gJacg_singleite) ...
                +"s assemble="+num2str(deb.chrono_assemble_coupled_singleite)  )
            
            disp("      jacdiag="+num2str(deb.chrono_Jacdiag_singleite)+"s jacce="+num2str(deb.chrono_Jacce_singleite) ...
                +"s jacpe="+num2str(deb.chrono_Jacpe_singleite) +"s jacps="+num2str(deb.chrono_Jacps_singleite) ...
                +"s jaccsn="+num2str(deb.chrono_Jaccsn_singleite) +"s jaccsp="+num2str(deb.chrono_Jaccsp_singleite)  +"s")


            disp("      jaccedce="+num2str(deb.chrono_Jaccedce_singleite)+" jacpedpe="+num2str(deb.chrono_Jacpedpe_singleite) ...
                +"s jacpsdps="+num2str(deb.chrono_Jacpsdps_singleite)+"s jaccsdcs="+num2str(deb.chrono_Jaccsdcs_singleite) ...
                +"s total="+num2str(deb.chrono_Jaccedce_singleite+deb.chrono_Jacpedpe_singleite+deb.chrono_Jacpsdps_singleite+deb.chrono_Jaccsdcs_singleite) )
        end

        deb.chrono_newtsol_setup = deb.chrono_newtsol_setup + deb.chrono_newtsol_setup_singleite;
        deb.chrono_matrix_inversion = deb.chrono_matrix_inversion + deb.chrono_matrix_inversion_singleite;
        deb.chrono_newtsol_update = deb.chrono_newtsol_update + deb.chrono_newtsol_update_singleite;
        deb.chrono_newtsol_setup_singleite=0;
        deb.chrono_matrix_inversion_singleite=0;
        deb.chrono_newtsol_update_singleite=0;

        deb.chrono_updt_dtso_singleite=0;
        deb.chrono_Calc_Jac_singleite=0;
        deb.chrono_Calc_f_singleite=0;
        deb.chrono_Calc_gJacg_singleite=0;
        deb.chrono_assemble_coupled_singleite=0;

        deb.chrono_Jacdiag_singleite=0;
        deb.chrono_Jacce_singleite=0;
        deb.chrono_Jacpe_singleite=0;
        deb.chrono_Jacps_singleite=0;
        deb.chrono_Jaccsn_singleite=0;
        deb.chrono_Jaccsp_singleite=0;

        if deb.timing_jacobian==1
            deb.chrono_Jaccedce_singleite=0;
            deb.chrono_Jacpedpe_singleite=0;
            deb.chrono_Jacpsdps_singleite=0;
            deb.chrono_Jaccsdcs_singleite=0;
        end

        if ((dV_dt<-sol.max_allowed_voltage_time_differential && fv.SOC_neg<0.2) || (dV_dt>sol.max_allowed_voltage_time_differential && fv.SOC_neg>0.8) ...
                    || ((min(fv.SOC_neg,fv.SOC_pos)<-0.00000000001 || max(fv.SOC_neg,fv.SOC_pos)>1.00000000001) ) ...
                    || (fv.V>sol.max_allowed_voltage || fv.V<sol.min_allowed_voltage)) ...
                    && not(ex.I_array(sol.time_ite)==0)

            cprintf('err'," Charge complete. Current intensity reduced to 0.\n")

            ex.I_array(sol.time_ite+1:length(sol.time_array)) = 0 ;
            hist.charge_time = sol.time ;
            hist.charge_ite  = sol.time_ite ;
            sol.dt=sol.max_dt;
        end

        newt_chrono=toc(newt_timer);
        
    end
    
    disp("Solver ends. Estimated charge time "+num2str(hist.charge_time)+" at ite "+num2str(hist.charge_ite));
end
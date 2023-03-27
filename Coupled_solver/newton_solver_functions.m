classdef newton_solver_functions
    methods

        function M = calculate_mass_matrix_electrolyte(obj,eps)
            M1=diag(eps);
            M=M1;
        end

        function M = calculate_mass_matrix(obj,max_ind,Volumes)
            M1=diag(3.0/4.0*ones(max_ind,1)) + diag(1.0/8.0*ones(max_ind-1,1),1) + diag(1.0/8.0*ones(max_ind-1,1),-1);
            M1(1,2)=1/4;
            M1(max_ind,max_ind-1)=1/4;
            M2=diag(Volumes);
            M=M1*M2;
        end

        function [ce_next,ce_next_save,pe_next,pe_next_save,ps_next,ps_next_save,csn_next,csp_next,norm_delta_coupled,newt_ite] = ...
                                            Newton_solver_coupled(obj,csn_next,csp_next,rn,...
                                                                rp,j,csn,csp,cse, ...
                                                                Mn,Mp,ce_next,pe_next,ps_next,...
                                                                kappa, sigma,dx,source_ce,...
                                                                source_csn,source_csp,source_pe,source_ps,current_source_psBC,...
                                                                ce,separator_index,newt_max_ite,M,...
                                                                newt_lim,Jac_method,relax_factor)
            % Declare global variables
            global sol
            global eq_build_fun
            global deb
            global p
            global fv
            global ini
            global hist
            global BV_fun


            norm_delta_coupled=0;

            % Set up length of various vectors
            len=length(ce_next);
            len_ps=length(ps_next);
            lenr=sol.part_nb_cell+1;
            size_csn=size(csn);
            size_csp=size(csp);

            % Resize solid concentration vecotrs (from a table of size (nb_cell_n+nb_cell_p,part_nb_cell) to (nb_cell_n+nb_cell_p)*part_nb_cell)
            
            %test_array=[1 5 9; 2 6 10 ; 3 7 11 ; 4 8 12];
            %resized_test     = reshape(test_array,sol.nb_cell_n*(sol.part_nb_cell+1),1);
            resized_csn       = reshape(csn,sol.nb_cell_n*(sol.part_nb_cell+1),1);
            resized_csn_next  = reshape(csn,sol.nb_cell_n*(sol.part_nb_cell+1),1);
            resized_csp       = reshape(csp,sol.nb_cell_p*(sol.part_nb_cell+1),1);
            resized_csp_next  = reshape(csp,sol.nb_cell_p*(sol.part_nb_cell+1),1);
            cse_next  = cse;
            
            if deb.prints<=0
                hist.residuals   = zeros(6,sol.newton_meth_max_ite);
            end

            %Update parameter functions 
            
            [p.De_eff,kappa,p.kappa_D_eff] = eq_build_fun.update_param_functions(ce,cse,resized_csn,resized_csp);
            

            %calculate initial right hand side vectors (necessary for eqution with a time derivative (ie. for cs and ce))
            f_csn= eq_build_fun.LHS_f_cs(resized_csn,p.Dsn_array,rn,source_csn);
            %f_csn_copy= eq_build_fun.LHS_f_cs(resized_csn,Dsn,rn,source_csn);
            f_csp= eq_build_fun.LHS_f_cs(resized_csp,p.Dsp_array,rp,source_csp);
            f_ce = eq_build_fun.LHS_f_ce(ce,p.De_eff,dx,source_ce,"f_ce");
            f_ps = eq_build_fun.LHS_f_ps(fv.ps,ce,sigma,dx,source_ps,separator_index,current_source_psBC,"f_ps",0);
            f_pe = eq_build_fun.LHS_f_pe(fv.pe,ce,kappa,p.kappa_D_eff,dx,source_pe,"f_pe");

            ce_next_save= ones(length(ce_next),newt_max_ite);
            pe_next_save= ones(length(pe_next),newt_max_ite);
            ps_next_save= ones(length(ps_next),newt_max_ite);
            %% Start the newton method to solve g(c(t+dt))=0  (g(c(t+dt))=M(c(t+dt)-c(t))/dt -1/2 (f(t+dt)+f(t)))


            if deb.prints>5
                disp("DEBUG BEN DFN_equations_coupled_solver j and sources for ps, pe and ce newt "+num2str(sol.time_ite))
                disp(transpose(fv.Ueq))
                disp(fv.csn)
                disp(transpose(fv.ps))
                disp(fv.cse)
            end

            for newt_ite = 1:1:newt_max_ite

                newt_solv_setup_timer=tic;
                update_dt_sources_timer=tic;

                %if newt_ite>=40 && mod(newt_ite,20)==0
                %    if sol.dt>=2
                %        sol.dt=floor(sol.dt/2);
                %    else
                %        sol.dt=sol.dt/2;
                %    end
                %end

                %Update parameter functions and Butler-Volmer equation
                
                [p.De_eff,kappa,p.kappa_D_eff] = eq_build_fun.update_param_functions(ce_next,cse_next,resized_csn_next,resized_csp_next);
                fv.j=BV_fun.butler_volmer_equation(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,"newton_solver_functions next ite");
                [source_pe,source_ps,source_ce,source_csn,source_csp] = eq_build_fun.update_sources();
                

                deb.chrono_updt_dtso_singleite = deb.chrono_updt_dtso_singleite + toc(update_dt_sources_timer);
                Calc_Jac_timer=tic;

                %%Calculate the jacobian matrix of the left hand side vector for the concentration at time t+dt: Jf(c)
                % The jacobian is calculated using finite differences.
                
                Calc_Jacdiag_timer=tic;
                if deb.timing_jacobian==1
                    Calc_Jaccedce_timer=tic;
                end            
    
                Jac_f_ce_next = eq_build_fun.LHS_Jac_f_Fdiff_ce(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,p.De_eff,dx,source_ce,1);

                if deb.timing_jacobian==1
                    deb.chrono_Jaccedce_singleite = deb.chrono_Jaccedce_singleite + toc(Calc_Jaccedce_timer);
                    Calc_Jacpedpe_timer=tic;
                end

                Jac_f_pe_next = eq_build_fun.LHS_Jac_f_Fdiff_pe(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,kappa,p.kappa_D_eff,dx,source_pe);

                if deb.timing_jacobian==1
                    deb.chrono_Jacpedpe_singleite = deb.chrono_Jacpedpe_singleite + toc(Calc_Jacpedpe_timer);
                    Calc_Jacpsdps_timer=tic;
                end

                Jac_f_ps_next = eq_build_fun.LHS_Jac_f_Fdiff_ps(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,sigma,dx,source_ps,separator_index,current_source_psBC);

                if deb.timing_jacobian==1
                    deb.chrono_Jacpsdps_singleite = deb.chrono_Jacpsdps_singleite + toc(Calc_Jacpsdps_timer);
                    Calc_Jaccsdcs_timer=tic;
                end

                Jac_f_csn_next= eq_build_fun.LHS_Jac_f_Fdiff_cs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csn_next,p.Dsn_array,rn,source_csn,"csn",1);
                Jac_f_csp_next= eq_build_fun.LHS_Jac_f_Fdiff_cs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csp_next,p.Dsp_array,rp,source_csp,"csp",1);

                if deb.timing_jacobian==1
                    deb.chrono_Jaccsdcs_singleite = deb.chrono_Jaccsdcs_singleite + toc(Calc_Jaccsdcs_timer);
                end
                deb.chrono_Jacdiag_singleite = deb.chrono_Jacdiag_singleite + toc(Calc_Jacdiag_timer);

                Calc_Jacce_timer=tic;
                Jac_f_ce_dpe_next = eq_build_fun.LHS_Jac_f_Fdiff_cedpe(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,p.De_eff,dx,source_ce);
                Jac_f_ce_dps_next = eq_build_fun.LHS_Jac_f_Fdiff_cedps(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,p.De_eff,dx,source_ce);
                Jac_f_ce_dcs_next = eq_build_fun.LHS_Jac_f_Fdiff_cedcs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,p.De_eff,dx,source_ce,1);
                deb.chrono_Jacce_singleite = deb.chrono_Jacce_singleite + toc(Calc_Jacce_timer);

                Calc_Jacpe_timer=tic;
                Jac_f_pe_dce_next = eq_build_fun.LHS_Jac_f_Fdiff_pedce(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,kappa,p.kappa_D_eff,dx,source_pe,1);
                Jac_f_pe_dps_next = eq_build_fun.LHS_Jac_f_Fdiff_pedps(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,kappa,p.kappa_D_eff,dx,source_pe);
                Jac_f_pe_dcs_next = eq_build_fun.LHS_Jac_f_Fdiff_pedcs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,kappa,p.kappa_D_eff,dx,source_pe,1);
                deb.chrono_Jacpe_singleite = deb.chrono_Jacpe_singleite + toc(Calc_Jacpe_timer);

                Calc_Jacps_timer=tic;
                Jac_f_ps_dce_next = eq_build_fun.LHS_Jac_f_Fdiff_psdce(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,sigma,dx,source_ps,separator_index,current_source_psBC);
                Jac_f_ps_dpe_next = eq_build_fun.LHS_Jac_f_Fdiff_psdpe(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,sigma,dx,source_ps,separator_index,current_source_psBC);
                Jac_f_ps_dcs_next = eq_build_fun.LHS_Jac_f_Fdiff_psdcs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,sigma,dx,source_ps,separator_index,current_source_psBC,1);
                deb.chrono_Jacps_singleite = deb.chrono_Jacps_singleite + toc(Calc_Jacps_timer);

                Calc_Jaccsn_timer=tic;
                Jac_f_csn_dce_next = eq_build_fun.LHS_Jac_f_Fdiff_csdce(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csn_next,p.Dsn_array,rn,source_csn,"csn");
                Jac_f_csn_dpe_next = eq_build_fun.LHS_Jac_f_Fdiff_csdpe(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csn_next,p.Dsn_array,rn,source_csn,"csn");
                Jac_f_csn_dps_next = eq_build_fun.LHS_Jac_f_Fdiff_csdps(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csn_next,p.Dsn_array,rn,source_csn,"csn");
                deb.chrono_Jaccsn_singleite = deb.chrono_Jaccsn_singleite + toc(Calc_Jaccsn_timer);

                Calc_Jaccsp_timer=tic;
                Jac_f_csp_dce_next = eq_build_fun.LHS_Jac_f_Fdiff_csdce(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csp_next,p.Dsp_array,rp,source_csp,"csp");
                Jac_f_csp_dpe_next = eq_build_fun.LHS_Jac_f_Fdiff_csdpe(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csp_next,p.Dsp_array,rp,source_csp,"csp");
                Jac_f_csp_dps_next = eq_build_fun.LHS_Jac_f_Fdiff_csdps(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csp_next,p.Dsp_array,rp,source_csp,"csp");
                deb.chrono_Jaccsp_singleite = deb.chrono_Jaccsp_singleite + toc(Calc_Jaccsp_timer);


                    

                deb.chrono_Calc_Jac_singleite = deb.chrono_Calc_Jac_singleite + toc(Calc_Jac_timer);
                Calc_f_timer=tic;

                f_ce_next       = eq_build_fun.LHS_f_ce(ce_next,p.De_eff,dx,source_ce,"f_ce_next");
                f_pe_next       = eq_build_fun.LHS_f_pe(pe_next,ce_next,kappa,p.kappa_D_eff,dx,source_pe,"f_pe_next");
                f_ps_next       = eq_build_fun.LHS_f_ps(ps_next,ce_next,sigma,dx,source_ps,separator_index,current_source_psBC,"f_ps_next",0);
                f_csn_next      = eq_build_fun.LHS_f_cs(resized_csn_next,p.Dsn_array,rn,source_csn);
                f_csp_next      = eq_build_fun.LHS_f_cs(resized_csp_next,p.Dsp_array,rp,source_csp);
                
                deb.chrono_Calc_f_singleite = deb.chrono_Calc_f_singleite + toc(Calc_f_timer);
                Calc_gJacg_timer=tic;


                %%Compilation of the Crank Nicholson method s function and jacobian
                %g(c(t+dt))=M(c(t+dt)-c(t))/dt -1/2 (f(t+dt)+f(t))
                %Jac_g(c(t+dt))=M/dt -1/2 jac_f(t+dt)
                %Jac_g*(delta(c(t+dt)))=g
                

                %Solve Li concentration in the solid
                gcsn=Mn*(resized_csn_next-resized_csn)/sol.dt - 0.5*(f_csn_next+f_csn);
                gcsp=Mp*(resized_csp_next-resized_csp)/sol.dt - 0.5*(f_csp_next+f_csp);
                Jac_gcsn=Mn/sol.dt -0.5*Jac_f_csn_next;
                Jac_gcsp=Mp/sol.dt -0.5*Jac_f_csp_next;

                Jac_gcsn_dpe=-0.5*Jac_f_csn_dpe_next;
                Jac_gcsp_dpe=-0.5*Jac_f_csp_dpe_next;
                Jac_gcsn_dce=-0.5*Jac_f_csn_dce_next;
                Jac_gcsp_dce=-0.5*Jac_f_csp_dce_next;
                Jac_gcsn_dps=-0.5*Jac_f_csn_dps_next;
                Jac_gcsp_dps=-0.5*Jac_f_csp_dps_next;

                %Solve Li concentration in electrolyte
                gc=M*(ce_next-ce)/sol.dt - 0.5*(f_ce_next+f_ce);
                Jac_gce=M/sol.dt -0.5*Jac_f_ce_next;
                Jac_gce_dpe= -0.5*Jac_f_ce_dpe_next;
                Jac_gce_dps= -0.5*Jac_f_ce_dps_next;
                Jac_gce_dcs= -0.5*Jac_f_ce_dcs_next;


                %Solve electric potential in electrolyte
                gp = f_pe_next;
                Jac_gp= Jac_f_pe_next;
                Jac_gpe_dce= Jac_f_pe_dce_next;
                Jac_gpe_dps= Jac_f_pe_dps_next;
                Jac_gpe_dcs= Jac_f_pe_dcs_next;

                %Solve electric potential in solid
                gps = f_ps_next;
                Jac_gps= Jac_f_ps_next;
                Jac_gps_dce= Jac_f_ps_dce_next;
                Jac_gps_dpe= Jac_f_ps_dpe_next;
                Jac_gps_dcs= Jac_f_ps_dcs_next;

                deb.chrono_Calc_gJacg_singleite = deb.chrono_Calc_gJacg_singleite + toc(Calc_gJacg_timer);


                % Form the coupled system 

                assemble_coupled_timer=tic;

                Jac_coupled = (blkdiag(Jac_gcsn,Jac_gcsp,Jac_gce,Jac_gp,Jac_gps)) ;

                begin_csn =1;
                end_csn   =(sol.part_nb_cell+1)*(sol.nb_cell_n);
                begin_csp =(sol.part_nb_cell+1)*(sol.nb_cell_n)+1;
                end_csp   =(sol.part_nb_cell+1)*(sol.nb_cell_n+sol.nb_cell_p);
                begin_cs  =1;
                end_cs    =(sol.part_nb_cell+1)*(sol.nb_cell_n+sol.nb_cell_p);
                begin_ce  =(sol.part_nb_cell+1)*(sol.nb_cell_n+sol.nb_cell_p)+1;
                end_ce    =(sol.part_nb_cell+1)*(sol.nb_cell_n+sol.nb_cell_p)+sol.nb_cell;
                begin_pe  =(sol.part_nb_cell+1)*(sol.nb_cell_n+sol.nb_cell_p)+sol.nb_cell+1;
                end_pe    =(sol.part_nb_cell+1)*(sol.nb_cell_n+sol.nb_cell_p)+sol.nb_cell+sol.nb_cell;
                begin_ps  =(sol.part_nb_cell+1)*(sol.nb_cell_n+sol.nb_cell_p)+sol.nb_cell+sol.nb_cell+1;
                end_ps    =(sol.part_nb_cell+1)*(sol.nb_cell_n+sol.nb_cell_p)+sol.nb_cell+sol.nb_cell+sol.nb_cell_n+sol.nb_cell_p;
                
                

                Jac_coupled(begin_ce:end_ce,begin_pe:end_pe)=Jac_gce_dpe;
                Jac_coupled(begin_ce:end_ce,begin_ps:end_ps)=Jac_gce_dps;
                Jac_coupled(begin_ce:end_ce,begin_cs:end_cs)=Jac_gce_dcs;

                Jac_coupled(begin_pe:end_pe,begin_ce:end_ce)=Jac_gpe_dce;
                Jac_coupled(begin_pe:end_pe,begin_ps:end_ps)=Jac_gpe_dps;
                Jac_coupled(begin_pe:end_pe,begin_cs:end_cs)=Jac_gpe_dcs;

                Jac_coupled(begin_ps:end_ps,begin_ce:end_ce)=Jac_gps_dce;
                Jac_coupled(begin_ps:end_ps,begin_pe:end_pe)=Jac_gps_dpe;
                Jac_coupled(begin_ps:end_ps,begin_cs:end_cs)=Jac_gps_dcs;

                Jac_coupled(begin_csn:end_csn,begin_ce:end_ce)=Jac_gcsn_dce;
                Jac_coupled(begin_csn:end_csn,begin_pe:end_pe)=Jac_gcsn_dpe;
                Jac_coupled(begin_csn:end_csn,begin_ps:end_ps)=Jac_gcsn_dps;

                Jac_coupled(begin_csp:end_csp,begin_ce:end_ce)=Jac_gcsp_dce;
                Jac_coupled(begin_csp:end_csp,begin_pe:end_pe)=Jac_gcsp_dpe;
                Jac_coupled(begin_csp:end_csp,begin_ps:end_ps)=Jac_gcsp_dps;

                Jac_coupled=sparse(Jac_coupled);

                %figure(12454545);
                %spy(Jac_coupled)

                g_coupled   = cat(1,gcsn,gcsp);
                g_coupled   = cat(1,g_coupled,gc);
                g_coupled   = cat(1,g_coupled,gp);
                g_coupled   = cat(1,g_coupled,gps);

                deb.chrono_assemble_coupled_singleite = deb.chrono_assemble_coupled_singleite + toc(assemble_coupled_timer);
                deb.chrono_newtsol_setup_singleite = deb.chrono_newtsol_setup_singleite + toc(newt_solv_setup_timer);



                % Solve the system

                mat_inv_timer=tic;

                Newton_method_version=1;
                if Newton_method_version==1 || (Newton_method_version==2 && newt_ite==1) % Original Newton method
                    delta_coupled = (Jac_coupled) \ (-g_coupled);

                elseif Newton_method_version==2  % McDougall 2019
                    delta_interm = (Jac_coupled) \ (-g_coupled);

                    Aest2_up    = 3*(g_coupled_m1 - g_coupled) * transpose(- ones(length(delta_coupled),1) ./ delta_coupled) - Jac_coupled_m1 - 2*Jac_coupled;
                    Aest2_down  = transpose(-delta_coupled) *Jac_coupled;
                    Aest2= Aest2_up * transpose(ones(1,length(Aest2_down))./Aest2_down) ;

                    F=0.5 + sqrt(max(0 , 0.25 - abs(dot(Aest2,transpose(delta_interm)))));
                    delta_coupled = delta_interm/F;
                end

                if Newton_method_version==2
                    g_coupled_m1    = g_coupled;
                    Jac_coupled_m1  = Jac_coupled;
                end

                deb.chrono_matrix_inversion_singleite = deb.chrono_matrix_inversion_singleite + toc(mat_inv_timer);



                % Update solutions, record history data and measure resudual

                newtsol_update_timer=tic;

                norm_delta_coupled = sqrt(sum(transpose(delta_coupled) * delta_coupled))/length(g_coupled);
                norm_delta_ps = sqrt(sum(transpose(delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len+len_ps)) * delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len+len_ps)))/len_ps;
                norm_delta_pe = sqrt(sum(transpose(delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len)) * delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len)))/len;
                norm_delta_ce = sqrt(sum(transpose(delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len)) * delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len)))/len;
                norm_delta_csn = sqrt(sum(transpose(delta_coupled(1 : sol.nb_cell_n*lenr)) * delta_coupled(1 : sol.nb_cell_n*lenr)))/(sol.nb_cell_n*lenr);
                norm_delta_csp = sqrt(sum(transpose(delta_coupled(sol.nb_cell_n*lenr+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr)) * delta_coupled(sol.nb_cell_n*lenr+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr)))/(sol.nb_cell_p*lenr);

                hist.residuals(1,newt_ite)=norm_delta_coupled;
                hist.residuals(2,newt_ite)=norm_delta_ps;
                hist.residuals(3,newt_ite)=norm_delta_pe;
                hist.residuals(4,newt_ite)=norm_delta_ce;
                hist.residuals(5,newt_ite)=norm_delta_csn;
                hist.residuals(6,newt_ite)=norm_delta_csp;

                if (deb.prints>=5) && (isnan(norm_delta_coupled)==1 | isreal(norm_delta_coupled)==0)
                    

                    disp("DEBUG BEN j ------------------------------ newt_ite="+num2str(newt_ite)+" , time ite="+num2str(sol.time_ite))
                    disp((fv.j))
                    pe_next_reduced = cat(1,pe_next(1:sol.nb_cell_n),pe_next(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell));
                    disp(transpose(ps_next-pe_next_reduced))


                    disp("DEBUG BEN cs newt_ite="+num2str(newt_ite)+" , time ite="+num2str(sol.time_ite))

                    disp(transpose(gcsn))
                    disp(transpose(resized_csn_next))
                    disp(transpose(resized_csp_next))
                    %disp((Jac_gcsn))
                    %disp((Jac_gcsp))
                    disp(transpose(delta_coupled(1 : sol.nb_cell_n*lenr)))
                    disp(transpose(delta_coupled(sol.nb_cell_n*lenr+1:sol.nb_cell_n*lenr+sol.nb_cell_p*lenr)))
                    
                    %figure(265);
                    %plot(1:sol.nb_cell_n*lenr,transpose(delta_coupled(1 : sol.nb_cell_n*lenr)),'LineWidth',2);
                    %grid on
                    %grid minor

                    

                    disp("DEBUG BEN ce newt_ite="+num2str(newt_ite)+" , time ite="+num2str(sol.time_ite))
                    disp(transpose(gc))
                    disp(transpose(ce_next))
                    disp((Jac_gce))
                    disp(transpose(delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len)))


                    disp("DEBUG BEN pe newt_ite="+num2str(newt_ite)+" , time ite="+num2str(sol.time_ite))

                    disp(transpose(gp))
                    disp(transpose(pe_next))
                    disp((Jac_gp))  
                    disp(transpose(delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len)))
                    
                    disp("DEBUG BEN ps newt_ite="+num2str(newt_ite)+" , time ite="+num2str(sol.time_ite))
                    disp(transpose(gps))
                    disp(source_ps)
                    disp(transpose(ps_next))
                    disp((Jac_gps))
                    disp(transpose(delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len+len_ps)))
                    

                    disp("DEBUG BEN ******************************")
                end

                
                % Update solutions 

                if isnan(norm_delta_coupled)==1 | isreal(norm_delta_coupled)==0
                    disp("Newton method for for coupled system did not converge: "+num2str(newt_ite)+" , "+num2str(sol.time_ite)+" , "+num2str(norm_delta_coupled)+" , "+num2str(det(Jac_coupled))+" , "+num2str(sol.dt))
                    newt_ite=newt_ite-1;
                    break
                end

                plateau=1;
                end_val=sol.nb_cell_n*lenr;
                resized_csn_next = resized_csn_next + relax_factor*delta_coupled(plateau:end_val);

                plateau=end_val+1;
                end_val = end_val + sol.nb_cell_p*lenr;
                resized_csp_next = resized_csp_next + relax_factor*delta_coupled(plateau:end_val);

                plateau=end_val+1;
                end_val = end_val + len;
                ce_next = ce_next + relax_factor*delta_coupled(plateau:end_val);

                plateau=end_val+1;
                end_val = end_val + len;
                pe_next = pe_next + relax_factor*delta_coupled(plateau:end_val);

                plateau=end_val+1;
                end_val = end_val + len_ps;
                ps_next = ps_next + relax_factor*delta_coupled(plateau:end_val);



                % Bring variables back to normal levels

                min_c=0.00000;
                
                resized_csn_next( resized_csn_next <= min_c ) = min_c;
                resized_csp_next( resized_csp_next <= min_c ) = min_c;
                ce_next( ce_next <= min_c ) = min_c;



                % Recompose cse (concentration at the surface of each particule.)

                for iiii =[1:sol.nb_cell_n]
                    cse_next(iiii)=resized_csn_next((sol.part_nb_cell+1)*iiii);
                end
                for iiii =[1:sol.nb_cell_p]
                    cse_next(iiii+sol.nb_cell_n)=resized_csp_next((sol.part_nb_cell+1)*(iiii));
                end

                %if deb.videos_generation==1
                %    if norm_delta_coupled>newt_lim
                %        ce_next_save(:,newt_ite)= ce_next;
                %        pe_next_save(:,newt_ite)= pe_next;
                %        ps_next_save(:,newt_ite)= ps_next;
                %    end 
                %end

                % Determine if solver must be stopped and display convergence data

                if deb.prints>=2
                    
                    disp("newt_ite="+num2str(newt_ite)+" sol.time_ite="+num2str(sol.time_ite)+" residuals: coupled="+num2str(norm_delta_coupled)+" , ps="+num2str(norm_delta_ps)+" , pe="+num2str(norm_delta_pe) ...
                                                                             +" , ce="+num2str(norm_delta_ce)+" , csn="+num2str(norm_delta_csn)+" , csp="+num2str(norm_delta_csp))
                end

                if (norm_delta_coupled<newt_lim) && newt_ite>1
                    
                    if deb.prints>=0
                        disp("      Newton method for the coupled system converged:"+num2str(newt_ite)+" , "+num2str(norm_delta_coupled) ...
                                                                            +" , "+num2str(norm_delta_ps)+" , "+num2str(norm_delta_pe) ...
                                                                            +" , "+num2str(norm_delta_ce)+" , "+num2str(norm_delta_csn)+" , "+num2str(norm_delta_csp) ...
                                                                            +" , and "+num2str(-(log10(max(hist.residuals(1,:)))-log10(min(hist.residuals(1,1:newt_ite))))) ...
                                                                            +" , "+num2str(max(hist.residuals(1,:))) +" , "+num2str(min(hist.residuals(1,1:newt_ite))))
                    end
                    break
                elseif newt_ite==newt_max_ite
                    if deb.prints>=1
                        disp("      Newton method for the coupled system did not converge:"+num2str(newt_ite)+" , "+num2str(norm_delta_coupled) ...
                                                                            +" , "+num2str(norm_delta_ps)+" , "+num2str(norm_delta_pe) ...
                                                                            +" , "+num2str(norm_delta_ce)+" , "+num2str(norm_delta_csn)+" , "+num2str(norm_delta_csp))
                        disp("      Limit of the solver: "+num2str(log10(newt_lim))+" current order of magnitude of residual: "+ num2str(log10(norm_delta_coupled))+" det(Jac_coupled)="+num2str(det(Jac_coupled)) )
                    end
                    newt_ite=newt_ite+1;
                    break
                end
                

                deb.chrono_newtsol_update_singleite = deb.chrono_newtsol_update_singleite + toc(newtsol_update_timer);


            end

            if sol.time_ite==sol.nb_steps && deb.prints>=4  % Write the jacobian for debuging purposes
                writetable(array2table(Jac_coupled),'Jac_coupled_.xlsx')
            end

            hist.sourceps(:,sol.time_ite)=source_ps;
            hist.sourcepsBC(sol.time_ite)=current_source_psBC(2);

            % Save the end residuals 
            if deb.prints>=0    
                hist.residuals_time(1,sol.time_ite)=norm_delta_coupled;
                hist.residuals_time(2,sol.time_ite)=norm_delta_ps;
                hist.residuals_time(3,sol.time_ite)=norm_delta_pe;
                hist.residuals_time(4,sol.time_ite)=norm_delta_ce;
                hist.residuals_time(5,sol.time_ite)=norm_delta_csn;
                hist.residuals_time(6,sol.time_ite)=norm_delta_csp;
                hist.newt_it_number(1,sol.time_ite)=newt_ite;

                if isnan(norm_delta_coupled)==0
                    hist.residuals_diff(1,sol.time_ite)=-(log10(max(hist.residuals(1,:)))-log10(hist.residuals(1,min(newt_ite,newt_max_ite))));
                    hist.residuals_diff(2,sol.time_ite)=-(log10(max(hist.residuals(2,:)))-log10(hist.residuals(2,min(newt_ite,newt_max_ite))));
                    hist.residuals_diff(3,sol.time_ite)=-(log10(max(hist.residuals(3,:)))-log10(hist.residuals(3,min(newt_ite,newt_max_ite))));
                    hist.residuals_diff(4,sol.time_ite)=-(log10(max(hist.residuals(4,:)))-log10(hist.residuals(4,min(newt_ite,newt_max_ite))));
                    hist.residuals_diff(5,sol.time_ite)=-(log10(max(hist.residuals(5,:)))-log10(hist.residuals(5,min(newt_ite,newt_max_ite))));
                    hist.residuals_diff(6,sol.time_ite)=-(log10(max(hist.residuals(6,:)))-log10(hist.residuals(6,min(newt_ite,newt_max_ite))));
                end
            end

            %Recompose cs vectors (concentration in the solid particles)
            csn_next  = reshape(resized_csn_next,sol.part_nb_cell+1,sol.nb_cell_n);
            csp_next  = reshape(resized_csp_next,sol.part_nb_cell+1,sol.nb_cell_n);
            
        end
    end
end
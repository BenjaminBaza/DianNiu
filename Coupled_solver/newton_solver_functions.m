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
                                            Newton_solver_coupled(obj,csn_next,csp_next,Dsn,Dsp,rn,...
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
                %disp(sol.nb_cell_n*sol.part_nb_cell)
                %disp(size_csn)
                %disp(ce)
                %disp(csn)
                %disp(resized_csn)
            end

            %Update parameter functions 
            if sol.newton_update_sources==1
                [p.De_eff,kappa,p.kappa_D_eff] = eq_build_fun.update_param_functions(ce,cse);
            end


            %calculate initial right hand side vectors (necessary for eqution with a time derivative (ie. for cs and ce))
            f_csn= eq_build_fun.LHS_f_cs(resized_csn,Dsn,rn,source_csn);
            f_csp= eq_build_fun.LHS_f_cs(resized_csp,Dsp,rp,source_csp);
            f_ce = eq_build_fun.LHS_f_ce(ce,p.De_eff,dx,source_ce,"f_ce");
            f_ps = eq_build_fun.LHS_f_ps(fv.ps,ce,sigma,dx,source_ps,separator_index,current_source_psBC,"f_ps",0);
            f_pe = eq_build_fun.LHS_f_pe(fv.pe,ce,kappa,p.kappa_D_eff,dx,source_pe,"f_pe");

            ce_next_save= ones(length(ce_next),newt_max_ite);
            pe_next_save= ones(length(pe_next),newt_max_ite);
            ps_next_save= ones(length(ps_next),newt_max_ite);
            %% Start the newton method to solve g(c(t+dt))=0  (g(c(t+dt))=M(c(t+dt)-c(t))/dt -1/2 (f(t+dt)+f(t)))


            if deb.prints>1
                disp("DEBUG BEN DFN_equations_coupled_solver j and sources for ps, pe and ce newt "+num2str(sol.time_ite))
                disp(transpose(fv.Ueq))
                disp(fv.csn)
                disp(transpose(fv.ps))
                disp(fv.cse)
            end

            for newt_ite = 1:1:newt_max_ite

                newt_solv_setup_timer=tic;
                update_dt_sources_timer=tic;

                if newt_ite>=40 && mod(newt_ite,20)==0
                    if sol.dt>=2
                        sol.dt=floor(sol.dt/2);
                    else
                        sol.dt=sol.dt/2;
                    end
                end

                %Update parameter functions and Butler-Volmer equation
                if sol.newton_update_sources==1
                    [p.De_eff,kappa,p.kappa_D_eff] = eq_build_fun.update_param_functions(ce_next,cse_next);
                    fv.j=BV_fun.butler_volmer_equation(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,"newton_solver_functions next ite");
                    [source_pe,source_ps,source_ce,source_csn,source_csp] = eq_build_fun.update_sources();
                end

                deb.chrono_updt_dtso_singleite = deb.chrono_updt_dtso_singleite + toc(update_dt_sources_timer);
                Calc_Jac_timer=tic;

                %%Calculate the jacobian matrix of the left hand side vector for the concentration at time t+dt: Jf(c)
                %if method=2, jacobian is calculated analytically, if method=1, jac
                %is calculated using finite differences.
                if Jac_method==1
                    
                    Calc_Jacdiag_timer=tic;
                    Calc_Jaccedce_timer=tic;
                    Jac_f_ce_next = eq_build_fun.LHS_Jac_f_Fdiff_ce(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,p.De_eff,dx,source_ce);
                    deb.chrono_Jaccedce_singleite = deb.chrono_Jaccedce_singleite + toc(Calc_Jaccedce_timer);

                    Jac_f_pe_next = eq_build_fun.LHS_Jac_f_Fdiff_pe(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,kappa,p.kappa_D_eff,dx,source_pe);
                    Jac_f_ps_next = eq_build_fun.LHS_Jac_f_Fdiff_ps(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,sigma,dx,source_ps,separator_index,current_source_psBC);
                    Jac_f_csn_next= eq_build_fun.LHS_Jac_f_Fdiff_cs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csn_next,Dsn,rn,source_csn,"csn");
                    Jac_f_csp_next= eq_build_fun.LHS_Jac_f_Fdiff_cs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csp_next,Dsp,rp,source_csp,"csp");
                    deb.chrono_Jacdiag_singleite = deb.chrono_Jacdiag_singleite + toc(Calc_Jacdiag_timer);

                    Calc_Jacce_timer=tic;
                    Jac_f_ce_dpe_next = eq_build_fun.LHS_Jac_f_Fdiff_cedpe(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,p.De_eff,dx,source_ce);
                    Jac_f_ce_dps_next = eq_build_fun.LHS_Jac_f_Fdiff_cedps(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,p.De_eff,dx,source_ce);
                    Jac_f_ce_dcs_next = eq_build_fun.LHS_Jac_f_Fdiff_cedcs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,p.De_eff,dx,source_ce);
                    deb.chrono_Jacce_singleite = deb.chrono_Jacce_singleite + toc(Calc_Jacce_timer);

                    Calc_Jacpe_timer=tic;
                    Jac_f_pe_dce_next = eq_build_fun.LHS_Jac_f_Fdiff_pedce(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,kappa,p.kappa_D_eff,dx,source_pe);
                    Jac_f_pe_dps_next = eq_build_fun.LHS_Jac_f_Fdiff_pedps(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,kappa,p.kappa_D_eff,dx,source_pe);
                    Jac_f_pe_dcs_next = eq_build_fun.LHS_Jac_f_Fdiff_pedcs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,kappa,p.kappa_D_eff,dx,source_pe);
                    deb.chrono_Jacpe_singleite = deb.chrono_Jacpe_singleite + toc(Calc_Jacpe_timer);

                    Calc_Jacps_timer=tic;
                    Jac_f_ps_dce_next = eq_build_fun.LHS_Jac_f_Fdiff_psdce(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,sigma,dx,source_ps,separator_index,current_source_psBC);
                    Jac_f_ps_dpe_next = eq_build_fun.LHS_Jac_f_Fdiff_psdpe(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,sigma,dx,source_ps,separator_index,current_source_psBC);
                    Jac_f_ps_dcs_next = eq_build_fun.LHS_Jac_f_Fdiff_psdcs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,sigma,dx,source_ps,separator_index,current_source_psBC);
                    deb.chrono_Jacps_singleite = deb.chrono_Jacps_singleite + toc(Calc_Jacps_timer);

                    Calc_Jaccsn_timer=tic;
                    Jac_f_csn_dce_next = eq_build_fun.LHS_Jac_f_Fdiff_csdce(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csn_next,Dsn,rn,source_csn,"csn");
                    Jac_f_csn_dpe_next = eq_build_fun.LHS_Jac_f_Fdiff_csdpe(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csn_next,Dsn,rn,source_csn,"csn");
                    Jac_f_csn_dps_next = eq_build_fun.LHS_Jac_f_Fdiff_csdps(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csn_next,Dsn,rn,source_csn,"csn");
                    deb.chrono_Jaccsn_singleite = deb.chrono_Jaccsn_singleite + toc(Calc_Jaccsn_timer);

                    Calc_Jaccsp_timer=tic;
                    Jac_f_csp_dce_next = eq_build_fun.LHS_Jac_f_Fdiff_csdce(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csp_next,Dsp,rp,source_csp,"csp");
                    Jac_f_csp_dpe_next = eq_build_fun.LHS_Jac_f_Fdiff_csdpe(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csp_next,Dsp,rp,source_csp,"csp");
                    Jac_f_csp_dps_next = eq_build_fun.LHS_Jac_f_Fdiff_csdps(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csp_next,Dsp,rp,source_csp,"csp");
                    deb.chrono_Jaccsp_singleite = deb.chrono_Jaccsp_singleite + toc(Calc_Jaccsp_timer);


                    %if (1==0 | det(Jac_f_pe_next)==0)
                    %    if (deb.prints>0)
                    %        disp("The jacobian for the potential in the electrolyte has a determinant of 0. The diagonal boost method is used to make the matrix invertible.")
                    %    end
                    %    for iiii=1:1:len
                    %        Jac_f_pe_next(iiii,iiii)=Jac_f_pe_next(iiii,iiii)*1.0000000001;
                    %    end
                    %end
                    %if (1==0 | det(Jac_f_ps_next(1:sol.nb_cell_n,1:sol.nb_cell_n))==0)
                    %    if (deb.prints>0)
                    %        disp("The jacobian for the potential in the neg. electrode has a determinant of 0. The diagonal boost method is used to make the matrix invertible.")
                    %    end
                    %    for iiii=1:1:sol.nb_cell_n
                    %        Jac_f_ps_next(iiii,iiii)=Jac_f_ps_next(iiii,iiii)*1.0000000001;
                    %    end
                    %end
                    %if (1==0 | (det(Jac_f_ps_next(sol.nb_cell_n+1:sol.nb_cell_p+sol.nb_cell_n,sol.nb_cell_n+1:sol.nb_cell_p+sol.nb_cell_n))==0))
                    %    if (deb.prints>0)
                    %        disp("The jacobian for the potential in the pos. electrode has a determinant of 0. The diagonal boost method is used to make the matrix invertible.")
                    %    end
                    %    for iiii=sol.nb_cell_n+1:1:sol.nb_cell_p+sol.nb_cell_n
                    %        Jac_f_ps_next(iiii,iiii)=Jac_f_ps_next(iiii,iiii)*1.0000000001;
                    %    end
                    %end

                end

                deb.chrono_Calc_Jac_singleite = deb.chrono_Calc_Jac_singleite + toc(Calc_Jac_timer);
                Calc_f_timer=tic;

                f_ce_next       = eq_build_fun.LHS_f_ce(ce_next,p.De_eff,dx,source_ce,"f_ce_next");
                f_pe_next       = eq_build_fun.LHS_f_pe(pe_next,ce_next,kappa,p.kappa_D_eff,dx,source_pe,"f_pe_next");
                f_ps_next       = eq_build_fun.LHS_f_ps(ps_next,ce_next,sigma,dx,source_ps,separator_index,current_source_psBC,"f_ps_next",0);
                f_csn_next      = eq_build_fun.LHS_f_cs(resized_csn_next,Dsn,rn,source_csn);
                f_csp_next      = eq_build_fun.LHS_f_cs(resized_csp_next,Dsp,rp,source_csp);
                
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

                Jac_coupled = blkdiag(Jac_gcsn,Jac_gcsp,Jac_gce,Jac_gp,Jac_gps) ;

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
                %Jac_coupled(begin_ce:end_ce,begin_cs:end_cs)=Jac_gce_dcs;

                Jac_coupled(begin_pe:end_pe,begin_ce:end_ce)=Jac_gpe_dce;
                Jac_coupled(begin_pe:end_pe,begin_ps:end_ps)=Jac_gpe_dps;
                %Jac_coupled(begin_pe:end_pe,begin_cs:end_cs)=Jac_gpe_dcs;


                Jac_coupled(begin_ps:end_ps,begin_ce:end_ce)=Jac_gps_dce;
                Jac_coupled(begin_ps:end_ps,begin_pe:end_pe)=Jac_gps_dpe;
                %Jac_coupled(begin_ps:end_ps,begin_cs:end_cs)=Jac_gps_dcs;

                Jac_coupled(begin_csn:end_csn,begin_ce:end_ce)=Jac_gcsn_dce;
                Jac_coupled(begin_csn:end_csn,begin_pe:end_pe)=Jac_gcsn_dpe;
                Jac_coupled(begin_csn:end_csn,begin_ps:end_ps)=Jac_gcsn_dps;

                Jac_coupled(begin_csp:end_csp,begin_ce:end_ce)=Jac_gcsp_dce;
                Jac_coupled(begin_csp:end_csp,begin_pe:end_pe)=Jac_gcsp_dpe;
                Jac_coupled(begin_csp:end_csp,begin_ps:end_ps)=Jac_gcsp_dps;


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

                if (deb.prints>=2) && (isnan(norm_delta_coupled)==1 | isreal(norm_delta_coupled)==0)
                    

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
                    

                    %if newt_ite==4
                    %    break
                    %end

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


                plot_csp_local_newton_results=0;
                if plot_csp_local_newton_results==1
                    figure(2622222);
                    fs = 16;
                    set(gcf,'Position',[50 50 1800 1000]);     

                    subplot(5,5,[1,2])
                    plot(1:sol.nb_cell_p*lenr,transpose(delta_coupled(sol.nb_cell_n*lenr+1:sol.nb_cell_n*lenr+sol.nb_cell_p*lenr)),'LineWidth',2);
                    title('delta coupled','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[4,5])
                    plot(1:sol.nb_cell_p*lenr,resized_csp_next,'LineWidth',2);
                    title('csp next','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[6,7])
                    plot(1:sol.nb_cell_p,fv.j(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell),'LineWidth',2);
                    title('j','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[9,10])
                    plot(1:sol.nb_cell_p,source_csp,'LineWidth',2);
                    title('source csp','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[11,12])
                    plot(1:sol.nb_cell_n,cse_next(1:sol.nb_cell_n),'LineWidth',2);
                    title('csen','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[14,15])
                    plot(sol.nb_cell_n+1:2*sol.nb_cell_n,cse_next(sol.nb_cell_n+1:2*sol.nb_cell_n),'LineWidth',2);
                    title('csep','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[16,17])
                    plot(1:sol.nb_cell_p*lenr,gcsp,'LineWidth',2);
                    title('gcsp','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[19,20])
                    plot(1:sol.nb_cell_p*lenr,Mp*(resized_csp_next-resized_csp)/sol.dt,'LineWidth',2);
                    title('M*dcsp/dt','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[21,22])
                    plot(1:sol.nb_cell_p*lenr,- 0.5*(f_csp_next+f_csp),'LineWidth',2);
                    title('-fcsp','fontsize',fs);
                    grid on
                    grid minor
                end


                plot_csn_local_newton_results=0;
                if plot_csn_local_newton_results==1 && sol.time_ite>=19
                    figure(266);
                    fs = 16;
                    set(gcf,'Position',[50 50 1800 1000]);     

                    subplot(5,5,[1,2])
                    plot(1:sol.nb_cell_n*lenr,transpose(delta_coupled(1:sol.nb_cell_n*lenr)),'LineWidth',2);
                    title('delta coupled','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[4,5])
                    plot(1:sol.nb_cell_n*lenr,resized_csn_next,'LineWidth',2);
                    title('csn next','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[6,7])
                    plot(1:sol.nb_cell_n,fv.j(1:sol.nb_cell_n),'LineWidth',2);
                    title('j','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[9,10])
                    plot(1:sol.nb_cell_n,source_csn,'LineWidth',2);
                    title('source csn','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[11,12])
                    plot(1:sol.nb_cell_n,cse_next(1:sol.nb_cell_n),'LineWidth',2);
                    title('csen','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[14,15])
                    plot(sol.nb_cell_n+1:2*sol.nb_cell_n,cse_next(sol.nb_cell_n+1:2*sol.nb_cell_n),'LineWidth',2);
                    title('csep','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[16,17])
                    plot(1:sol.nb_cell_n*lenr,gcsn,'LineWidth',2);
                    title('gcsn','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[19,20])
                    plot(1:sol.nb_cell_n*lenr,Mn*(resized_csn_next-resized_csn)/sol.dt,'LineWidth',2);
                    title('M*dcsn/dt','fontsize',fs);
                    grid on
                    grid minor

                    subplot(5,5,[21,22])
                    plot(1:sol.nb_cell_n*lenr,- 0.5*(f_csn_next+f_csn),'LineWidth',2);
                    title('-fcsn','fontsize',fs);
                    grid on
                    grid minor
                end
                 


                % Bring variables back to normal levels

                min_c=0.00000;


                for iiii=1:1:length(ce_next)
                    if ce_next(iiii)<=0
                        ce_next(iiii)=min_c;
                    end
                end

                % Recompose cse (concentration at the surface of each particule.)

                for iiii =[1:sol.nb_cell_n]
                    cse_next(iiii)=resized_csn_next((sol.part_nb_cell+1)*iiii);
                end
                for iiii =[1:sol.nb_cell_p]
                    cse_next(iiii+sol.nb_cell_n)=resized_csp_next((sol.part_nb_cell+1)*(iiii));
                end

                if deb.videos_generation==1
                    if norm_delta_coupled>newt_lim
                        ce_next_save(:,newt_ite)= ce_next;
                        pe_next_save(:,newt_ite)= pe_next;
                        ps_next_save(:,newt_ite)= ps_next;
                    end 
                end

                % Determine if solver must be stopped and display convergence data

                if deb.prints>=2
                    
                    disp("newt_ite="+num2str(newt_ite)+" sol.time_ite="+num2str(sol.time_ite)+" residuals: coupled="+num2str(norm_delta_coupled)+" , ps="+num2str(norm_delta_ps)+" , pe="+num2str(norm_delta_pe) ...
                                                                             +" , ce="+num2str(norm_delta_ce)+" , csn="+num2str(norm_delta_csn)+" , csp="+num2str(norm_delta_csp))
                end

                %if (-(log10(max(hist.residuals(1,:)))-log10(hist.residuals(1,newt_ite)))<log10(newt_lim)) && newt_ite>1
                %if (-(log10(max(hist.residuals(1,:)))-log10(min(hist.residuals(1,1:newt_ite))))<log10(newt_lim)) && newt_ite>1
                %if (norm_delta_coupled<newt_lim || (-(log10(max(hist.residuals(1,:)))-log10(min(hist.residuals(1,1:newt_ite))))<log10(newt_lim)-1)) && newt_ite>1
                if (norm_delta_coupled<newt_lim) && newt_ite>1
                    
                    if deb.prints>=1
                        disp("Newton method for the coupled system converged:"+num2str(newt_ite)+" , "+num2str(norm_delta_coupled) ...
                                                                            +" , "+num2str(norm_delta_ps)+" , "+num2str(norm_delta_pe) ...
                                                                            +" , "+num2str(norm_delta_ce)+" , "+num2str(norm_delta_csn)+" , "+num2str(norm_delta_csp) ...
                                                                            +" , and "+num2str(-(log10(max(hist.residuals(1,:)))-log10(min(hist.residuals(1,1:newt_ite))))) ...
                                                                            +" , "+num2str(max(hist.residuals(1,:))) +" , "+num2str(min(hist.residuals(1,1:newt_ite))))
                    end
                    break
                elseif newt_ite==newt_max_ite
                    if deb.prints>=1
                        disp("Newton method for the coupled system did not converge:"+num2str(newt_ite)+" , "+num2str(norm_delta_coupled) ...
                                                                            +" , "+num2str(norm_delta_ps)+" , "+num2str(norm_delta_pe) ...
                                                                            +" , "+num2str(norm_delta_ce)+" , "+num2str(norm_delta_csn)+" , "+num2str(norm_delta_csp))
                        disp("Limit of the solver: "+num2str(log10(newt_lim))+" current order of magnitude of residual: "+ num2str(log10(norm_delta_coupled))+" det(Jac_coupled)="+num2str(det(Jac_coupled)) )
                    end
                    newt_ite=newt_ite+1;
                    break
                end
                

                deb.chrono_newtsol_update_singleite = deb.chrono_newtsol_update_singleite + toc(newtsol_update_timer);


            end

            if sol.time_ite==sol.nb_steps && deb.prints>=4  % Write the jacobian for debuging purposes
                writetable(array2table(Jac_coupled),'Jac_coupled_.xlsx')
            end

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
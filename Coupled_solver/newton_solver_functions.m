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

        function [ce_next,ce_next_save,pe_next,pe_next_save,ps_next,ps_next_save,csn_next,csp_next,norm_delta_coupled] = ...
                                            Newton_solver_coupled(obj,csn_next,csp_next,Dsn,Dsp,rn,...
                                                                rp,j,csn,csp,cse, ...
                                                                Mn,Mp,ce_next,pe_next,ps_next,...
                                                                De,kappa, sigma,dx,source_ce,...
                                                                source_csn,source_csp,source_pe,source_ps,current_source_psBC,...
                                                                ce,separator_index,newt_max_ite,M,dt,...
                                                                newt_lim,Jac_method,relax_factor,time_ite)
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
                [De,kappa,p.kappa_D_eff] = eq_build_fun.update_param_functions(ce,cse);
            end


            %calculate initial right hand side vectors (necessary for eqution with a time derivative (ie. for cs and ce))
            f_csn= eq_build_fun.LHS_f_cs(resized_csn,Dsn,rn,source_csn);
            f_csp= eq_build_fun.LHS_f_cs(resized_csp,Dsp,rp,source_csp);
            f_ce = eq_build_fun.LHS_f_ce(ce,De,dx,source_ce,"f_ce");
            f_ps = eq_build_fun.LHS_f_ps(fv.ps,ce,sigma,dx,source_ps,separator_index,current_source_psBC,"f_ps",0);
            f_pe = eq_build_fun.LHS_f_pe(fv.pe,ce,kappa,p.kappa_D_eff,dx,source_pe,"f_pe");

            ce_next_save= ones(length(ce_next),newt_max_ite);
            pe_next_save= ones(length(pe_next),newt_max_ite);
            ps_next_save= ones(length(ps_next),newt_max_ite);
            %% Start the newton method to solve g(c(t+dt))=0  (g(c(t+dt))=M(c(t+dt)-c(t))/dt -1/2 (f(t+dt)+f(t)))
            for newt_ite = 1:1:newt_max_ite

                %Update parameter functions and Butler-Volmer equation

                if sol.newton_update_sources==1
                    [De,kappa,p.kappa_D_eff] = eq_build_fun.update_param_functions(ce_next,cse_next);
                    fv.j=BV_fun.butler_volmer_equation(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,"newton_solver_functions next ite");

                   
                    [source_pe,source_ps,source_ce,source_csn,source_csp] = eq_build_fun.update_sources();
                    
                end

                %%Calculate the jacobian matrix of the left hand side vector for the concentration at time t+dt: Jf(c)
                %if method=1, jacobian is calculated analytically, if method=2, jac
                %is calculated using finite differences.
                if Jac_method==1
                    %Jac_f_ce_next = eq_build_fun.LHS_Jac_f_Fdiff_ce(ce_next,De,dx,source_ce);
                    %Jac_f_pe_next = eq_build_fun.LHS_Jac_f_Fdiff_pe(pe_next,ce_next,kappa,p.kappa_D_eff,dx,source_pe);
                    %Jac_f_ps_next = eq_build_fun.LHS_Jac_f_Fdiff_ps(ps_next,ce_next,sigma,dx,source_ps,separator_index,current_source_psBC);
                    %Jac_f_csn_next= eq_build_fun.LHS_Jac_f_Fdiff_cs(resized_csn_next,Dsn,rn,source_csn);
                    %Jac_f_csp_next= eq_build_fun.LHS_Jac_f_Fdiff_cs(resized_csp_next,Dsp,rp,source_csp);

                    Jac_f_ce_next = eq_build_fun.LHS_Jac_f_Fdiff_ce(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,De,dx,source_ce);
                    Jac_f_pe_next = eq_build_fun.LHS_Jac_f_Fdiff_pe(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,kappa,p.kappa_D_eff,dx,source_pe);
                    Jac_f_ps_next = eq_build_fun.LHS_Jac_f_Fdiff_ps(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,sigma,dx,source_ps,separator_index,current_source_psBC);
                    Jac_f_csn_next= eq_build_fun.LHS_Jac_f_Fdiff_cs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csn_next,Dsn,rn,source_csn);
                    Jac_f_csp_next= eq_build_fun.LHS_Jac_f_Fdiff_cs(pe_next,ps_next,ce_next,cse_next,ini.T0,fv.Ueq,resized_csp_next,Dsp,rp,source_csp);



                    if (1==1 |det(Jac_f_pe_next)==0)
                        if (deb.prints>0)
                            disp("The jacobian for the potential in the electrolyte has a determinant of 0. The diagonal boost method is used to make the matrix invertible.")
                        end
                        for iiii=1:1:len
                            Jac_f_pe_next(iiii,iiii)=Jac_f_pe_next(iiii,iiii)*1.0000000001;
                        end
                    end
                    if (det(Jac_f_ps_next(1:sol.nb_cell_n,1:sol.nb_cell_n))==0)
                        if (deb.prints>0)
                            disp("The jacobian for the potential in the neg. electrode has a determinant of 0. The diagonal boost method is used to make the matrix invertible.")
                        end
                        for iiii=1:1:sol.nb_cell_n
                            Jac_f_ps_next(iiii,iiii)=Jac_f_ps_next(iiii,iiii)*1.0000000001;
                        end
                    end
                    if 1==1 | (det(Jac_f_ps_next(sol.nb_cell_n+1:sol.nb_cell_p+sol.nb_cell_n,sol.nb_cell_n+1:sol.nb_cell_p+sol.nb_cell_n))==0)
                        if (deb.prints>0)
                            disp("The jacobian for the potential in the pos. electrode has a determinant of 0. The diagonal boost method is used to make the matrix invertible.")
                        end
                        for iiii=sol.nb_cell_n+1:1:sol.nb_cell_p+sol.nb_cell_n
                            Jac_f_ps_next(iiii,iiii)=Jac_f_ps_next(iiii,iiii)*1.0000000001;
                        end
                    end
                end


                f_ce_next       = eq_build_fun.LHS_f_ce(ce_next,De,dx,source_ce,"f_ce_next");
                f_pe_next       = eq_build_fun.LHS_f_pe(pe_next,ce_next,kappa,p.kappa_D_eff,dx,source_pe,"f_pe_next");
                f_ps_next       = eq_build_fun.LHS_f_ps(ps_next,ce_next,sigma,dx,source_ps,separator_index,current_source_psBC,"f_ps_next",0);
                f_csn_next      = eq_build_fun.LHS_f_cs(resized_csn_next,Dsn,rn,source_csn);
                f_csp_next      = eq_build_fun.LHS_f_cs(resized_csp_next,Dsp,rp,source_csp);
                
                %%Compilation of the Crank Nicholson method s function and jacobian
                %g(c(t+dt))=M(c(t+dt)-c(t))/dt -1/2 (f(t+dt)+f(t))
                %Jac_g(c(t+dt))=M/dt -1/2 jac_f(t+dt)
                %Jac_g*(delta(c(t+dt)))=g
                

                %Solve Li concentration in the solid
                gcsn=Mn*(resized_csn_next-resized_csn)/dt - 0.5*(f_csn_next+f_csn);
                gcsp=Mp*(resized_csp_next-resized_csp)/dt - 0.5*(f_csp_next+f_csp);
                Jac_gcsn=Mn/dt -0.5*Jac_f_csn_next;
                Jac_gcsp=Mp/dt -0.5*Jac_f_csp_next;

                %Solve Li concentration in electrolyte
                gc=M*(ce_next-ce)/dt - 0.5*(f_ce_next+f_ce);
                Jac_gc=M/dt -0.5*Jac_f_ce_next;


                %Solve electric potential in electrolyte
                gp = f_pe_next;
                Jac_gp= Jac_f_pe_next;

                %Solve electric potential in solid
                gps = f_ps_next;
                Jac_gps= Jac_f_ps_next;

                % Form the coupled system 
                Jac_coupled = blkdiag(Jac_gcsn,Jac_gcsp,Jac_gc,Jac_gp,Jac_gps) ;

                
                g_coupled   = cat(1,gcsn,gcsp);
                g_coupled   = cat(1,g_coupled,gc);
                g_coupled   = cat(1,g_coupled,gp);
                g_coupled   = cat(1,g_coupled,gps);

                % Solve the system
                    
                %delta_coupled = (Jac_coupled) \ (-g_coupled);
                delta_coupled = linsolve(Jac_coupled,-g_coupled);
                
                
                delta_csn = linsolve(Jac_gcsn,-gcsn);
                delta_csp = linsolve(Jac_gcsp,-gcsp);
                delta_c = linsolve(Jac_gc,-gc);
                delta_p = linsolve(Jac_gp,-gp);
                delta_ps = linsolve(Jac_gps,-gps);
                delta_ps1 = linsolve(Jac_gps(1:3,1:3),-gps(1:3));
                delta_ps2 = linsolve(Jac_gps(4:6,4:6),0*gps(4:6));

                %delta_coupled = cat(1,delta_csn,delta_csp);
                %delta_coupled = cat(1,delta_coupled,delta_c);
                %delta_coupled = cat(1,delta_coupled,delta_p);
                %delta_coupled = cat(1,delta_coupled,delta_ps);

                norm_delta_coupled = sqrt(sum(transpose(delta_coupled) * delta_coupled))/length(g_coupled);
                norm_delta_ps = sqrt(sum(transpose(delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len+len_ps)) * delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len+len_ps)))/len_ps;
                norm_delta_pe = sqrt(sum(transpose(delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len)) * delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len)))/len;
                norm_delta_ce = sqrt(sum(transpose(delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len)) * delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len)))/len;
                norm_delta_csp = sqrt(sum(transpose(delta_coupled(1 : sol.nb_cell_n*lenr)) * delta_coupled(1 : sol.nb_cell_n*lenr)))/(sol.nb_cell_n*lenr);
                norm_delta_csn = sqrt(sum(transpose(delta_coupled(sol.nb_cell_n*lenr+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr)) * delta_coupled(sol.nb_cell_n*lenr+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr)))/(sol.nb_cell_p*lenr);

                hist.residuals(1,newt_ite)=norm_delta_coupled;
                hist.residuals(2,newt_ite)=norm_delta_ps;
                hist.residuals(3,newt_ite)=norm_delta_pe;
                hist.residuals(4,newt_ite)=norm_delta_ce;
                hist.residuals(5,newt_ite)=norm_delta_csn;
                hist.residuals(6,newt_ite)=norm_delta_csp;

                if deb.prints>=2
                    

                    disp("DEBUG BEN j ------------------------------ newt_ite="+num2str(newt_ite)+" , time ite="+num2str(time_ite))
                    disp((fv.j))
                    pe_next_reduced = cat(1,pe_next(1:sol.nb_cell_n),pe_next(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell));
                    disp(transpose(ps_next-pe_next_reduced))


                    disp("DEBUG BEN cs newt_ite="+num2str(newt_ite)+" , time ite="+num2str(time_ite))

                    disp(transpose(gcsn))
                    disp(transpose(resized_csn_next))
                    disp(transpose(resized_csp_next))
                    %disp((Jac_gcsn))
                    disp(transpose(delta_coupled(1 : sol.nb_cell_n*lenr)))
                    disp(transpose(delta_coupled(sol.nb_cell_n*lenr+1:sol.nb_cell_n*lenr+sol.nb_cell_p*lenr)))
                    %disp(source_csn)
                    %disp(source_csp)
                    %disp(p.A_s_n)

                    disp("DEBUG BEN ce newt_ite="+num2str(newt_ite)+" , time ite="+num2str(time_ite))
                    disp(transpose(gc))
                    disp(transpose(ce_next))
                    disp((Jac_gc))
                    disp(transpose(delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len)))


                    disp("DEBUG BEN pe newt_ite="+num2str(newt_ite)+" , time ite="+num2str(time_ite))

                    disp(transpose(gp))
                    disp(transpose(pe_next))
                    disp((Jac_gp))
                    disp(transpose(delta_coupled(sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+1 : sol.nb_cell_n*lenr+sol.nb_cell_p*lenr+len+len)))
                    
                    disp("DEBUG BEN ps newt_ite="+num2str(newt_ite)+" , time ite="+num2str(time_ite))
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
                    disp("Newton method for lithium concentration in electrlyte did not converge: "+num2str(newt_ite)+" , "+num2str(norm_delta_coupled))
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

                min_c=0.001;

                for iiii=1:1:length(resized_csn_next)
                    if resized_csn_next(iiii)<min_c
                        resized_csn_next(iiii)=min_c;
                    elseif  resized_csn_next(iiii)>p.csn_max
                        resized_csn_next(iiii)=p.csn_max;
                    end
                end

                for iiii=1:1:length(resized_csp_next)
                    if resized_csp_next(iiii)<min_c
                        resized_csp_next(iiii)=min_c;
                    elseif  resized_csp_next(iiii)>p.csp_max
                        resized_csp_next(iiii)=p.csp_max;
                    end
                end

                for iiii=1:1:length(ce_next)
                    if ce_next(iiii)<=0
                        ce_next(iiii)=0.000001;
                    end
                end

                for iiii=1:1:length(ps_next)
                    cap=100000;
                    if abs(ps_next(iiii))>cap
                        ps_next(iiii)=cap*ps_next(iiii)/abs(ps_next(iiii));
                    end
                end
                
                for iiii =[1:sol.nb_cell_n+sol.nb_cell_p]
                    if iiii<=sol.nb_cell_n
                        cse_next(iiii)=resized_csn_next(sol.part_nb_cell*iiii);
                    else
                        cse_next(iiii)=resized_csn_next(sol.part_nb_cell*(iiii-sol.nb_cell_n));
                    end
                end

                if deb.videos_generation==1
                    if norm_delta_coupled>newt_lim
                        ce_next_save(:,newt_ite)= ce_next;
                        pe_next_save(:,newt_ite)= pe_next;
                        ps_next_save(:,newt_ite)= ps_next;
                    end 
                end

                if deb.prints>=1
                    
                    disp("the ite number and residual are:"+num2str(newt_ite)+" , "+num2str(norm_delta_coupled)+" , "+num2str(norm_delta_ps)+" , "+num2str(norm_delta_pe) ...
                                                                             +" , "+num2str(norm_delta_ce)+" , "+num2str(norm_delta_csn)+" , "+num2str(norm_delta_csp))
                end


                %if (-(log10(max(hist.residuals(1,:)))-log10(hist.residuals(1,newt_ite)))<log10(newt_lim)) && newt_ite>1
                %if (-(log10(max(hist.residuals(1,:)))-log10(min(hist.residuals(1,1:newt_ite))))<log10(newt_lim)) && newt_ite>1
                %if (norm_delta_coupled<newt_lim || (-(log10(max(hist.residuals(1,:)))-log10(min(hist.residuals(1,1:newt_ite))))<log10(newt_lim)-1)) && newt_ite>1
                if (norm_delta_coupled<newt_lim) && newt_ite>1
                    
                    if deb.prints>=0
                        disp("Newton method for lithium concentration in electrlyte has converged:"+num2str(newt_ite)+" , "+num2str(norm_delta_coupled) ...
                                                                            +" , "+num2str(norm_delta_ps)+" , "+num2str(norm_delta_pe) ...
                                                                            +" , "+num2str(norm_delta_ce)+" , "+num2str(norm_delta_csn)+" , "+num2str(norm_delta_csp) ...
                                                                            +" , "+num2str(-(log10(max(hist.residuals(1,:)))-log10(min(hist.residuals(1,1:newt_ite))))) ...
                                                                            +" , "+num2str(max(hist.residuals(1,:))) +" , "+num2str(min(hist.residuals(1,1:newt_ite))))
                    end
                    break
                elseif newt_ite==newt_max_ite
                    
                    disp("Newton method for lithium concentration in electrlyte did not converge:"+num2str(newt_ite)+" , "+num2str(norm_delta_coupled) ...
                                                                            +" , "+num2str(norm_delta_ps)+" , "+num2str(norm_delta_pe) ...
                                                                            +" , "+num2str(norm_delta_ce)+" , "+num2str(norm_delta_csn)+" , "+num2str(norm_delta_csp))
                    disp("Limit of the solver: "+num2str(log10(newt_lim))+" current order of magnitude of residual: "+ num2str(-(log10(max(hist.residuals(1,:)))-log10(min(hist.residuals(1,:))))) ...
                                                                                                                +"  "+ num2str(log10(norm_delta_coupled)) )
                    break
                end
                
            end

            if deb.prints>=0
                hist.residuals_time(1,time_ite)=norm_delta_coupled;
                hist.residuals_time(2,time_ite)=norm_delta_ps;
                hist.residuals_time(3,time_ite)=norm_delta_pe;
                hist.residuals_time(4,time_ite)=norm_delta_ce;
                hist.residuals_time(5,time_ite)=norm_delta_csn;
                hist.residuals_time(6,time_ite)=norm_delta_csp;
                hist.newt_it_number(1,time_ite)=newt_ite;

                if isnan(norm_delta_coupled)==0
                    hist.residuals_diff(1,time_ite)=-(log10(max(hist.residuals(1,:)))-log10(hist.residuals(1,newt_ite)));
                    hist.residuals_diff(2,time_ite)=-(log10(max(hist.residuals(2,:)))-log10(hist.residuals(2,newt_ite)));
                    hist.residuals_diff(3,time_ite)=-(log10(max(hist.residuals(3,:)))-log10(hist.residuals(3,newt_ite)));
                    hist.residuals_diff(4,time_ite)=-(log10(max(hist.residuals(4,:)))-log10(hist.residuals(4,newt_ite)));
                    hist.residuals_diff(5,time_ite)=-(log10(max(hist.residuals(5,:)))-log10(hist.residuals(5,newt_ite)));
                    hist.residuals_diff(6,time_ite)=-(log10(max(hist.residuals(6,:)))-log10(hist.residuals(6,newt_ite)));
                end
            end

            if (deb.prints>0)
                %disp(transpose(resized_csn_next))
                %disp(transpose(resized_csp_next))
                %disp(csn_next)
                %disp(csp_next)
                %disp(test_array)
                %disp(resized_test)
                %disp(test_array_2)
                %disp(csn_next)
                %disp(csn_next(1,1)-100)
                %disp(csp_next)
                %test_array_2 = reshape(resized_test,sol.part_nb_cell+1,sol.nb_cell_n);
            end

            csn_next  = reshape(resized_csn_next,sol.part_nb_cell+1,sol.nb_cell_n);
            csp_next  = reshape(resized_csp_next,sol.part_nb_cell+1,sol.nb_cell_n);
            
        end
    end
end
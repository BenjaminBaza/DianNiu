% Following the method of Xia et al. 2017 and Zheng et al. 2013

function DFN_equations_coupled_solver()

    global sol_fun
    global BV_fun
    global ini
    global deb
    global sol
    global ex
    global p
    global fv
    global eq_build_fun


    %%create left hand side matrices M

    if exist('M_electrlyte_diff','var') == 0
        epsilon_e = cat(1,p.eps_e_n*sol.dxn*ones(sol.nb_cell_n,1), p.eps_e_s*sol.dxs*ones(sol.nb_cell_s,1));
        epsilon_e = cat(1, epsilon_e , p.eps_e_p*sol.dxp*ones(sol.nb_cell_p,1));
        M_electrlyte_diff=sol_fun.calculate_mass_matrix_electrolyte(epsilon_e);
        clear epsilon_e
    end

    if exist('Mp_part_diff','var') == 0
        max_ind=sol.part_nb_cell+1;
        drn=ones(1,max_ind);
        drn(1)=(sol.part_coord_n(2)-sol.part_coord_n(1));
        drn(max_ind)=(sol.part_coord_n(max_ind)-sol.part_coord_n(max_ind-1));
        drp=ones(1,max_ind);
        drp(1)=(sol.part_coord_p(2)-sol.part_coord_p(1));
        drp(max_ind)=(sol.part_coord_p(max_ind)-sol.part_coord_p(max_ind-1));
        for ind=2:1:max_ind-1
            drn(ind)=(sol.part_coord_n(ind+1)-sol.part_coord_n(ind-1))/2;
            drp(ind)=(sol.part_coord_p(ind+1)-sol.part_coord_p(ind-1))/2;
        end

        Volumes_n=ones(1,max_ind);
        Volumes_p=ones(1,max_ind);
        for ind=1:1:max_ind
            Volumes_n(ind)=sol.part_coord_n(ind)^2*drn(ind) + drn(ind)^3/12.0;
            Volumes_p(ind)=sol.part_coord_p(ind)^2*drp(ind) + drp(ind)^3/12.0;
        end

        Mp_part_diff=sol_fun.calculate_mass_matrix(max_ind,Volumes_p);
        Mn_part_diff=sol_fun.calculate_mass_matrix(max_ind,Volumes_n);
        Mn_dummy=Mn_part_diff;
        Mp_dummy=Mp_part_diff;

        for iiii = 1:1:sol.nb_cell_n-1
            Mn_part_diff=blkdiag(Mn_part_diff,Mn_dummy);
        end
        for iiii = 1:1:sol.nb_cell_p-1
            Mp_part_diff=blkdiag(Mp_part_diff,Mp_dummy);
        end
    end

    solver_converged=0;
    solver_attempts =0;
    while solver_converged==0
        solver_attempts = solver_attempts + 1;

        csn_reshaped=reshape(fv.csn,sol.nb_cell_n*(sol.part_nb_cell+1),1);
        csp_reshaped=reshape(fv.csp,sol.nb_cell_p*(sol.part_nb_cell+1),1);

        [p.De_eff,p.kappa,p.kappa_D_eff] = eq_build_fun.update_param_functions(fv.ce,fv.cse,csn_reshaped,csp_reshaped);

        %% Build sources from values of the last iteration
        current_source_n=-ex.I_array(sol.time_ite)/(p.coll_A_n*p.sig_eff_n);
        current_source_p=-ex.I_array(sol.time_ite)/(p.coll_A_p*p.sig_eff_p);
        current_source_psBC=[current_source_n,current_source_p];

        fv.j=BV_fun.butler_volmer_equation(fv.pe,fv.ps,fv.ce,fv.cse,ini.T0,fv.Ueq,"DFN_equations_coupled_solver original");

        [source_pe,source_ps,source_ce,source_csn,source_csp] = eq_build_fun.update_sources();

        if deb.prints>1
            disp("DEBUG BEN DFN_equations_coupled_solver j and sources for ps, pe and ce DFNeq "+num2str(sol.time_ite)+" "+num2str(solver_attempts))
            disp(transpose(fv.Ueq))
            disp(fv.csn)
            disp(transpose(fv.ps))
            disp(fv.cse)
            

        end


        %% Initialize the vector for Li concentration and electric potential in the electrolyte at the next time step: (c(t+dt)=c(t))
        ce_next=fv.ce;
        pe_next=fv.pe;
        ps_next=fv.ps;
        csn_next=fv.csn;
        csp_next=fv.csp;
        %%Calculate the left hand side vector for the concentration at time t: f(c) and solve the system

        [ce_next,ce_next_save,pe_next,pe_next_save,ps_next,ps_next_save,csn_next,csp_next,rms,sol.newt_ite] = ...
                                                sol_fun.Newton_solver_coupled(  csn_next,csp_next,sol.part_coord_n,...
                                                                                sol.part_coord_p,fv.j, fv.csn,fv.csp,fv.cse,...
                                                                                Mn_part_diff,Mp_part_diff,ce_next,pe_next,ps_next,...
                                                                                p.kappa_eff,p.sig_eff,sol.cell_dx,source_ce, ...
                                                                                source_csn, source_csp,source_pe,source_ps,current_source_psBC,...
                                                                                fv.ce,sol.nb_cell_n,sol.newton_meth_max_ite,M_electrlyte_diff, ...
                                                                                sol.newton_meth_res_threshold,1,sol.newton_relax_factor);


        sol.time_ite_save=sol.time_ite;
        if isnan(rms)~=0 | isreal(rms) == 0 | sol.newt_ite ==sol.newton_meth_max_ite+1
            disp("Newton solver at time ite "+num2str(sol.time_ite)+" failed to converge at attempt "+num2str(solver_attempts))
            if solver_attempts==30
                deb.break_time_loop=1;
                sol.time_ite=sol.time_ite-1;
                break
            else
                %sol.dt=sol.max_dt/(2^solver_attempts);
                sol.dt=sol.dt/(2^solver_attempts);
            end
        else
            solver_converged=1;
        end

    end

    fv.ce=ce_next;
    fv.pe=pe_next;
    fv.ps=ps_next;
    fv.csn=csn_next;
    fv.csp=csp_next;
    fv.cse=cat(2,fv.csn(length(sol.part_coord_n),:),fv.csp(length(sol.part_coord_p),:));
    


end
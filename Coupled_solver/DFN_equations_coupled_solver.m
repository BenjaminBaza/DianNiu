% Following the method of Xia et al. 2017 and Zheng et al. 2013

global sol_fun
global BV_fun
global ini
global deb

Current_intensity_temp=0.0001;
current_source_n=-Current_intensity_temp/(p.coll_A_n*p.sig_eff_n);
current_source_p=-Current_intensity_temp/(p.coll_A_p*p.sig_eff_p);
current_source_psBC=[current_source_n,current_source_p];



fv.j=BV_fun.butler_volmer_equation(fv.pe,fv.ps,fv.ce,fv.cse,p.k0,p.alpha,p.Faraday,p.Rg, ...
                    ini.T0,fv.Ueq,p.Rfilm,p.csn_max,p.csp_max,sol.nb_cell_n,sol.nb_cell_s);


source_n=p.A_s_n*fv.j(1:sol.nb_cell_n);
source_s=0*ones(1,sol.nb_cell_s);
source_p=p.A_s_p*fv.j(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell);

source_ce=cat(2,(1-p.t_plus)*source_n,source_s);
source_ce=cat(2,source_ce,(1-p.t_plus)*source_p);

source_pe=cat(2,p.Faraday*source_n,source_s);
source_pe=cat(2,source_pe,p.Faraday*source_p);

source_ps=-p.Faraday*source_n;
source_ps=cat(2,source_ps,-p.Faraday*source_p);

source_csn=-fv.j(1:sol.nb_cell_n)/p.Dsn;
source_csp=-fv.j(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell)/p.Dsp;

disp("sources!!!!1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 ")
disp(source_pe)
disp(source_ps)
disp(source_ce)
disp(source_csn)
disp(source_csp)

[source_pe,source_ps,source_ce,source_csn,source_csp] = eq_build_fun.update_sources();
disp("sources!!!!2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222 ")
disp(source_pe)
disp(source_ps)
disp(source_ce)
disp(source_csn)
disp(source_csp)
disp("sources!!!!333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333 ")



if deb.prints>1
    disp("DEBUG BEN j and sources for ps, pe and ce")
    disp(fv.j)
    disp(source_pe)
    disp(source_ps)
    disp(transpose(p.kappa_eff))
    disp(transpose(p.sig_eff))
end

%%create left hand side matrices M

if exist('M_electrlyte_diff','var') == 0
    epsilon_e = cat(1,p.eps_e_n*sol.dxn*ones(sol.nb_cell_n,1), p.eps_e_s*sol.dxs*ones(sol.nb_cell_s,1));
    epsilon_e = cat(1, epsilon_e , p.eps_e_p*sol.dxp*ones(sol.nb_cell_p,1));
    M_electrlyte_diff=sol_fun.calculate_mass_matrix_electrolyte(epsilon_e);
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

%% Initialize the vector for Li concentration and electric potential in the electrolyte at the next time step: (c(t+dt)=c(t))
ce_next=fv.ce;
pe_next=fv.pe;
ps_next=fv.ps;
csn_next=fv.csn;
csp_next=fv.csp;
%%Calculate the left hand side vector for the concentration at time t: f(c) and solve the system

segregation=0;
[ce_next,ce_next_save,pe_next,pe_next_save,ps_next,ps_next_save,csn_next,csp_next,rms] = ...
                                        sol_fun.Newton_solver_coupled(  csn_next,csp_next,p.Dsn,p.Dsp,sol.part_coord_n,...
                                                                        sol.part_coord_p,fv.j, fv.csn,fv.csp,fv.cse,...
                                                                        Mn_part_diff,Mp_part_diff,ce_next,pe_next,ps_next,...
                                                                        p.De_eff,p.kappa_eff,p.sig_eff,sol.cell_dx,source_ce, ...
                                                                        source_csn, source_csp,source_pe,source_ps,current_source_psBC,...
                                                                        fv.ce,sol.nb_cell_n,sol.newton_meth_max_ite,M_electrlyte_diff,sol.dt, ...
                                                                        sol.newton_meth_res_threshold,1,sol.newton_relax_factor,ite);

if isnan(rms)~=0 | isreal(rms) == 0
    deb.break_time_loop=1;
end

fv.ce=ce_next;
fv.pe=pe_next;
fv.ps=ps_next;
fv.csn=csn_next;
fv.csp=csp_next;
fv.cse=cat(2,fv.csn(length(sol.part_coord_n),:),fv.csp(length(sol.part_coord_p),:));


if deb.videos_generation==1
    global vis_fun
    vis_fun.animate_data('ce.avi',sol.cell_center_coord,ce_next_save,1,'x','c_e','Li concentration in electrolyte');
    vis_fun.animate_data('pe.avi',sol.cell_center_coord,pe_next_save,2,'x','phi_e','Potential in electrolyte');
    vis_fun.animate_data_solid('ps.avi',sol.cell_center_coord,ps_next_save(1:sol.nb_cell_n,:), ...
                                ps_next_save(sol.nb_cell_n+1:sol.nb_cell_n+sol.nb_cell_p,:),3,'x','phi_s','Potential in the solid' ...
                                , sol.nb_cell_n, sol.nb_cell_s, sol.nb_cell_p);
end




%disp(transpose(ce_next));

clear ce_next j_temp temporary_j;
clear max_ind_sep max_ind_ano max_ind_cat min_ind_sep min_ind_ano min_ind_cat ;
clear source_n source_s source_p source_ce source_ps source_pe temporary_j relax_factor;
clear pe_next ps_next ce_next ; %ce_next_save pe_next_save ps_next_save;
clear drn drp Volumes_n Volumes_p max_ind drn drp;
clear Mn_dummy Mp_dummy iiii ind;
% Following the method of Xia et al. 2017

global BV_fun

if deb.prints==1
    disp('Solving diffusion in electrolyte');
end 

Current_intensity_temp=1;
current_source_n=-Current_intensity_temp/(p.coll_A_n*p.sig_eff_n);
current_source_p=-Current_intensity_temp/(p.coll_A_p*p.sig_eff_p);
current_source_psBC=[current_source_n,current_source_p];



relax_factor=1.;

temporary_j=0;
if temporary_j==1
    jn=5;%-5.35e-5; %DEBUG BEN temporary hard coded value for development purposes
    jp=5;%-5.35e-5; %DEBUG BEN temporary hard coded value for development purposes
    source_n=p.A_s_n*jn*ones(1,sol.nb_cell_n);
    source_s=0*ones(1,sol.nb_cell_s);
    source_p=p.A_s_p*jp*ones(1,sol.nb_cell_p);
else

    fv.j=butler_volmer_eq(fv.pe,fv.ps,fv.ce,fv.cse,p.k0,p.alpha,p.Faraday,p.Rg, ...
                        ini.T0,fv.Ueq,p.Rfilm,p.csn_max,p.csp_max,sol.nb_cell_n,sol.nb_cell_s);
    source_n=p.A_s_n*fv.j(1:sol.nb_cell_n);
    source_s=0*ones(1,sol.nb_cell_s);
    source_p=p.A_s_p*fv.j(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell);
end

source_ce=cat(2,(1-p.t_plus)*source_n,source_s);
source_ce=cat(2,source_ce,(1-p.t_plus)*source_p);

source_pe=cat(2,p.Faraday*source_n,source_s);
source_pe=cat(2,source_pe,p.Faraday*source_p);

source_ps=-p.Faraday*source_n;
source_ps=cat(2,source_ps,-p.Faraday*source_p);

if deb.prints==1
    disp("DEBUG BEN j and sources for ps, pe and ce")
    disp(fv.j)
    disp(source_pe)
    disp(p.kappa_eff)
    disp(p.sig_eff)
end

%%create left hand side matrices M
if exist('epsilon_e','var') == 0
    epsilon_e = cat(1,p.eps_e_n*sol.dxn*ones(sol.nb_cell_n,1), p.eps_e_s*sol.dxs*ones(sol.nb_cell_s,1));
    epsilon_e = cat(1, epsilon_e , p.eps_e_p*sol.dxp*ones(sol.nb_cell_p,1));
end
M_electrlyte_diff=calculate_mass_matrix_electrolyte(epsilon_e);

%% Initialize the vector for Li concentration and electric potential in the electrolyte at the next time step: (c(t+dt)=c(t))
ce_next=fv.ce;
pe_next=fv.pe;
ps_next=fv.ps;
%%Calculate the left hand side vector for the concentration at time t: f(c) and solve the system

segregation=1;
[ce_next,ce_next_save,pe_next,pe_next_save,ps_next,ps_next_save] = ...
                                            Newton_concentration_electrolyte(ce_next,pe_next,ps_next,p.De_eff,p.kappa_eff, ...
                                            p.sig_eff,sol.cell_dx,source_ce,source_pe,source_ps,current_source_psBC,fv.ce, ...
                                            sol.nb_cell_n,sol.newton_meth_max_ite,M_electrlyte_diff,sol.dt, ...
                                            sol.newton_meth_res_threshold,1,relax_factor, segregation);


fv.ce=ce_next;
fv.pe=pe_next;
fv.ps=ps_next;


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
clear pe_next ps_next ce_next ; %ce_next_save pe_next_save ps_next_save


%% Functions

function M = calculate_mass_matrix_electrolyte(eps)
    M1=diag(eps);
    M=M1;
end

function [ce_next,ce_next_save,pe_next,pe_next_save,ps_next,ps_next_save] = ...
                                            Newton_concentration_electrolyte(ce_next,pe_next,ps_next,D,kappa, ...
                                            sigma,dx,source_ce,source_pe,source_ps,current_source_psBC,ce, ...
                                            separator_index,newt_max_ite,M,dt,newt_lim,Jac_method,relax_factor,segregation_switch)
    global deb
    len=length(ce_next);
    len_ps=length(ps_next);

    f_ce=LHS_f_ce(ce,D,dx,source_ce);
    ce_next_save= ones(length(ce_next),newt_max_ite);
    pe_next_save= ones(length(pe_next),newt_max_ite);
    ps_next_save= ones(length(ps_next),newt_max_ite);
    %% Start the newton method to solve g(c(t+dt))=0  (g(c(t+dt))=M(c(t+dt)-c(t))/dt -1/2 (f(t+dt)+f(t)))
    for newt_ite = 1:1:newt_max_ite

        %%Calculate the jacobian matrix of the left hand side vector for the concentration at time t+dt: Jf(c)
        %if method=1, jacobian is calculated analytically, if method=2, jac
        %is calculated using finite differences.
        if Jac_method==1
            Jac_f_ce_next=LHS_Jac_f_Fdiff_ce(ce_next,D,dx,source_ce);
            Jac_f_pe_next=LHS_Jac_f_Fdiff_pe(pe_next,ce_next,kappa,dx,source_pe);
            Jac_f_ps_next=LHS_Jac_f_Fdiff_ps(ps_next,ce_next,sigma,dx,source_ps,separator_index,current_source_psBC);
        end
        f_ce_next=LHS_f_ce(ce_next,D,dx,source_ce);
        f_pe_next=LHS_f_pe(pe_next,ce_next,kappa,dx,source_pe);
        f_ps_next=LHS_f_ps(ps_next,ce_next,sigma,dx,source_ps,separator_index,current_source_psBC);

        
        %disp("in Newton method")
        %disp(transpose(f_c_next))
        %disp(transpose(Jac_f_c_next))
        
        %%Compilation of the Crank Nicholson method s function and jacobian
        %g(c(t+dt))=M(c(t+dt)-c(t))/dt -1/2 (f(t+dt)+f(t))
        %Jac_g(c(t+dt))=M/dt -1/2 jac_f(t+dt)
        %Jac_g*(delta(c(t+dt)))=g
        
        if segregation_switch==1
            %Solve Li concentration in electrolyte
            gc=M*(ce_next-ce)/dt - 0.5*(f_ce_next+f_ce);
            Jac_gc=M/dt -0.5*Jac_f_ce_next;

            delta_c = linsolve(Jac_gc,-gc);
            norm_delta_c = sum(transpose(delta_c) * delta_c)/len;
            ce_next=ce_next + relax_factor*delta_c;

            %Solve electric potential in electrolyte
            gp = f_pe_next;
            Jac_gp= Jac_f_pe_next;

            delta_p = linsolve(Jac_gp,-gp);
            norm_delta_p = sum(transpose(delta_p) * delta_p)/len;
            pe_next=pe_next + relax_factor*delta_p;

            %Solve electric potential in solid
            gps = f_ps_next;
            Jac_gps= Jac_f_ps_next;

            delta_ps = linsolve(Jac_gps,-gps);
            norm_delta_ps = sum(transpose(delta_ps) * delta_ps)/len_ps;
            ps_next=ps_next + relax_factor*delta_ps;


            if deb.prints==1
                disp("the ite number and residual are:"+num2str(newt_ite)+" , "+num2str(norm_delta_c)+" , "+num2str(norm_delta_p)+" and "+num2str(norm_delta_ps))
            end

            if max(norm_delta_c,max(norm_delta_p,norm_delta_ps))<newt_lim
                if deb.prints==1
                    disp("Newton method for lithium concentration in electrlyte has converged")
                end
                break
            elseif newt_ite==newt_max_ite
                disp("Newton method for lithium concentration in electrlyte did not converge")
            end
            if deb.videos_generation==1
                if norm_delta_c>newt_lim
                    ce_next_save(:,newt_ite)= ce_next;
                end
                if norm_delta_p>newt_lim
                    pe_next_save(:,newt_ite)= pe_next;
                end
                if norm_delta_ps>newt_lim
                    ps_next_save(:,newt_ite)= ps_next;
                end
            end
        else
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
            Jac_coupled = blkdiag(Jac_gc,Jac_gp,Jac_gps) ;
            g_coupled   = cat(1,gc,gp);
            g_coupled   = cat(1,g_coupled,gps);

            % Solve the system
            delta_coupled = linsolve(Jac_coupled,-g_coupled);
            norm_delta_coupled = sum(transpose(delta_coupled) * delta_coupled)/length(g_coupled);

            % Update solutions 
            ce_next=ce_next + relax_factor*delta_coupled(1:len);
            pe_next=pe_next + relax_factor*delta_coupled(len+1:len+len);
            ps_next=ps_next + relax_factor*delta_coupled(2*len+1:2*len+len_ps);


            if deb.prints==1
                disp("the ite number and residual are:"+num2str(newt_ite)+" , "+num2str(norm_delta_coupled))
            end

            if norm_delta_coupled<newt_lim
                if deb.prints==1
                    disp("Newton method for lithium concentration in electrlyte has converged:"+num2str(newt_ite)+" , "+num2str(norm_delta_coupled))
                end
                break
            elseif newt_ite==newt_max_ite
                disp("Newton method for lithium concentration in electrlyte did not converge:"+num2str(newt_ite)+" , "+num2str(norm_delta_coupled))
            end
            if deb.videos_generation==1
                if norm_delta_coupled>newt_lim
                    ce_next_save(:,newt_ite)= ce_next;
                    pe_next_save(:,newt_ite)= pe_next;
                    ps_next_save(:,newt_ite)= ps_next;
                end
                
            end
        end
        
    end
end

function [dxm,dxc,dxp,pm,pc,pp,cm,cc,cp,Dm,Dc,Dp]= concentration_electrolyte_BC(p,c,dx,i,D)
    len=length(c);
    if i==len
        dxm=dx(i-1);dxc=dx(i);dxp=dx(i);
        cm=c(i-1);cc=c(i);cp=c(i);
        pm=p(i-1);pc=p(i);pp=p(i);
        Dm=D(i-1);Dc=D(i);Dp=D(i);
    elseif i==1
        dxm=dx(i);dxc=dx(i);dxp=dx(i+1);
        cm=c(i);cc=c(i);cp=c(i+1);
        pm=p(i);pc=p(i);pp=p(i+1);
        Dm=D(i);Dc=D(i);Dp=D(i+1);
    else
        dxm=dx(i-1);dxc=dx(i);dxp=dx(i+1);
        cm=c(i-1);cc=c(i);cp=c(i+1);
        pm=p(i-1);pc=p(i);pp=p(i+1);
        Dm=D(i-1);Dc=D(i);Dp=D(i+1);
    end
end

function [dxm,dxc,dxp,pm,pc,pp,cm,cc,cp,Dm,Dc,Dp]= Potential_electrolyte_BC(p,c,dx,i,D)
    len=length(c);
    if i==len
        dxm=dx(i-1);dxc=dx(i);dxp=dx(i);
        cm=c(i-1);cc=c(i);cp=c(i);
        pm=p(i-1);pc=p(i);pp=p(i);
        Dm=D(i-1);Dc=D(i);Dp=D(i);
    elseif i==1
        dxm=dx(i);dxc=dx(i);dxp=dx(i+1);
        cm=c(i);cc=c(i);cp=c(i+1);
        %pm=p(i);pc=p(i);pp=p(i+1);
        pm=-p(i);pc=p(i);pp=p(i+1); % DEBUG BEN This BC forces potential to 0 at x=0. It is incorrect but completes closure of the system. This is temporary
        Dm=D(i);Dc=D(i);Dp=D(i+1);
    else
        dxm=dx(i-1);dxc=dx(i);dxp=dx(i+1);
        cm=c(i-1);cc=c(i);cp=c(i+1);
        pm=p(i-1);pc=p(i);pp=p(i+1);
        Dm=D(i-1);Dc=D(i);Dp=D(i+1);
    end
end

function [dxm,dxc,dxp,pm,pc,pp,Dm,Dc,Dp]= Potential_solid_BC(p,dx,i,D,separator_index,source_BC)
    len=length(p);
    if i==len
        dxm=dx(i-1);dxc=dx(i);dxp=dx(i);
        pm=p(i-1);pc=p(i);pp=p(i)+source_BC(2)*dxc;
        Dm=D(i-1);Dc=D(i);Dp=D(i);
    elseif i==1
        dxm=dx(i);dxc=dx(i);dxp=dx(i+1);
        %pm=p(i)-source_BC(1)*dxc;pc=p(i);pp=p(i+1); % This boundary is correct but apparently redundant. 
        pm=-p(i);pc=p(i);pp=p(i+1);   % DEBUG BEN This BC forces potential to 0 at x=0. Which appears to be correct (see Plett and Qiao).
        Dm=D(i);Dc=D(i);Dp=D(i+1);
    elseif i==separator_index
        dxm=dx(i-1);dxc=dx(i);dxp=dx(i);
        pm=p(i-1);pc=p(i);pp=p(i);
        Dm=D(i-1);Dc=D(i);Dp=D(i);
    elseif i==separator_index+1
        dxm=dx(i);dxc=dx(i);dxp=dx(i+1);
        %pm=p(i);pc=p(i);pp=p(i+1);
        pm=p(i-1);pc=p(i);pp=p(i+1); %DEBUG BEN This BC connects potential in the positive electrode to the negative electrode
        Dm=D(i);Dc=D(i);Dp=D(i+1);
    else
        dxm=dx(i-1);dxc=dx(i);dxp=dx(i+1);
        pm=p(i-1);pc=p(i);pp=p(i+1);
        Dm=D(i-1);Dc=D(i);Dp=D(i+1);
    end
end

function fc = LHS_f_ce (c,D,dx,source)
    len=length(c);
    fc=zeros(len,1);
    
    for i = 1:1:len
        fc(i)=LHS_f_ce_single_cell (c,D,dx,source,i);     
    end
end

function f = LHS_f_ce_single_cell (c,D,dx,source,i)
    [dxm,dxc,dxp,dummym,dummyc,dummyp,cm,cc,cp,Dm,Dc,Dp]=concentration_electrolyte_BC(c,c,dx,i,D);
            
    f =  2*(Dp*dxp+Dc*dxc)*(cp-cc)/((dxc+dxp)^2) +  ...
        2*(Dm*dxm+Dc*dxc)*(cm-cc)/((dxc+dxm)^2)+ dxc*source(i);        
end


function Jac_fc = LHS_Jac_f_Fdiff_ce(c,D,dx,source)
    dc_prop=1.0001;
    len=length(c);
    Jac_fc=zeros(len);
    
    for i = 1:1:len        
        fc   =  LHS_f_ce_single_cell (c,D,dx,source,i);
        for j = max(1,i-1):1:min(len,i+1)
            cpdc=c;
            cpdc(j)=cpdc(j)*dc_prop;
            fc_diff   =  LHS_f_ce_single_cell (cpdc,D,dx,source,i);
            %disp(fc_diff)
            Jac_fc(i,j)= (fc_diff-fc)/((c(j)*(dc_prop-1)));
        end
        
    end
        
end


function f = LHS_f_pe_single_cell (p,c,kappa,dx,source,i)
    [dxm,dxc,dxp,pm,pc,pp,cm,cc,cp,km,kc,kp]=Potential_electrolyte_BC(p,c,dx,i,kappa);
            
    f =  2*(kp*dxp+kc*dxc)*(pp-pc)/((dxc+dxp)^2) +  ...
        2*(km*dxm+kc*dxc)*(pm-pc)/((dxc+dxm)^2) + ...
        2*(kp*dxp+kc*dxc)*(log(cp)-log(cc))/((dxc+dxp)^2) +  ...
        2*(km*dxm+kc*dxc)*(log(cm)-log(cc))/((dxc+dxm)^2) + ...
        dxc*source(i); 

        %if i==2
        %    disp(num2str(f)+", "+num2str((km*dxm+kc*dxc)*(pm-pc))+", "+num2str((km*dxm+kc*dxc))+", "+num2str((km))+", "+num2str((dxc)))
        %end       
end

function fp = LHS_f_pe (p,c,kappa,dx,source)
    len=length(p);
    fp=zeros(len,1);
    
    for i = 1:1:len
        fp(i)=LHS_f_pe_single_cell (p,c,kappa,dx,source,i);     
    end
end

function Jac_fp = LHS_Jac_f_Fdiff_pe(p,c,D,dx,source)
    dp_prop=1.0001;
    len=length(p);
    Jac_fp=zeros(len);
    
    for i = 1:1:len        
        fp   =  LHS_f_pe_single_cell (p,c,D,dx,source,i);
        for j = max(1,i-1):1:min(len,i+1)
            ppdp=p;
            ppdp(j)=ppdp(j)*dp_prop;
            fp_diff   =  LHS_f_pe_single_cell (ppdp,c,D,dx,source,i);
            %disp(fp_diff)
            Jac_fp(i,j)= (fp_diff-fp)/((p(j)*(dp_prop-1)));
        end
        
    end
        
end


function f = LHS_f_ps_single_cell (p,c,sigma,dx,source,i,separator_index,source_BC)
    [dxm,dxc,dxp,pm,pc,pp,sigm,sigc,sigp]=Potential_solid_BC(p,dx,i,sigma,separator_index,source_BC);
            
    f =  2*(sigp*dxp+sigc*dxc)*(pp-pc)/((dxc+dxp)^2) +  ...
        2*(sigm*dxm+sigc*dxc)*(pm-pc)/((dxc+dxm)^2) + ...
        dxc*source(i); 
     
end

function fp = LHS_f_ps (p,c,sigma,dx,source,separator_index,source_BC)
    len=length(p);
    fp=zeros(len,1);
    
    for i = 1:1:len
        fp(i)=LHS_f_ps_single_cell (p,c,sigma,dx,source,i,separator_index,source_BC);     
    end
end

function Jac_fp = LHS_Jac_f_Fdiff_ps(p,c,D,dx,source,separator_index,source_BC)
    dp_prop=1.0001;
    len=length(p);
    Jac_fp=zeros(len);
    
    for i = 1:1:len        
        fp   =  LHS_f_ps_single_cell (p,c,D,dx,source,i,separator_index,source_BC);
        for j = max(1,i-1):1:min(len,i+1)
            ppdp=p;
            ppdp(j)=ppdp(j)*dp_prop;
            fp_diff   =  LHS_f_ps_single_cell (ppdp,c,D,dx,source,i,separator_index,source_BC);
            %disp(fp_diff)
            Jac_fp(i,j)= (fp_diff-fp)/((p(j)*(dp_prop-1)));
        end
        
    end
        
end



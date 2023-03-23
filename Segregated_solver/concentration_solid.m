global BV_fun

if deb.prints==1
    disp('Solving diffusion in particles');
end
% Following the method described in the article by Zeng, Y., P. Albertus, 
%R. Klein, N. Chaturvedi, A. Kojic, M. Z. Bazant,
%and J. Christensen. “Efficient Conservative Numerical Schemes
%for 1D Nonlinear Spherical Diffusion Equations with Applications
%in Battery Modeling.” Journal of the Electrochemical Society

j_temp=-5;%-5.35e-5; %DEBUG BEN temporary hard coded value for development purposes

temporary_j=0;
if temporary_j==1
    jn=5;%-5.35e-5; %DEBUG BEN temporary hard coded value for development purposes
    jp=5;%-5.35e-5; %DEBUG BEN temporary hard coded value for development purposes
    source_n=p.A_s_n*jn*ones(1,sol.nb_cell_n);
    source_s=0*ones(1,sol.nb_cell_s);
    source_p=p.A_s_p*jp*ones(1,sol.nb_cell_p);
else

    fv.j=butler_volmer_eq(fv.pe,fv.ps,fv.ce,fv.cse,p.k0,p.alpha,p.Faraday,p.Rg, ini.T0,fv.Ueq,p.Rfilm,p.csn_max,p.csp_max,sol.nb_cell_n,sol.nb_cell_s);
    source_n=p.A_s_n*fv.j(1:sol.nb_cell_n);
    source_s=0*ones(1,sol.nb_cell_s);
    source_p=p.A_s_p*fv.j(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell);
end

if deb.prints==1
    disp("DEBUG BEN j for csn csp")
    disp(fv.j)
end

%%create left hand side matrices M
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

if exist('Mp_part_diff','var') == 0
    Volumes_n=ones(1,max_ind);
    Volumes_p=ones(1,max_ind);
    for ind=1:1:max_ind
        Volumes_n(ind)=sol.part_coord_n(ind)^2*drn(ind) + drn(ind)^3/12.0;
        Volumes_p(ind)=sol.part_coord_p(ind)^2*drp(ind) + drp(ind)^3/12.0;
    end

    Mp_part_diff=calculate_mass_matrix(max_ind,Volumes_p);
    Mn_part_diff=calculate_mass_matrix(max_ind,Volumes_n);
end

%% Initialize the vector for concentration at the next time step: (c(t+dt)=c(t))
csn_next=fv.csn;
csp_next=fv.csp;

%%Calculate the left hand side vector for the concentration at time t: f(c) and solve the system

for cross_cell_ind = 1:1:sol.nb_cell_n
    [csn_next(:,cross_cell_ind),csn_next_save] = Newton_concentration(csn_next(:,cross_cell_ind),p.Dsn,sol.part_coord_n, ...
                            fv.j(cross_cell_ind),fv.csn(:,cross_cell_ind),sol.newton_meth_max_ite,Mn_part_diff,sol.dt,sol.newton_meth_res_threshold,2);
end
for cross_cell_ind = 1:1:sol.nb_cell_p
    [csp_next(:,cross_cell_ind),csp_next_save] = Newton_concentration(csp_next(:,cross_cell_ind),p.Dsp,sol.part_coord_p, ...
                            fv.j(sol.nb_cell_n+sol.nb_cell_s+cross_cell_ind),fv.csp(:,cross_cell_ind),sol.newton_meth_max_ite,Mp_part_diff,sol.dt,sol.newton_meth_res_threshold,2);
end


fv.csn=csn_next;
fv.csp=csp_next;
fv.cse=cat(2,fv.csn(length(sol.part_coord_n),:),fv.csp(length(sol.part_coord_p),:));




clear outcome_anim csn_next csp_next csn_next_save csp_next_save j_temp ;
clear drn drp Volumes_n Volumes_p ind temporary_j;



%% Functions

function M = calculate_mass_matrix(max_ind,Volumes)
    M1=diag(3.0/4.0*ones(max_ind,1)) + diag(1.0/8.0*ones(max_ind-1,1),1) + diag(1.0/8.0*ones(max_ind-1,1),-1);
    M1(1,2)=1/4;
    M1(max_ind,max_ind-1)=1/4;
    M2=diag(Volumes);
    M=M1*M2;
end

function [csn_next,csn_next_save] = Newton_concentration(csn_next,D,r,j,csn,newt_max_ite,Mn,dt,newt_lim,Jac_method)
    
    fn_c=LHS_f(csn,D,r,j,csn);
    csn_next_save= ones(length(csn_next),newt_max_ite);
    %% Start the newton method to solve g(c(t+dt))=0  (g(c(t+dt))=M(c(t+dt)-c(t))/dt -1/2 (f(t+dt)+f(t)))
    for newt_ite = 1:1:newt_max_ite
        len=length(csn_next);

        %%Calculate the jacobian matrix of the left hand side vector for the concentration at time t+dt: Jf(c)
        %if method=1, jacobian is calculated analytically, if method=2, jac
        %is calculated using finite differences.
        if Jac_method==1
            Jac_fn_c_next=LHS_Jac_f_analytic(csn_next,D,r,j,csn);
        elseif Jac_method==2
            Jac_fn_c_next=LHS_Jac_f_Fdiff(csn_next,D,r,j,csn);
        end
        fn_c_next=LHS_f(csn_next,D,r,j,csn);

        %disp(Jac_fn_c_next);
        %disp(Jac_fn_c_next_2-Jac_fn_c_next);
        %disp(transpose(fn_c_next));

        %%Compilation of the Crank Nicholson method s function and jacobian
        %g(c(t+dt))=M(c(t+dt)-c(t))/dt -1/2 (f(t+dt)+f(t))
        %Jac_g(c(t+dt))=M/dt -1/2 jac_f(t+dt)
        %Jac_g*(delta(c(t+dt)))=g

        gn=Mn*(csn_next-csn)/dt - 0.5*(fn_c_next+fn_c);

        Jac_gn=Mn/dt -0.5*Jac_fn_c_next;

        delta_cn = linsolve(Jac_gn,-gn);
        norm_delta_cn = sum(transpose(delta_cn) * delta_cn)/len;
        csn_next=csn_next+delta_cn;

        if norm_delta_cn<newt_lim
            %disp("Newton method for lithium concentration in solid has converged. ite: "+int2str(newt_ite)+" , residual: "+num2str(norm_delta_cn))
            break
        elseif newt_ite==newt_max_ite
            disp("Newton method for lithium concentration in solid did not converge. ite: "+int2str(newt_ite)+" , residual: "+num2str(norm_delta_cn))
        end
        csn_next_save(:,newt_ite)= csn_next;
    end
end

function [rm,rp,cm,cp]= particle_diff_BC(c,r,i,j,D,cminus)
    len=length(c);
    if i==len
        rm=r(i-1);
        rp=r(len)+ r(len) -r(len-1);
        cm=c(i-1);
        cp=cminus(len) - j/D*(r(len)-r(len-1));
    elseif i==1
        rm=-r(i+1);
        rp=r(i+1);
        cm=c(i+1);
        cp=c(i+1);
    else
        rm=r(i-1);
        rp=r(i+1);
        cm=c(i-1);
        cp=c(i+1);
    end
end

function fc = LHS_f(c,D,r,j,cminus)
    len=length(c);
    fc=zeros(len,1);
    
    for i = 1:1:len
        [rm,rp,cm,cp]=particle_diff_BC(c,r,i,j,D,cminus);
        
        factor_plus = D * ((rp+r(i))/2)^2 /(rp-r(i));
        factor_minus= D * ((r(i)+rm)/2)^2 /(r(i)-rm);
            
        fc(i)=     cp*factor_plus + cm*factor_minus - c(i)*(factor_plus+factor_minus);        
    end
    %disp([cp,(r(len)-rm),factor_plus]);
end


function Jac_fc = LHS_Jac_f_analytic(c,D,r,j,cminus)
    len=length(c);
    Jac_fc=zeros(len);
    
    for i=1:len
        [rm,rp,cm,cp]=particle_diff_BC(c,r,i,j,D,cminus);
        
        factor_plus = D * ((rp+r(i))/2)^2 /(rp-r(i));
        factor_minus= D * ((r(i)+rm)/2)^2 /(r(i)-rm);
        Jac_fc(i,i)= - (factor_plus+factor_minus);
        if i~= len
            Jac_fc(i,i+1)= (factor_plus);
        end
        if i~=1
            Jac_fc(i,i-1)= (factor_minus);
        end
    end
   
end

function Jac_fc = LHS_Jac_f_Fdiff(c,D,r,j,cminus)
    dc_prop=1.0001;
    len=length(c);
    Jac_fc=zeros(len);
    
    for i = 1:1:len
        [rm,rp,cm,cp]=particle_diff_BC(c,r,i,j,D,cminus);
        
        factor_plus = D * ((rp+r(i))/2)^2 /(rp-r(i));
        factor_minus= D * ((r(i)+rm)/2)^2 /(r(i)-rm);
            
        fc=     cp*factor_plus + cm*factor_minus - c(i)*(factor_plus+factor_minus);
        fcpdc=  cp*factor_plus + cm*factor_minus - (c(i)*dc_prop)*(factor_plus+factor_minus);    
        Jac_fc(i,i)= (fcpdc-fc)/((c(i)*(dc_prop-1)));
        
        if i~=len
            fcpdcp= (cp*dc_prop)*factor_plus +cm*factor_minus - c(i)*(factor_plus+factor_minus);
            Jac_fc(i,i+1)= (fcpdcp-fc)/((cp*(dc_prop-1)));
        end
        if i~=1
            fcpdcm= cp*factor_plus +(cm*dc_prop)*factor_minus -c(i)*(factor_plus+factor_minus);
            Jac_fc(i,i-1)= (fcpdcm-fc)/((cm*(dc_prop-1)));
        end
        
    end
        
end



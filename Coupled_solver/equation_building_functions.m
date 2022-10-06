classdef equation_building_functions
    methods

        function fc = LHS_f_cs (obj,c,D,r,source)
            len=length(c);
            lenr=length(r);
            fc=zeros(len,1);
            
            for i = 1:1:len
                i_rcell=mod(i,lenr);
                if i_rcell==0
                    i_rcell=lenr;
                end
                i_xcell=ceil(i/lenr);
                %disp("DEBUG BEN "+num2str(i_xcell)+"  and "+num2str(i)+"  and "+num2str(i_rcell)+"  and "+num2str(lenr))
                sour=source(i_xcell);
                fc(i)=obj.LHS_f_cs_single_cell(c,D,r,sour,i,i_rcell);     
            end
        end

        function f = LHS_f_cs_single_cell (obj,c,D,r,source,i,i_rcell)
            global BC_fun
            [rm,rc,rp,cm,cc,cp]=BC_fun.concentration_solid_BC(c,r,i,source,D,i_rcell);
                    
            factor_plus = D * ((rp+rc)/2)^2 /(rp-rc);
            factor_minus= D * ((rc+rm)/2)^2 /(rc-rm);
            f= cp*factor_plus + cm*factor_minus - cc*(factor_plus+factor_minus);       
        end

        function Jac_fc = LHS_Jac_f_Fdiff_cs(obj,c,D,r,source)
            dc_prop=1.0001;
            len=length(c);
            lenr=length(r);
            lenc=len/lenr;
            Jac_fc=zeros(len);
            
            for i = 1:1:len   
                i_rcell=mod(i,lenr);
                if i_rcell==0
                    i_rcell=lenr;
                end
                i_xcell=ceil(i/lenr);
                sour=source(i_xcell);

                fc   =  obj.LHS_f_cs_single_cell (c,D,r,sour,i,i_rcell);
                for ii = max(1,i-1):1:min(len,i+1)
                    cpdc=c;
                    %cpdc(ii)=cpdc(ii)*dc_prop;
                    %fc_diff   =  obj.LHS_f_cs_single_cell (cpdc,D,r,sour,i,i_rcell);
                    %Jac_fc(i,ii)= (fc_diff-fc)/((c(ii)*(dc_prop-1)));

                    if not(c(ii)==0)
                        cpdc(ii)=cpdc(ii)*dc_prop;
                        fc_diff   =  obj.LHS_f_cs_single_cell (cpdc,D,r,sour,i,i_rcell);
                        Jac_fc(i,ii)= (fc_diff-fc)/((c(ii)*(dc_prop-1)));
                    else
                        cpdc(ii)=0.000001;
                        fc_diff   =  obj.LHS_f_cs_single_cell (cpdc,D,r,sour,i,i_rcell);
                        Jac_fc(i,ii)= (fc_diff-fc)/(0.000001);
                    end

                end
                
            end
                
        end




        function fc = LHS_f_ce (obj,c,D,dx,source,ID)
            len=length(c);
            fc=zeros(len,1);
            
            for i = 1:1:len
                fc(i)=obj.LHS_f_ce_single_cell (c,D,dx,source,i);     
            end
        end

        function f = LHS_f_ce_single_cell (obj,c,D,dx,source,i)
            global BC_fun
            [dxm,dxc,dxp,dummym,dummyc,dummyp,cm,cc,cp,Dm,Dc,Dp]=BC_fun.concentration_electrolyte_BC(c,c,dx,i,D);
                    
            f =  2*(Dp*dxp+Dc*dxc)*(cp-cc)/((dxc+dxp)^2) +  ...
                2*(Dm*dxm+Dc*dxc)*(cm-cc)/((dxc+dxm)^2)+ dxc*source(i);        
        end

        function Jac_fc = LHS_Jac_f_Fdiff_ce(obj,c,D,dx,source)
            dc_prop=1.0001;
            len=length(c);
            Jac_fc=zeros(len);
            
            for i = 1:1:len        
                fc   =  obj.LHS_f_ce_single_cell (c,D,dx,source,i);
                for j = max(1,i-1):1:min(len,i+1)
                    cpdc=c;

                    if not(c(j)==0)
                        cpdc(j)=cpdc(j)*dc_prop;
                        fc_diff   =  obj.LHS_f_ce_single_cell (cpdc,D,dx,source,i);
                        Jac_fc(i,j)= (fc_diff-fc)/((c(j)*(dc_prop-1)));
                    else
                        cpdc(j)=0.000001;
                        fc_diff   =  obj.LHS_f_ce_single_cell (cpdc,D,dx,source,i);
                        Jac_fc(i,j)= (fc_diff-fc)/(0.000001);
                    end

                end
                
            end
                
        end



        function f = LHS_f_pe_single_cell (obj,p,c,kappa,kappa_D,dx,source,i)
            global BC_fun
            [dxm,dxc,dxp,pm,pc,pp,cm,cc,cp,km,kc,kp,kdm,kdc,kdp]=BC_fun.Potential_electrolyte_BC(p,c,dx,i,kappa,kappa_D);
                    
            f =  2*(kp*dxp+kc*dxc)*(pp-pc)/((dxc+dxp)^2) +  ...
                2*(km*dxm+kc*dxc)*(pm-pc)/((dxc+dxm)^2) + ...
                2*(kdp*dxp/cp+kdc*dxc/cc)*(cp-cc)/((dxc+dxp)^2) +  ...
                2*(kdm*dxm/cm+kdc*dxc/cc)*(cm-cc)/((dxc+dxm)^2) + ...
                dxc*source(i);      

                %2*(kdp*dxp+kdc*dxc)*(log(max(cp,0.0000001))-log(max(cc,0.0000001)))/((dxc+dxp)^2) +  ...
                %2*(kdm*dxm+kdc*dxc)*(log(max(cm,0.0000001))-log(max(cc,0.0000001)))/((dxc+dxm)^2) + ...

        end

        function fp = LHS_f_pe (obj,p,c,kappa,kappa_D,dx,source,ID)
            len=length(p);
            fp=zeros(len,1);
            
            for i = 1:1:len
                fp(i)=obj.LHS_f_pe_single_cell (p,c,kappa,kappa_D,dx,source,i);     
            end
        end

        function Jac_fp = LHS_Jac_f_Fdiff_pe(obj,p,c,kappa,kappa_D,dx,source)
            dp_prop=1.0001;
            len=length(p);
            Jac_fp=zeros(len);
            
            for i = 1:1:len        
                fp   =  obj.LHS_f_pe_single_cell (p,c,kappa,kappa_D,dx,source,i);
                for j = max(1,i-1):1:min(len,i+1)
                    ppdp=p;
                    
                    if not(p(j)==0)
                        ppdp(j)=ppdp(j)*dp_prop;
                        fp_diff   =  obj.LHS_f_pe_single_cell (ppdp,c,kappa,kappa_D,dx,source,i);
                        Jac_fp(i,j)= (fp_diff-fp)/((p(j)*(dp_prop-1)));
                    else
                        ppdp(j)=0.000001;
                        fp_diff   =  obj.LHS_f_pe_single_cell (ppdp,c,kappa,kappa_D,dx,source,i);
                        Jac_fp(i,j)= (fp_diff-fp)/(0.000001);
                    end
                end
                
            end
                
        end


        function f = LHS_f_ps_single_cell (obj,p,c,sigma,dx,source,i,separator_index,source_BC,print_info)
            global BC_fun
            [dxm,dxc,dxp,pm,pc,pp,sigm,sigc,sigp]=BC_fun.Potential_solid_BC(p,dx,i,sigma,separator_index,source_BC);
                    
            f =  2*(sigp*dxp+sigc*dxc)*(pp-pc)/((dxc+dxp)^2) +  ...
                2*(sigm*dxm+sigc*dxc)*(pm-pc)/((dxc+dxm)^2) + ...
                dxc*source(i); 
            if i==1 && print_info==1
                disp("DEBUG BEN LHS_f_ps_single_cell "+num2str(f)+"  "+num2str(dxc*source(i))+"  "+num2str(2*(sigp*dxp+sigc*dxc)*(pp-pc)/((dxc+dxp)^2) +2*(sigm*dxm+sigc*dxc)*(pm-pc)/((dxc+dxm)^2)))
                disp("DEBUG BEN LHS_f_ps_single_cell "+num2str(pc)+"  "+num2str(pp)+"  "+num2str(pm))
                disp("DEBUG BEN LHS_f_ps_single_cell "+num2str(2*(sigm*dxm+sigc*dxc)*(pm-pc)/((dxc+dxm)^2))+"  "+num2str((sigm*dxm+sigc*dxc))+"  "+num2str((dxc+dxm)^2))
            end
        end

        function fp = LHS_f_ps (obj,p,c,sigma,dx,source,separator_index,source_BC,ID)
            len=length(p);
            fp=zeros(len,1);
            
            for i = 1:1:len
                fp(i)=obj.LHS_f_ps_single_cell (p,c,sigma,dx,source,i,separator_index,source_BC,1);     
            end
            disp("DEBUG BEN eq_build_fncts "+num2str(fp)+"  "+ID)
            disp("DEBUG BEN eq_build_fncts "+num2str(p)+"  "+num2str(transpose(source)))
            disp("DEBUG BEN eq_build_fncts "+"  "+num2str(sigma)+"  "+num2str(source_BC))
        end

        function Jac_fp = LHS_Jac_f_Fdiff_ps(obj,p,c,D,dx,source,separator_index,source_BC)
            dp_prop=1.0001;
            len=length(p);
            Jac_fp=zeros(len);
            
            for i = 1:1:len        
                fp   =  obj.LHS_f_ps_single_cell (p,c,D,dx,source,i,separator_index,source_BC,0);
                for j = max(1,i-1):1:min(len,i+1)
                    ppdp=p;
                    if not(p(j)==0)
                        ppdp(j)=ppdp(j)*dp_prop;
                        fp_diff   =  obj.LHS_f_ps_single_cell (ppdp,c,D,dx,source,i,separator_index,source_BC,0);
                        Jac_fp(i,j)= (fp_diff-fp)/((p(j)*(dp_prop-1)));
                    else
                        ppdp(j)=0.000001;
                        fp_diff   =  obj.LHS_f_ps_single_cell (ppdp,c,D,dx,source,i,separator_index,source_BC,0);
                        Jac_fp(i,j)= (fp_diff-fp)/(0.000001);
                    end

                    if i == j & Jac_fp(i,j)==0
                        Jac_fp(i,j)=0.000000001;
                    end
                end
                
            end
                
        end

        function [De,kappa_eff, kappa_D_eff] = update_param_functions (obj,ce,cse) 
            global p
            global sol
            global fv
            global ini
            global param_functions

            De=zeros(sol.nb_cell,1);
            kappa_eff=zeros(sol.nb_cell,1);
            fv.Ueq=0*fv.Ueq;

            if p.De_function_mode==1
                count_loc=0;
                for cccc = 1:1:length(ce)
                    count_loc=count_loc+1;
                    De(count_loc) = param_functions.electrolyte_diffusivity(ce(cccc))*p.De_eff_coeff(count_loc);
                end
            end

            if p.kappa_function_mode==1
                count_loc=0;
                for cccc = 1:1:length(ce)
                    count_loc=count_loc+1;
                    kappa_eff(count_loc) = param_functions.electrolyte_conductivity(ce(cccc))*p.De_eff_coeff(count_loc);
                    kappa_D_eff(count_loc) = kappa_eff(count_loc)*2* p.Rg * ini.T0/p.Faraday * (p.t_plus-1);
                end 
            end

            if p.kappa_function_mode==1
                count_loc=0;
                for cccc = 1:1:sol.nb_cell_n
                    fv.Ueq(cccc) = param_functions.neg_electrode_Ueq(cse(cccc),cccc);
                end
                for cccc = sol.nb_cell_n+1:1:sol.nb_cell_n+sol.nb_cell_p
                    fv.Ueq(cccc+sol.nb_cell_s) = param_functions.pos_electrode_Ueq(cse(cccc),cccc);
                end
            end


        end

        function source = update_source_single_cell (obj,i,ID) 
            global p
            global sol
            global fv
            global ini
            global param_functions

            A_s=p.A_s_n;
            Ds= p.Dsn;
            if i>sol.nb_cell_n
                A_s=p.A_s_p;
                Ds= p.Dsp;
            end

            if ID=="ps"
                source=-p.Faraday*A_s*fv.j(i);
            elseif ID=="pe"
                source= p.Faraday*A_s*fv.j(i);
            elseif ID=="ce"
                source= (1-p.t_plus)*A_s*fv.j(i);
            elseif ID=="cs"
                source= -fv.j(i)/Ds;
            else
                disp("ERROR update_source_single_cell : ID string is incorrect.")
            end
        end


        function [source_pe,source_ps,source_ce,source_csn,source_csp] = update_sources (obj) 
            global p
            global sol
            global fv
            global ini
            global param_functions

            source_ps=zeros(1,sol.nb_cell-sol.nb_cell_s);
            source_pe=zeros(1,sol.nb_cell);
            source_ce=zeros(1,sol.nb_cell);
            source_cs=zeros(1,sol.nb_cell);
            

            for i = 1:1:sol.nb_cell_n
                source_ps(i)=obj.update_source_single_cell (i,"ps") ;
                source_pe(i)=obj.update_source_single_cell (i,"pe") ;
                source_ce(i)=obj.update_source_single_cell (i,"ce") ;
                source_cs(i)=obj.update_source_single_cell (i,"cs") ;
            end

            for i = sol.nb_cell_n+sol.nb_cell_s+1:1:sol.nb_cell
                source_ps(i-sol.nb_cell_s)=obj.update_source_single_cell (i,"ps") ;
                source_pe(i)=obj.update_source_single_cell (i,"pe") ;
                source_ce(i)=obj.update_source_single_cell (i,"ce") ;
                source_cs(i)=obj.update_source_single_cell (i,"cs") ;
            end

            source_csn=source_cs(1:sol.nb_cell_n);
            source_csp=source_cs(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell);

        end
    end
end



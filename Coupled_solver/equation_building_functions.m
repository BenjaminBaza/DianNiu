classdef equation_building_functions
    methods

        function ind = pe2ps_index(obj,i)
            global p
            global sol
            
            if i>sol.nb_cell_n
                ind=i-sol.nb_cell_s;
            else
                ind=i;
            end
        end

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
                sour=source(i_xcell);
                fc(i)=obj.LHS_f_cs_single_cell(c,D,r,sour,i,i_rcell,0);     
            end
        end

        function f = LHS_f_cs_single_cell (obj,c,D,r,source,i,i_rcell,print_info)
            global BC_fun
            [rm,rc,rp,cm,cc,cp]=BC_fun.concentration_solid_BC(c,r,i,source,D,i_rcell);
                    
            factor_plus = D * ((rp+rc)/2)^2 /(rp-rc);
            factor_minus= D * ((rc+rm)/2)^2 /(rc-rm);
            f= cp*factor_plus + cm*factor_minus - cc*(factor_plus+factor_minus);  
            if print_info ==1
                disp("DEBUG BEN LHS_f_cs_single_cell "+num2str(cp)+"  "+num2str(factor_plus)+"  "+num2str(cp*factor_plus)+"  "+num2str(cm*factor_minus)+"  "+num2str(cc*(factor_plus+factor_minus))+"  "+num2str(f))
            end     
        end

        function Jac_fcs = LHS_Jac_f_Fdiff_cs(obj,pe,ps,ce,cse,T,Ueq,resize_cs,D,r,source,electrode)
            global fv
            global sol
            global p
            global BV_fun
            mock_j = fv.j;

            dc_prop=1.0001;
            len=length(resize_cs);
            lenr=length(r);
            lenc=len/lenr;
            Jac_fcs=zeros(len);
            
            mock_source = source;
            mock_j = fv.j;
            mock_cse = cse;

            for i = 1:1:len   
                i_rcell=mod(i,lenr);
                if i_rcell==0
                    i_rcell=lenr;
                end
                i_xcell=ceil(i/lenr);
                sour=source(i_xcell);

                fcs   =  obj.LHS_f_cs_single_cell (resize_cs,D,r,sour,i,i_rcell,0);
                for ii = max(1,i-1):1:min(len,i+1)
                    if (i_rcell==lenr && ii==i+1)||(i_rcell==1 && ii==i-1)
                        continue
                    end

                    cpdc=resize_cs;
                    if not(cpdc(ii)==0)
                        divider=(cpdc(ii)*(dc_prop-1));
                        cpdc(ii)=cpdc(ii)*dc_prop;
                    else
                        cpdc(ii)=0.000001;
                        divider=0.000001;
                    end

                    if i==ii && i_rcell==lenr 
                        if electrode=="csn"
                            i_xfullcell=i_xcell;
                            i_xfullercell=i_xcell;
                            mock_cse(i_xfullcell)=min(cpdc(i),p.csn_max);
                        else
                            i_xfullcell=i_xcell+sol.nb_cell_n;
                            i_xfullercell=i_xcell+sol.nb_cell_n+sol.nb_cell_s;
                            mock_cse(i_xfullcell)=min(cpdc(i),p.csp_max);
                        end

                        mock_j=fv.j;
                        mock_j(i_xfullercell)=BV_fun.butler_volmer_singlecell_standalone(i_xfullercell,pe,ps,ce,mock_cse,Ueq);

                        %mock_j=BV_fun.butler_volmer_equation(pe,ps,ce,mock_cse,T,Ueq,"LHS_Jac_f_Fdiff_cs");
                        
                        sour=obj.update_source_single_cell(i_xfullercell,mock_j,"cs");     
                    end

                    fcs_diff   =  obj.LHS_f_cs_single_cell (cpdc,D,r,sour,i,i_rcell,0);
                    Jac_fcs(i,ii)= (fcs_diff-fcs)/(divider);
                    
                end
            end
        end


        function Jac_fcs = LHS_Jac_f_Fdiff_csdce(obj,pe,ps,ce,cse,T,Ueq,resize_cs,D,r,source,electrode)
            global fv
            global sol
            global p
            global BV_fun
            mock_j = fv.j;

            dc_prop=1.0001;
            len=length(resize_cs);
            lenr=length(r);
            lenc=len/lenr;
            Jac_fcs=zeros(len,sol.nb_cell);
            
            mock_source = source;
            mock_j = fv.j;
            mock_cse = cse;

            for i = 1:1:len   
                i_rcell=mod(i,lenr);
                if i_rcell==0
                    i_rcell=lenr;
                end
                if not(i_rcell==lenr)
                    continue
                end 

                i_xcell=ceil(i/lenr);
                if electrode=="csn"
                    i_xfullcell=i_xcell;
                else
                    i_xfullcell=i_xcell+sol.nb_cell_n+sol.nb_cell_s;
                end
                sour=source(i_xcell);
                sour_save=sour;

                fcs   =  obj.LHS_f_cs_single_cell (resize_cs,D,r,sour,i,i_rcell,0);

                cpdc=ce;
                if not(cpdc(i_xfullcell)==0)
                    divider=(cpdc(i_xfullcell)*(dc_prop-1));
                    cpdc(i_xfullcell)=cpdc(i_xfullcell)*dc_prop;
                else
                    cpdc(i_xfullcell)=0.000001;
                    divider=0.000001;
                end

                if  i_rcell==lenr 
                    mock_j=fv.j;
                    mock_j(i_xfullcell)=BV_fun.butler_volmer_singlecell_standalone(i_xfullcell,pe,ps,cpdc,cse,Ueq);

                    %mock_j=BV_fun.butler_volmer_equation(pe,ps,cpdc,cse,T,Ueq,"LHS_Jac_f_Fdiff_csdce"); 
                    sour=obj.update_source_single_cell(i_xfullcell,mock_j,"cs"); 
                end

                fcs_diff   =  obj.LHS_f_cs_single_cell (resize_cs,D,r,sour,i,i_rcell,0);
                Jac_fcs(i,i_xfullcell)= (fcs_diff-fcs)/(divider);
                %disp("DEBUG BEN LHS_Jac_f_Fdiff_cs_dce "+num2str(mock_j(i_xfullcell)-fv.j(i_xfullcell))+"  "+num2str(fv.j(i_xfullcell))+"  "+num2str(sour-sour_save)+"  "+num2str(sour)+"  "+num2str((fcs_diff-fcs)))

            
            end
        end

        function Jac_fcs = LHS_Jac_f_Fdiff_csdpe(obj,pe,ps,ce,cse,T,Ueq,resize_cs,D,r,source,electrode)
            global fv
            global sol
            global p
            global BV_fun
            mock_j = fv.j;

            dc_prop=1.0001;
            len=length(resize_cs);
            lenr=length(r);
            lenc=len/lenr;
            Jac_fcs=zeros(len,sol.nb_cell);
            
            mock_source = source;
            mock_j = fv.j;
            mock_cse = cse;

            for i = 1:1:len 
                i_rcell=mod(i,lenr);
                if i_rcell==0
                    i_rcell=lenr;
                end
                if not(i_rcell==lenr)
                    continue
                end  

                
                i_xcell=ceil(i/lenr);
                if electrode=="csn"
                    i_xfullcell=i_xcell;
                else
                    i_xfullcell=i_xcell+sol.nb_cell_n+sol.nb_cell_s;
                end
                sour=source(i_xcell);
                sour_save=sour;


                fcs   =  obj.LHS_f_cs_single_cell (resize_cs,D,r,sour,i,i_rcell,0);
                
                

                cpdc=pe;
                if not(cpdc(i_xfullcell)==0)
                    divider=(cpdc(i_xfullcell)*(dc_prop-1));
                    cpdc(i_xfullcell)=cpdc(i_xfullcell)*dc_prop;
                else
                    cpdc(i_xfullcell)=0.000001;
                    divider=0.000001;
                end

                if  i_rcell==lenr 
                    mock_j=fv.j;
                    mock_j(i_xfullcell)=BV_fun.butler_volmer_singlecell_standalone(i_xfullcell,cpdc,ps,ce,cse,Ueq);

                    %mock_j=BV_fun.butler_volmer_equation(cpdc,ps,ce,cse,T,Ueq,"LHS_Jac_f_Fdiff_csdce"); 
                    sour=obj.update_source_single_cell(i_xfullcell,mock_j,"cs"); 
                end

                fcs_diff   =  obj.LHS_f_cs_single_cell (resize_cs,D,r,sour,i,i_rcell,0);
                Jac_fcs(i,i_xfullcell)= (fcs_diff-fcs)/(divider);
            
            end
        end


        function Jac_fcs = LHS_Jac_f_Fdiff_csdps(obj,pe,ps,ce,cse,T,Ueq,resize_cs,D,r,source,electrode)
            global fv
            global sol
            global p
            global BV_fun
            mock_j = fv.j;

            dc_prop=1.0001;
            len=length(resize_cs);
            lenr=length(r);
            lenc=len/lenr;
            Jac_fcs=zeros(len,sol.nb_cell_n+sol.nb_cell_p);
            
            mock_source = source;
            mock_j = fv.j;
            mock_cse = cse;

            for i = 1:1:len 
                i_rcell=mod(i,lenr);
                if i_rcell==0
                    i_rcell=lenr;
                end
                if not(i_rcell==lenr)
                    continue
                end  

                
                i_xcell=ceil(i/lenr);
                if electrode=="csn"
                    i_xfullcell=i_xcell;
                    i_xfullercell=i_xcell;
                else
                    i_xfullcell=i_xcell+sol.nb_cell_n;
                    i_xfullercell=i_xcell+sol.nb_cell_n+sol.nb_cell_s;
                end
                sour=source(i_xcell);
                sour_save=sour;


                fcs   =  obj.LHS_f_cs_single_cell (resize_cs,D,r,sour,i,i_rcell,0);
                
                

                cpdc=ps;
                if not(cpdc(i_xfullcell)==0)
                    divider=(cpdc(i_xfullcell)*(dc_prop-1));
                    cpdc(i_xfullcell)=cpdc(i_xfullcell)*dc_prop;
                else
                    cpdc(i_xfullcell)=0.000001;
                    divider=0.000001;
                end

                if  i_rcell==lenr 
                    mock_j=fv.j;
                    mock_j(i_xfullercell)=BV_fun.butler_volmer_singlecell_standalone(i_xfullercell,pe,cpdc,ce,cse,Ueq);

                    %mock_j=BV_fun.butler_volmer_equation(pe,cpdc,ce,cse,T,Ueq,"LHS_Jac_f_Fdiff_csdce"); 
                    
                    sour=obj.update_source_single_cell(i_xfullercell,mock_j,"cs"); 
                end

                fcs_diff   =  obj.LHS_f_cs_single_cell (resize_cs,D,r,sour,i,i_rcell,0);
                Jac_fcs(i,i_xfullcell)= (fcs_diff-fcs)/(divider);
            
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

        function Jac_fce = LHS_Jac_f_Fdiff_ce(obj,pe,ps,ce,cse,T,Ueq,D,dx,source)
            dc_prop=1.0001;
            len=length(ce);
            Jac_fce=zeros(len);
            
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            for i = 1:1:len        
                fce   =  obj.LHS_f_ce_single_cell (ce,D,dx,source,i);
                for j = max(1,i-1):1:min(len,i+1)
                    cpdc=ce;
                    if not(ce(j)==0)
                        cpdc(j)=cpdc(j)*dc_prop;
                        divider=(ce(j)*(dc_prop-1));
                    else
                        cpdc(j)=0.000001;
                        divider=0.000001;
                    end

                    if i==j 
                        if i<=sol.nb_cell_n
                            csmax=p.csn_max;
                        else
                            csmax=p.csp_max;
                        end
                        %DEBUG BEN new calculation of mock_j

                        mock_j=fv.j;
                        mock_j(i)=BV_fun.butler_volmer_singlecell_standalone(i,pe,ps,cpdc,cse,Ueq);

                        %mock_j=BV_fun.butler_volmer_equation(pe,ps,cpdc,cse,T,Ueq,"LHS_Jac_f_Fdiff_ce");
                        %DEBUG BEN new calculation of mock_j

                        mock_source = source;
                        mock_source(i)=obj.update_source_single_cell(i,mock_j,"ce");
                    
                        source_used=mock_source;
                    else
                        source_used=source;
                    end

                    fce_diff   =  obj.LHS_f_ce_single_cell (cpdc,D,dx,source_used,i);
                    Jac_fce(i,j)= (fce_diff-fce)/(divider);

                    if i == j & Jac_fce(i,j)==0
                        Jac_fce(i,j)=0.000000001;
                    end
                end
            end
        end


        function Jac_dfce_dpe = LHS_Jac_f_Fdiff_cedpe(obj,pe,ps,ce,cse,T,Ueq,D,dx,source)
            dc_prop=1.0001;
            len=length(ce);
            Jac_dfce_dpe=zeros(len);
            
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            for i = 1:1:len        
                fce   =  obj.LHS_f_ce_single_cell (ce,D,dx,source,i);
                for j = max(1,i):1:min(len,i)
                    cpdc=pe;
                    if not(pe(j)==0)
                        cpdc(j)=cpdc(j)*dc_prop;
                        divider=(pe(j)*(dc_prop-1));
                    else
                        cpdc(j)=0.000001;
                        divider=0.000001;
                    end

                    if i==j 
                        if i<=sol.nb_cell_n
                            csmax=p.csn_max;
                        else
                            csmax=p.csp_max;
                        end

                        mock_j=fv.j;
                        mock_j(i)=BV_fun.butler_volmer_singlecell_standalone(i,cpdc,ps,ce,cse,Ueq);

                        %mock_j=BV_fun.butler_volmer_equation(cpdc,ps,ce,cse,T,Ueq,"LHS_Jac_f_Fdiff_ce");
                        
                        mock_source = source;
                        mock_source(i)=obj.update_source_single_cell(i,mock_j,"ce");
                    
                        source_used=mock_source;
                    else
                        source_used=source;
                    end

                    fce_diff   =  obj.LHS_f_ce_single_cell (ce,D,dx,source_used,i);
                    Jac_dfce_dpe(i,j)= (fce_diff-fce)/(divider);
                end
            end
        end


        function Jac_dfce_dps = LHS_Jac_f_Fdiff_cedps(obj,pe,ps,ce,cse,T,Ueq,D,dx,source)            
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            dc_prop=1.0001;
            len=length(ce);
            Jac_dfce_dps=zeros(len,sol.nb_cell_n+sol.nb_cell_p);

            for i = 1:1:len   
                if i>sol.nb_cell_n && i<=(sol.nb_cell_n+sol.nb_cell_s)
                    continue
                end  

                fce   =  obj.LHS_f_ce_single_cell (ce,D,dx,source,i);

                for j = max(1,i):1:min(len,i)
                    j_eff=j;
                    if j>sol.nb_cell_n
                        j_eff=j-sol.nb_cell_s;
                    end
 
                    cpdc=ps;
                    if not(ps(j_eff)==0)
                        cpdc(j_eff)=cpdc(j_eff)*dc_prop;
                        divider=(ps(j_eff)*(dc_prop-1));
                    else
                        cpdc(j_eff)=0.000001;
                        divider=0.000001;
                    end

                    if i<=sol.nb_cell_n
                        csmax=p.csn_max;
                    else
                        csmax=p.csp_max;
                    end
                    mock_j=fv.j;
                    mock_j(i)=BV_fun.butler_volmer_singlecell_standalone(i,pe,cpdc,ce,cse,Ueq);
                    %mock_j=BV_fun.butler_volmer_equation(pe,cpdc,ce,cse,T,Ueq,"LHS_Jac_f_Fdiff_cedps");
                    mock_source = source;
                    mock_source(i)=obj.update_source_single_cell(i,mock_j,"ce");
                    
                    source_used=mock_source;
                    
                    fce_diff   =  obj.LHS_f_ce_single_cell (ce,D,dx,source_used,i);
                    Jac_dfce_dps(i,j_eff)= (fce_diff-fce)/(divider);
                end
            end
        end


        function Jac_dfce_dcs = LHS_Jac_f_Fdiff_cedcs(obj,pe,ps,ce,cse,T,Ueq,D,dx,source)            
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            dc_prop=1.0001;
            len=length(ce);
            Jac_dfce_dcs=zeros(len,(sol.nb_cell_n+sol.nb_cell_p)*(sol.part_nb_cell+1));

            for i = 1:1:len   
                if i>sol.nb_cell_n && i<=(sol.nb_cell_n+sol.nb_cell_s)
                    continue
                end  

                fce   =  obj.LHS_f_ce_single_cell (ce,D,dx,source,i);

                for j = max(1,i):1:min(len,i)
                    csmax=p.csn_max;
                    j_eff=j;
                    if j>sol.nb_cell_n
                        csmax=p.csp_max;
                        j_eff=j-sol.nb_cell_s;
                    end
 
                    cpdc=cse;
                    if not(cse(j_eff)==0)
                        cpdc(j_eff)=min(cpdc(j_eff)*dc_prop,csmax);
                        divider=(cse(j_eff)*(dc_prop-1));
                    else
                        cpdc(j_eff)=0.000001;
                        divider=0.000001;
                    end

                    mock_j=fv.j;
                    mock_j(i)=BV_fun.butler_volmer_singlecell_standalone(i,pe,ps,ce,cpdc,Ueq);
                    
                    %mock_j=BV_fun.butler_volmer_equation(pe,ps,ce,cpdc,T,Ueq,"LHS_Jac_f_Fdiff_cedcs");
                    
                    mock_source = source;
                    mock_source(i)=obj.update_source_single_cell(i,mock_j,"ce");
                    
                    source_used=mock_source;
                    
                    fce_diff   =  obj.LHS_f_ce_single_cell (ce,D,dx,source_used,i);
                    Jac_dfce_dcs(i,j_eff*(sol.part_nb_cell+1))= (fce_diff-fce)/(divider);
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

        function Jac_fpe = LHS_Jac_f_Fdiff_pe(obj,pe,ps,ce,cse,T,Ueq,kappa,kappa_D,dx,source)
            dp_prop=1.0001;
            len=length(pe);
            Jac_fpe=zeros(len);
            
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            for i = 1:1:len        
                fpe   =  obj.LHS_f_pe_single_cell (pe,ce,kappa,kappa_D,dx,source,i);
                
                for j = max(1,i-1):1:min(len,i+1)
                    ppdp=pe;
                    if not(pe(j)==0)
                        ppdp(j)=ppdp(j)*dp_prop;
                        divider=(pe(j)*(dp_prop-1));
                    else
                        ppdp(j)=0.000001;
                        divider=0.000001;
                    end

                    if i==j 
                        if i<=sol.nb_cell_n
                            csmax=p.csn_max;
                        else
                            csmax=p.csp_max;
                        end
                        
                        mock_j=fv.j;
                        mock_j(i)=BV_fun.butler_volmer_singlecell_standalone(i,ppdp,ps,ce,cse,Ueq);

                        %mock_j=BV_fun.butler_volmer_equation(ppdp,ps,ce,cse,T,Ueq,"LHS_Jac_f_Fdiff_pe");
                        mock_source = source;
                        mock_source(i)=obj.update_source_single_cell(i,mock_j,"pe");
                    
                        source_used=mock_source;
                    else
                        source_used=source;
                    end

                    fpe_diff   =  obj.LHS_f_pe_single_cell (ppdp,ce,kappa,kappa_D,dx,source_used,i);
                    Jac_fpe(i,j)= (fpe_diff-fpe)/(divider);

                    if i == j & Jac_fpe(i,j)==0
                        Jac_fpe(i,j)=0.000000001;
                    end
                end
            end
        end


        function Jac_fpe = LHS_Jac_f_Fdiff_pedce(obj,pe,ps,ce,cse,T,Ueq,kappa,kappa_D,dx,source)
            dp_prop=1.0001;
            len=length(pe);
            Jac_fpe=zeros(len);
            
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            for i = 1:1:len        
                fpe   =  obj.LHS_f_pe_single_cell (pe,ce,kappa,kappa_D,dx,source,i);
                
                for j = max(1,i-1):1:min(len,i+1)
                    ppdp=ce;
                    if not(ce(j)==0)
                        ppdp(j)=ppdp(j)*dp_prop;
                        divider=(ce(j)*(dp_prop-1));
                    else
                        ppdp(j)=0.000001;
                        divider=0.000001;
                    end

                    if i==j 
                        if i<=sol.nb_cell_n
                            csmax=p.csn_max;
                        else
                            csmax=p.csp_max;
                        end

                        mock_j=fv.j;
                        mock_j(i)=BV_fun.butler_volmer_singlecell_standalone(i,pe,ps,ppdp,cse,Ueq);

                        %mock_j=BV_fun.butler_volmer_equation(pe,ps,ppdp,cse,T,Ueq,"LHS_Jac_fpe_Fdiff_ce");
                        mock_source = source;
                        mock_source(i)=obj.update_source_single_cell(i,mock_j,"pe");
                    
                        source_used=mock_source;
                    else
                        source_used=source;
                    end

                    fpe_diff   =  obj.LHS_f_pe_single_cell (pe,ppdp,kappa,kappa_D,dx,source_used,i);
                    Jac_fpe(i,j)= (fpe_diff-fpe)/(divider);
                end
            end
        end


        function Jac_fpe = LHS_Jac_f_Fdiff_pedps(obj,pe,ps,ce,cse,T,Ueq,kappa,kappa_D,dx,source)      
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            dp_prop=1.0001;
            len=length(pe);
            Jac_fpe=zeros(len,sol.nb_cell_n+sol.nb_cell_p);


            for i = 1:1:len        
                if i>sol.nb_cell_n && i<=(sol.nb_cell_n+sol.nb_cell_s)
                    continue
                end

                fpe   =  obj.LHS_f_pe_single_cell (pe,ce,kappa,kappa_D,dx,source,i);
                
                for j = max(1,i):1:min(len,i)
                    j_eff=j;
                    if j>sol.nb_cell_n
                        j_eff=j-sol.nb_cell_s;
                    end
 
                    ppdp=ps;
                    if not(ps(j_eff)==0)
                        ppdp(j_eff)=ppdp(j_eff)*dp_prop;
                        divider=(ps(j_eff)*(dp_prop-1));
                    else
                        ppdp(j_eff)=0.000001;
                        divider=0.000001;
                    end

                    if i<=sol.nb_cell_n
                        csmax=p.csn_max;
                    else
                        csmax=p.csp_max;
                    end

                    mock_j=fv.j;
                    mock_j(i)=BV_fun.butler_volmer_singlecell_standalone(i,pe,ppdp,ce,cse,Ueq);

                    %mock_j=BV_fun.butler_volmer_equation(pe,ppdp,ce,cse,T,Ueq,"LHS_Jac_fpe_Fdiff_pedps");
                    mock_source = source;
                    mock_source(i)=obj.update_source_single_cell(i,mock_j,"pe");
                
                    source_used=mock_source;

                    fpe_diff   =  obj.LHS_f_pe_single_cell (pe,ce,kappa,kappa_D,dx,source_used,i);
                    Jac_fpe(i,j_eff)= (fpe_diff-fpe)/(divider);
                end
            end
        end


        function Jac_dfpe_dcs = LHS_Jac_f_Fdiff_pedcs(obj,pe,ps,ce,cse,T,Ueq,kappa,kappa_D,dx,source)      
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            dp_prop=1.0001;
            len=length(pe);
            Jac_dfpe_dcs=zeros(len,(sol.nb_cell_n+sol.nb_cell_p)*(sol.part_nb_cell+1));

            for i = 1:1:len        
                if i>sol.nb_cell_n && i<=(sol.nb_cell_n+sol.nb_cell_s)
                    continue
                end

                fpe   =  obj.LHS_f_pe_single_cell (pe,ce,kappa,kappa_D,dx,source,i);
                
                for j = max(1,i):1:min(len,i)
                    csmax=p.csn_max;
                    j_eff=j;
                    if j>sol.nb_cell_n
                        csmax=p.csp_max;
                        j_eff=j-sol.nb_cell_s;
                    end
 
                    ppdp=cse;
                    if not(cse(j_eff)==0)
                        ppdp(j_eff)=min(ppdp(j_eff)*dp_prop,csmax);
                        divider=(cse(j_eff)*(dp_prop-1));
                    else
                        ppdp(j_eff)=0.000001;
                        divider=0.000001;
                    end

                    mock_j=fv.j;
                    mock_j(i)=BV_fun.butler_volmer_singlecell_standalone(i,pe,ps,ce,ppdp,Ueq);

                    %mock_j=BV_fun.butler_volmer_equation(pe,ps,ce,ppdp,T,Ueq,"LHS_Jac_fpe_Fdiff_pedcs");
                    mock_source = source;
                    mock_source(i)=obj.update_source_single_cell(i,mock_j,"pe");
                
                    source_used=mock_source;

                    fpe_diff   =  obj.LHS_f_pe_single_cell (pe,ce,kappa,kappa_D,dx,source_used,i);
                    Jac_dfpe_dcs(i,j_eff*(sol.part_nb_cell+1))= (fpe_diff-fpe)/(divider);
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

        function fp = LHS_f_ps (obj,p,c,sigma,dx,source,separator_index,source_BC,ID,print_info)
            len=length(p);
            fp=zeros(len,1);
            
            for i = 1:1:len
                fp(i)=obj.LHS_f_ps_single_cell (p,c,sigma,dx,source,i,separator_index,source_BC,0);     
            end
            if print_info==1
                disp("DEBUG BEN eq_build_fncts "+num2str(fp)+"  "+ID)
                disp("DEBUG BEN eq_build_fncts "+num2str(p)+"  "+num2str(transpose(source)))
                disp("DEBUG BEN eq_build_fncts "+"  "+num2str(sigma)+"  "+num2str(source_BC))
            end
        end

        function Jac_fps = LHS_Jac_f_Fdiff_ps(obj,pe,ps,ce,cse,T,Ueq,D,dx,source,separator_index,source_BC)
            dp_prop=1.0001;
            len=length(ps);
            Jac_fps=zeros(len);
            
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            for i = 1:1:len        
                fps   =  obj.LHS_f_ps_single_cell (ps,ce,D,dx,source,i,separator_index,source_BC,0);

                for j = max(1,i-1):1:min(len,i+1)
                    ppdp=ps;
                    if not(ps(j)==0)
                        ppdp(j)=ppdp(j)*dp_prop;
                        divider=(ps(j)*(dp_prop-1));
                    else
                        ppdp(j)=0.000001;
                        divider=0.000001;
                    end

                    if i==j 
                        if i<=sol.nb_cell_n
                            csmax=p.csn_max;
                            indexjacps=i;
                        else
                            csmax=p.csp_max;
                            indexjacps=i+sol.nb_cell_s;
                        end

                        mock_j=fv.j;
                        mock_j(indexjacps)=BV_fun.butler_volmer_singlecell_standalone(indexjacps,pe,ppdp,ce,cse,Ueq);

                        %mock_j=BV_fun.butler_volmer_equation(pe,ppdp,ce,cse,T,Ueq,"LHS_Jac_f_Fdiff_ps");
                        mock_source = source;
                        mock_source(i)=obj.update_source_single_cell(indexjacps,mock_j,"ps");
                    
                        source_used=mock_source;
                    else
                        source_used=source;
                    end

                    fps_diff   =  obj.LHS_f_ps_single_cell (ppdp,ce,D,dx,source_used,i,separator_index,source_BC,0);
                    Jac_fps(i,j)= (fps_diff-fps)/(divider);

                    if i == j & Jac_fps(i,j)==0
                        Jac_fps(i,j)=0.000000001;
                    end
                end 
            end 
        end

        function Jac_fps = LHS_Jac_f_Fdiff_psdpe(obj,pe,ps,ce,cse,T,Ueq,D,dx,source,separator_index,source_BC)
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            dp_prop=1.0001;
            len=length(ps);
            Jac_fps=zeros(len,sol.nb_cell);

            for i = 1:1:len        
                fps   =  obj.LHS_f_ps_single_cell (ps,ce,D,dx,source,i,separator_index,source_BC,0);

                for j = max(1,i):1:min(len,i)
                    if i<=sol.nb_cell_n
                        csmax=p.csn_max;
                        indexjacps=i;
                    else
                        csmax=p.csp_max;
                        indexjacps=i+sol.nb_cell_s;
                    end

                    ppdp=pe;
                    if not(pe(indexjacps)==0)
                        ppdp(indexjacps)=ppdp(indexjacps)*dp_prop;
                        divider=(pe(indexjacps)*(dp_prop-1));
                    else
                        ppdp(indexjacps)=0.000001;
                        divider=0.000001;
                    end

                    if i==j 

                        mock_j=fv.j;
                        mock_j(indexjacps)=BV_fun.butler_volmer_singlecell_standalone(indexjacps,ppdp,ps,ce,cse,Ueq);

                        %mock_j=BV_fun.butler_volmer_equation(ppdp,ps,ce,cse,T,Ueq,"LHS_Jac_f_Fdiff_psdpe");
                        mock_source = source;
                        mock_source(i)=obj.update_source_single_cell(indexjacps,mock_j,"ps");
                    
                        source_used=mock_source;
                    else
                        source_used=source;
                    end

                    fps_diff   =  obj.LHS_f_ps_single_cell (ps,ce,D,dx,source_used,i,separator_index,source_BC,0);
                    Jac_fps(i,j)= (fps_diff-fps)/(divider);

                    if i == j & Jac_fps(i,j)==0
                        Jac_fps(i,j)=0.000000001;
                    end
                end 
            end 
        end

        function Jac_fps = LHS_Jac_f_Fdiff_psdce(obj,pe,ps,ce,cse,T,Ueq,D,dx,source,separator_index,source_BC)
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            dp_prop=1.0001;
            len=length(ps);
            Jac_fps=zeros(len,sol.nb_cell);

            for i = 1:1:len        
                fps   =  obj.LHS_f_ps_single_cell (ps,ce,D,dx,source,i,separator_index,source_BC,0);

                for j = max(1,i):1:min(len,i)
                    if i<=sol.nb_cell_n
                        csmax=p.csn_max;
                        indexjacps=i;
                    else
                        csmax=p.csp_max;
                        indexjacps=i+sol.nb_cell_s;
                    end

                    ppdp=ce;
                    if not(ce(indexjacps)==0)
                        ppdp(indexjacps)=ppdp(indexjacps)*dp_prop;
                        divider=(ce(indexjacps)*(dp_prop-1));
                    else
                        ppdp(indexjacps)=0.000001;
                        divider=0.000001;
                    end

                    if i==j 
                        mock_j=fv.j;
                        mock_j(indexjacps)=BV_fun.butler_volmer_singlecell_standalone(indexjacps,pe,ps,ppdp,cse,Ueq);

                        %mock_j=BV_fun.butler_volmer_equation(pe,ps,ppdp,cse,T,Ueq,"LHS_Jac_f_Fdiff_psdce");
                        mock_source = source;
                        mock_source(i)=obj.update_source_single_cell(indexjacps,mock_j,"ps");
                    
                        source_used=mock_source;
                    else
                        source_used=source;
                    end

                    fps_diff   =  obj.LHS_f_ps_single_cell (ps,ce,D,dx,source_used,i,separator_index,source_BC,0);
                    Jac_fps(i,j)= (fps_diff-fps)/(divider);

                    if i == j & Jac_fps(i,j)==0
                        Jac_fps(i,j)=0.000000001;
                    end
                end 
            end 
        end


        function Jac_dfps_dcs = LHS_Jac_f_Fdiff_psdcs(obj,pe,ps,ce,cse,T,Ueq,D,dx,source,separator_index,source_BC)      
            global fv
            global sol
            global p
            global BV_fun
            mock_source = source;
            mock_j = fv.j;

            dp_prop=1.0001;
            len=length(ps);
            len2=length(pe);
            Jac_dfps_dcs=zeros(len,(sol.nb_cell_n+sol.nb_cell_p)*(sol.part_nb_cell+1));

            for i = 1:1:len2        
                if i>sol.nb_cell_n && i<=(sol.nb_cell_n+sol.nb_cell_s)
                    continue
                end
                
                for j = max(1,i):1:min(len,i)
                    csmax=p.csn_max;
                    j_eff=j;
                    if j>sol.nb_cell_n
                        csmax=p.csp_max;
                        j_eff=j-sol.nb_cell_s;
                    end

                    fps   =  obj.LHS_f_ps_single_cell (ps,ce,D,dx,source,j_eff,separator_index,source_BC,0);

                    ppdp=cse;
                    if not(cse(j_eff)==0)
                        ppdp(j_eff)=min(ppdp(j_eff)*dp_prop,csmax);
                        divider=(cse(j_eff)*(dp_prop-1));
                    else
                        ppdp(j_eff)=0.000001;
                        divider=0.000001;
                    end

                    mock_j=fv.j;
                    mock_j(i)=BV_fun.butler_volmer_singlecell_standalone(i,pe,ps,ce,ppdp,Ueq);

                    %mock_j=BV_fun.butler_volmer_equation(pe,ps,ce,ppdp,T,Ueq,"LHS_Jac_fps_Fdiff_psdcs");
                    mock_source = source;
                    mock_source(j_eff)=obj.update_source_single_cell(i,mock_j,"ps");
                
                    source_used=mock_source;

                    fps_diff   =  obj.LHS_f_ps_single_cell (ps,ce,D,dx,source_used,j_eff,separator_index,source_BC,0);
                    Jac_dfps_dcs(j_eff,j_eff*(sol.part_nb_cell+1))= (fps_diff-fps)/(divider);
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

        function source = update_source_single_cell (obj,i,j,ID) 
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
                source=-p.Faraday*A_s*j(i);
            elseif ID=="pe"
                source= p.Faraday*A_s*j(i);
            elseif ID=="ce"
                source= (1-p.t_plus)*A_s*j(i);
            elseif ID=="cs"
                source= -j(i)/Ds;
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
                source_ps(i)=obj.update_source_single_cell (i,fv.j,"ps") ;
                source_pe(i)=obj.update_source_single_cell (i,fv.j,"pe") ;
                source_ce(i)=obj.update_source_single_cell (i,fv.j,"ce") ;
                source_cs(i)=obj.update_source_single_cell (i,fv.j,"cs") ;
            end

            for i = sol.nb_cell_n+sol.nb_cell_s+1:1:sol.nb_cell
                source_ps(i-sol.nb_cell_s)=obj.update_source_single_cell (i,fv.j,"ps") ;
                source_pe(i)=obj.update_source_single_cell (i,fv.j,"pe") ;
                source_ce(i)=obj.update_source_single_cell (i,fv.j,"ce") ;
                source_cs(i)=obj.update_source_single_cell (i,fv.j,"cs") ;
            end

            source_csn=source_cs(1:sol.nb_cell_n);
            source_csp=source_cs(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell);

        end
    end
end



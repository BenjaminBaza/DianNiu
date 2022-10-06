classdef boundary_condition_functions
    methods

        function [rm,rc,rp,cm,cc,cp]= concentration_solid_BC(obj,c,r,i,source,D,i_rcell)
            len=length(c);
            lenr=length(r);
            if i_rcell==lenr
                rm=r(i_rcell-1);    rc=r(i_rcell);  rp=r(lenr)+ r(lenr) -r(lenr-1);
                cm=c(i-1);          cc=c(i);        cp=cc + source*(r(lenr)-r(lenr-1));
                %disp("DEBUG BEN solid BC "+num2str(cc)+"  "+num2str(cp)+"  "+num2str(source)+"  "+num2str(r(lenr)-r(lenr-1)))
            elseif i_rcell==1
                rm=-r(i_rcell+1);   rc=r(i_rcell);  rp=r(i_rcell+1);
                cm=c(i+1);          cc=c(i);        cp=c(i+1);
            else
                rm=r(i_rcell-1);    rc=r(i_rcell);  rp=r(i_rcell+1);
                cm=c(i-1);          cc=c(i);        cp=c(i+1);
            end
        end

        function [dxm,dxc,dxp,pm,pc,pp,cm,cc,cp,Dm,Dc,Dp]= concentration_electrolyte_BC(obj,p,c,dx,i,D)
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

        function [dxm,dxc,dxp,pm,pc,pp,cm,cc,cp,Dm,Dc,Dp,kDm,kDc,kDp]= Potential_electrolyte_BC(obj,p,c,dx,i,D,kD)
            global deb
            len=length(c);
            if i==len
                dxm=dx(i-1);dxc=dx(i);dxp=dx(i);
                cm=c(i-1);cc=c(i);cp=c(i);
                pm=p(i-1);pc=p(i);pp=p(i);
                Dm=D(i-1);Dc=D(i);Dp=D(i);
                kDm=kD(i-1);kDc=kD(i);kDp=kD(i);
            elseif i==1
                dxm=dx(i);dxc=dx(i);dxp=dx(i+1);
                cm=c(i);cc=c(i);cp=c(i+1);
                if deb.safe_BC_mode==1
                    pm=-p(i);pc=p(i);pp=p(i+1); % DEBUG BEN This BC forces potential to 0 at x=0. It is incorrect but completes closure of the system. This is temporary
                else
                    pm=p(i);pc=p(i);pp=p(i+1);  %the BC at the edge of the cell is zero gradient
                end
                Dm=D(i);Dc=D(i);Dp=D(i+1);
                kDm=kD(i);kDc=kD(i);kDp=kD(i+1);
            else
                dxm=dx(i-1);dxc=dx(i);dxp=dx(i+1);
                cm=c(i-1);cc=c(i);cp=c(i+1);
                pm=p(i-1);pc=p(i);pp=p(i+1);
                Dm=D(i-1);Dc=D(i);Dp=D(i+1);
                kDm=kD(i-1);kDc=kD(i);kDp=kD(i+1);
            end
        end

        function [dxm,dxc,dxp,pm,pc,pp,Dm,Dc,Dp]= Potential_solid_BC(obj,p,dx,i,D,separator_index,source_BC)
            global deb
            len=length(p);
            if i==len
                dxm=dx(i-1);dxc=dx(i);dxp=dx(i);
                pm=p(i-1);pc=p(i);pp=p(i)+source_BC(2)*dxc;
                Dm=D(i-1);Dc=D(i);Dp=D(i);
            elseif i==1
                dxm=dx(i);dxc=dx(i);dxp=dx(i+1);
                %pm=p(i)-source_BC(1)*dxc;pc=p(i);pp=p(i+1); % This boundary is correct but apparently redundant. 
                pm=0;pc=p(i);pp=p(i+1); % This boundary is correct but apparently redundant. 
                %pm=-p(i);pc=p(i);pp=p(i+1);   % DEBUG BEN This BC forces potential to 0 at x=0. Which appears to be correct (see Plett and Qiao).
                Dm=D(i);Dc=D(i);Dp=D(i+1);
            elseif i==separator_index
                dxm=dx(i-1);dxc=dx(i);dxp=dx(i);
                pm=p(i-1);pc=p(i);pp=p(i);
                Dm=D(i-1);Dc=D(i);Dp=D(i);
            elseif i==separator_index+1
                dxm=dx(i);dxc=dx(i);dxp=dx(i+1);
                if deb.safe_BC_mode==1
                    pm=p(i-1);pc=p(i);pp=p(i+1); %DEBUG BEN This BC connects potential in the positive electrode to the negative electrode
                else
                    pm=p(i);pc=p(i);pp=p(i+1); %the BC at the edge of the electrode is zero gradient
                end
                Dm=D(i);Dc=D(i);Dp=D(i+1);
            else
                dxm=dx(i-1);dxc=dx(i);dxp=dx(i+1);
                pm=p(i-1);pc=p(i);pp=p(i+1);
                Dm=D(i-1);Dc=D(i);Dp=D(i+1);
            end
        end

    end
end
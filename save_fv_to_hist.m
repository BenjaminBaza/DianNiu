function save_fv_to_hist(ite)
    global hist
    global sol
    global fv
    
    hist.cse(:,ite+1) = fv.cse;
    hist.csn(:,ite+1) = reshape(fv.csn,sol.nb_cell_n*(sol.part_nb_cell+1),1);
    hist.csp(:,ite+1) = reshape(fv.csp,sol.nb_cell_p*(sol.part_nb_cell+1),1);
    hist.ce(:,ite+1)  = fv.ce;
    hist.pe(:,ite+1)  = fv.pe;
    hist.ps(:,ite+1)  = fv.ps;
    hist.Ueq(:,ite+1) = fv.Ueq;
    hist.j(:,ite+1)   = fv.j;
    hist.V(ite)       = fv.V;
    hist.SOC_neg(ite)       = fv.SOC_neg;
    hist.SOC_pos(ite)       = fv.SOC_pos;

    hist.maxjn=max(max(fv.j(1:sol.nb_cell_n)),hist.maxjn);
    hist.minjn=min(min(fv.j(1:sol.nb_cell_n)),hist.minjn);
    hist.maxjp=max(max(fv.j(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell)),hist.maxjp);
    hist.minjp=min(min(fv.j(sol.nb_cell_n+sol.nb_cell_s+1:sol.nb_cell)),hist.minjp);

    hist.maxdt=max(sol.dt,hist.maxdt);
    hist.mindt=min(sol.dt,hist.mindt);

    hist.maxcsn=max(max(max(fv.csn)),hist.maxcsn);
    hist.mincsn=min(min(min(fv.csn)),hist.mincsn);
    hist.maxcsp=max(max(max(fv.csp)),hist.maxcsp);
    hist.mincsp=min(min(min(fv.csp)),hist.mincsp);


    dcsn_R=fv.csn(length(sol.part_coord_n),:)-fv.csn(1,:);
    dcsp_R=fv.csp(length(sol.part_coord_p),:)-fv.csp(1,:);
    hist.maxdcsn=max(max(dcsn_R),hist.maxdcsn);
    hist.mindcsn=min(min(dcsn_R),hist.mindcsn);
    hist.maxdcsp=max(max(dcsp_R),hist.maxdcsp);
    hist.mindcsp=min(min(dcsp_R),hist.mindcsp);

end
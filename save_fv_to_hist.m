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


end
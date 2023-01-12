function save_data()
	global p
	global fv
	global hist
	global ex
	global sol
    global deb


    disp('Writing data')	

    mkdir(deb.folder_name) 
    copyfile(deb.read_ctrl_file_name, deb.folder_name);

    %writestruct(ex , deb.folder_name+'/ex.xls');
    %writestruct(fv , deb.folder_name+'/fv.xls');
    %writestruct(hist , deb.folder_name+'/hist.xls');
    %writestruct(p , deb.folder_name+'/p.xls');
    %writestruct(sol , deb.folder_name+'/sol.xls');


    value =[deb.case_name, sol.nb_cell_n, sol.nb_cell_s, sol.nb_cell_p, sol.part_nb_cell, fv.V, fv.SOC_neg, fv.SOC_pos, hist.charge_time, sol.time_tot];
    value2=[deb.case_name, sol.nb_cell_n, sol.nb_cell_s, sol.nb_cell_p, sol.part_nb_cell, hist.maxjn, hist.minjn, hist.maxjp, hist.minjp, hist.maxdt, hist.mindt, ...
            hist.maxcsn, hist.mincsn, hist.maxcsp, hist.mincsp, hist.maxdcsn, hist.mindcsn, hist.maxdcsp, hist.mindcsp];
    [num,txt,raw] = xlsread('../Saved_data_DFN_DianNiu/Data_recap.xlsx');
    rawsize=size(raw);
    xlRange = "A"+num2str(rawsize(1)+1);
    xlswrite('../Saved_data_DFN_DianNiu/Data_recap.xlsx',value,'Sheet1',xlRange);
    xlswrite('../Saved_data_DFN_DianNiu/Data_recap.xlsx',value2,'Sheet2',xlRange);

    if deb.write_output_data>1

        %Clean data
        % Remove zero columns
        hist.csn( :, ~any(hist.csn,1) ) = [];  %columns
        hist.csp( :, ~any(hist.csp,1) ) = [];  %columns
        hist.cse( :, ~any(hist.cse,1) ) = [];  %columns
        hist.ce( :, ~any(hist.ce,1) ) = [];  %columns
        hist.pe(:,1)=hist.pe(:,1)-0.0000001;
        hist.pe( :, ~any(hist.pe,1) ) = [];  %columns
        hist.ps( :, ~any(hist.ps,1) ) = [];  %columns
        hist.j( :, ~any(hist.j,1) ) = [];  %columns
        hist.V( :, ~any(hist.V,1) ) = [];  %columns
        hist.Ueq( :, ~any(hist.Ueq,1) ) = [];  %columns
        hist.SOC_neg( :, ~any(hist.SOC_neg,1) ) = [];  %columns
        hist.SOC_pos( :, ~any(hist.SOC_pos,1) ) = [];  %columns
        hist.residuals( :, ~any(hist.residuals,1) ) = [];  %columns
        hist.residuals_time( :, ~any(hist.residuals_time,1) ) = [];  %columns
        hist.residuals_diff( :, ~any(hist.residuals_diff,1) ) = [];  %columns
        hist.newt_it_number( :, ~any(hist.newt_it_number,1) ) = [];  %columns
        hist.delta_coupled( :, ~any(hist.delta_coupled,1) ) = [];  %columns

        disp('Writing csn')
        writematrix(transpose(hist.csn), deb.folder_name+'/hist_csn.csv')  
        disp('Writing csp')
        writematrix(transpose(hist.csp), deb.folder_name+'/hist_csp.csv')  
        disp('Writing cse')
        writematrix(transpose(hist.cse), deb.folder_name+'/hist_cse.csv')  
        disp('Writing ce')
        writematrix(transpose(hist.ce), deb.folder_name+'/hist_ce.csv')  
        disp('Writing pe')
        writematrix(transpose(hist.pe), deb.folder_name+'/hist_pe.csv')  
        disp('Writing ps')
        writematrix(transpose(hist.ps), deb.folder_name+'/hist_ps.csv')  
        disp('Writing j')
        writematrix(transpose(hist.j), deb.folder_name+'/hist_j.csv')  
        disp('Writing V')
        writematrix(transpose(hist.V), deb.folder_name+'/hist_V.csv')  
        disp('Writing Ueq')
        writematrix(transpose(hist.Ueq), deb.folder_name+'/hist_Ueq.csv')  
        disp('Writing SOCn')
        writematrix(transpose(hist.SOC_neg), deb.folder_name+'/hist_SOC_neg.csv')  
        disp('Writing SOCp')
        writematrix(transpose(hist.SOC_pos), deb.folder_name+'/hist_SOC_pos.csv')  
        disp('Writing residuals')
        writematrix(transpose(hist.residuals), deb.folder_name+'/hist_residuals.csv')  
        disp('Writing residuals time')
        writematrix(transpose(hist.residuals_time), deb.folder_name+'/hist_residuals_time.csv')  
        disp('Writing residuals diff')
        writematrix(transpose(hist.residuals_diff), deb.folder_name+'/hist_residuals_diff.csv')  
        disp('Writing Newt ite')
        writematrix(transpose(hist.newt_it_number), deb.folder_name+'/hist_newt_it_number.csv')  
        disp('Writing delta coupled')
        writematrix(transpose(hist.delta_coupled), deb.folder_name+'/hist_delta_coupled.csv')  
    end
    disp('Done writing data')    

end
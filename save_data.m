function save_data()
	global p
	global fv
	global hist
	global ex
	global sol
    global deb


    disp('Writing data')	

    mkdir(deb.folder_name) 
    %writestruct(ex , deb.folder_name+'/ex.xls');
    %writestruct(fv , deb.folder_name+'/fv.xls');
    %writestruct(hist , deb.folder_name+'/hist.xls');
    %writestruct(p , deb.folder_name+'/p.xls');
    %writestruct(sol , deb.folder_name+'/sol.xls');

    writematrix(transpose(hist.csn), deb.folder_name+'/hist_csn.xlsx')  
    writematrix(transpose(hist.csp), deb.folder_name+'/hist_csp.xlsx')  
    writematrix(transpose(hist.cse), deb.folder_name+'/hist_cse.xlsx')  
    writematrix(transpose(hist.ce), deb.folder_name+'/hist_ce.xlsx')  
    writematrix(transpose(hist.pe), deb.folder_name+'/hist_ps.xlsx')  
    writematrix(transpose(hist.ps), deb.folder_name+'/hist_ps.xlsx')  
    writematrix(transpose(hist.j), deb.folder_name+'/hist_j.xlsx')  
    writematrix(transpose(hist.V), deb.folder_name+'/hist_V.xlsx')  
    writematrix(transpose(hist.Ueq), deb.folder_name+'/hist_Ueq.xlsx')  
    writematrix(transpose(hist.SOC_neg), deb.folder_name+'/hist_SOC_neg.xlsx')  
    writematrix(transpose(hist.SOC_pos), deb.folder_name+'/hist_SOC_pos.xlsx')  
    writematrix(transpose(hist.residuals), deb.folder_name+'/hist_residuals.xlsx')  
    writematrix(transpose(hist.residuals_time), deb.folder_name+'/hist_residuals_time.xlsx')  
    writematrix(transpose(hist.residuals_diff), deb.folder_name+'/hist_residuals_diff.xlsx')  
    writematrix(transpose(hist.newt_it_number), deb.folder_name+'/hist_newt_it_number.xlsx')  
    writematrix(transpose(hist.delta_coupled), deb.folder_name+'/hist_delta_coupled.xlsx')  

    disp('Done writing data')    

end
% Voltage calculation

function [V,SOC_neg,SOC_pos] = voltage_calc(ps, R,I)
	global deb
	global p
	global fv

	len=length(ps);
	V=ps(len)-ps(1) -R*I;

	SOC_neg=((mean(mean(fv.csn))/p.csn_max)-p.neg_stoichiometry_min)/(p.neg_stoichiometry_max-p.neg_stoichiometry_min);
	SOC_pos=((mean(mean(fv.csp))/p.csp_max)-p.pos_stoichiometry_max)/(p.pos_stoichiometry_min-p.pos_stoichiometry_max);

	if deb.prints>=5
		disp("DEBUG BEN volt calc")
		disp(transpose(ps)) 
		disp(num2str(R)+ "  " + num2str(I)+ "  " + num2str(R*I)) 
	end
end
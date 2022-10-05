% Voltage calculation

function V= voltage_calc(ps, R,I)
	global deb

	len=length(ps);
	V=ps(len)-ps(1) -R*I;
	if deb.prints>=2
		disp("DEBUG BEN volt calc")
		disp(transpose(ps)) 
		disp(num2str(R)+ "  " + num2str(I)+ "  " + num2str(R*I)) 
	end
end
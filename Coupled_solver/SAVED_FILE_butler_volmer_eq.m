% Butler Volmer equation

function j= butler_volmer_eq(pe,ps,ce,cse,k0,alpha,F,R,T,Ueq,Rfilm,csnmax,cspmax,neg_nb,sep_nb)
	global deb
	nu=1;
	maxnu=10^6;
	len=length(ce);
	j=zeros(1,len);

	for i = 1:1:len
		if i<=neg_nb
			cseloc=cse(i);
			indloc=i;
			i0=F*k0(1)*ce(i)^(1-alpha)*(csnmax-cseloc)^(1-alpha)*cseloc^alpha;
			nu=ps(i)-pe(i)-Ueq(i); %Could add the ionic resistance of the film layer (+F*Rfilm*j(i)) and iterate to find the combined vaues of nu and j
			nu=max(-maxnu,min(nu,maxnu));
			j(i)=i0/F* (exp((1-alpha)*F*nu/(R*T)) - exp((-alpha)*F*nu/(R*T)));
		elseif i<=neg_nb+sep_nb
			j(i)=0;
        else
        	cseloc=cse(i-sep_nb);
        	indloc=i-sep_nb;
			i0=F*k0(2)*ce(i)^(1-alpha)*(cspmax-cseloc)^(1-alpha)*cseloc^alpha;
			nu=ps(indloc)-pe(i)-Ueq(i); %Could add the ionic resistance of the film layer (+F*Rfilm*j(i)) and iterate to find the combined vaues of nu and j
			nu=max(-maxnu,min(nu,maxnu));
			j(i)=i0/F* (exp((1-alpha)*F*nu/(R*T)) - exp((-alpha)*F*nu/(R*T)));
			
			if deb.prints>=2 & isreal(j(i))==0
				disp("DEBUG BEN BV equation")
				disp(num2str(i0)+"   "+num2str(nu)+"   "+num2str(j(i))+"   "+num2str(- exp((-alpha)*F*nu/(R*T)))+"   "+num2str(exp((1-alpha)*F*nu/(R*T))));
				disp(num2str(ps(i-sep_nb))+"   "+num2str(-pe(i))+"   "+num2str(-Ueq(i)));
				disp(num2str(ce(i))+"   "+num2str(cspmax-cse(i-sep_nb))+"   "+num2str(cse(i-sep_nb)));
			end
		end


		if i>neg_nb+sep_nb || i<=neg_nb
			if abs(j(i))>10^7 || isnan(j(i))==1 || deb.prints>=3
				disp("DEBUG BEN BV eq")
				disp(num2str(i)+"  "+num2str(nu)+"  "+num2str(j(i))+"  "+num2str((exp((1-alpha)*F*nu/(R*T)) - exp((-alpha)*F*nu/(R*T)))) ...
					+"  "+num2str((exp((1-alpha)*F*nu/(R*T)) ))+"  "+num2str((-exp((-alpha)*F*nu/(R*T))))+"   "+ ...
					num2str((exp((1-alpha)*F*(3*10^1)/(R*T)) - exp((-alpha)*F*(3*10^1)/(R*T))))+"   "+ ...
					num2str((exp((1-alpha)*F*(-3*10^1)/(R*T)) - exp((-alpha)*F*(-3*10^1)/(R*T)))))

				disp(num2str((1-alpha)*F/(R*T))+"  "+num2str((R*T))+"  "+num2str((1-alpha)*nu*F/(R*T)))

				disp(num2str(i)+"  "+num2str(i0)+"   "+num2str(ce(i)^(1-alpha))+"   "+num2str(ce(i))+"   "+ num2str((cspmax-cseloc)^(1-alpha)*cseloc^alpha)+"  "+num2str(k0)+"  "+num2str(F))

				disp(num2str(i)+"  "+num2str(nu)+"  "+num2str(Ueq(i))+"  "+num2str(ps(indloc))+"  "+num2str(pe(i)))

				disp(transpose(Ueq))
				disp(transpose(ps))
				disp(transpose(pe))

				if abs(j(i))>10^10
					j(i)=10^10;
				end
			end
		end
		
	end

end
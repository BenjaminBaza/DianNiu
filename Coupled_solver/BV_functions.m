% Butler Volmer equation
classdef BV_functions
    methods

		function [j,nu,i0]= butler_volmer_singlecell(obj,pe,ps,ce,cse,k0,alpha,F,R,T,Ueq,Rfilm,csmax)
			global deb
			nu=1;
			maxnu=10^6;
			len=length(ce);
			j=zeros(1,len);

			cseloc=cse;
			%cseloc=min(cse, csmax);
			if cse>csmax
				disp("DEBUG BEN BV eq cse>csmax "+num2str(cse)+"  "+num2str(csmax))
			end
			
			i0=F*k0*ce^(1-alpha)*(csmax-cseloc)^(1-alpha)*cseloc^alpha;
			nu=ps-pe-Ueq; %Could add the ionic resistance of the film layer (+F*Rfilm*j(i)) and iterate to find the combined vaues of nu and j
			nu=max(-maxnu,min(nu,maxnu));
			j=i0/F* (exp((1-alpha)*F*nu/(R*T)) - exp((-alpha)*F*nu/(R*T)));
		end
	


		function j= butler_volmer_equation(obj,pe,ps,ce,cse,T,Ueq,ID)
			global deb
			global p
			global sol
			nu=1;
			maxnu=10^6;
			len=length(ce);
			j=zeros(1,len);

			for i = 1:1:len
				if i<=sol.nb_cell_n
					cseloc=cse(i);
					indloc=i;

					[j(i),nu,i0]=obj.butler_volmer_singlecell(pe(i),ps(indloc),ce(i),cseloc,p.k0(1),p.alpha,p.Faraday,p.Rg,T,Ueq(i),p.Rfilm,p.csn_max);

				elseif i<=sol.nb_cell_n+sol.nb_cell_s
					j(i)=0;
		        else
		        	cseloc=cse(i-sol.nb_cell_s);
		        	indloc=i-sol.nb_cell_s;

		        	[j(i),nu,i0]=obj.butler_volmer_singlecell(pe(i),ps(indloc),ce(i),cseloc,p.k0(2),p.alpha,p.Faraday,p.Rg,T,Ueq(i),p.Rfilm,p.csp_max);

					if deb.prints>=2 & isreal(j(i))==0
						disp("DEBUG BEN BV equation cell "+num2str(i)+" ID="+ID)
						disp(num2str(i0)+"   "+num2str(nu)+"   "+num2str(j(i))+"   "+num2str(- exp((-p.alpha)*p.Faraday*nu/(p.Rg*T)))+"   "+num2str(exp((1-p.alpha)*p.Faraday*nu/(p.Rg*T))));
						disp(num2str(ps(i-sol.nb_cell_s))+"   "+num2str(-pe(i))+"   "+num2str(-Ueq(i)));
						disp(num2str(ce(i))+"   "+num2str(p.csp_max-cse(i-sol.nb_cell_s))+"   "+num2str(cse(i-sol.nb_cell_s)));
					end
				end

				if i>sol.nb_cell_n+sol.nb_cell_s || i<=sol.nb_cell_n
					if ((abs(j(i))>10^10 || isnan(j(i))==1 || deb.prints>=5 )&& not(contains(ID,"LHS_Jac_")))
						disp("DEBUG BEN BV eq ID="+ID)
						disp(num2str(i)+"  "+num2str(nu)+"  "+num2str(j(i))+"  "+num2str((exp((1-p.alpha)*p.Faraday*nu/(p.Rg*T)) - exp((-p.alpha)*p.Faraday*nu/(p.Rg*T)))) ...
							+"  "+num2str((exp((1-p.alpha)*p.Faraday*nu/(p.Rg*T)) ))+"  "+num2str((-exp((-p.alpha)*p.Faraday*nu/(p.Rg*T))))+"   "+ ...
							num2str((exp((1-p.alpha)*p.Faraday*(3*10^1)/(p.Rg*T)) - exp((-p.alpha)*p.Faraday*(3*10^1)/(p.Rg*T))))+"   "+ ...
							num2str((exp((1-p.alpha)*p.Faraday*(-3*10^1)/(p.Rg*T)) - exp((-p.alpha)*p.Faraday*(-3*10^1)/(p.Rg*T)))))

						disp(num2str((1-p.alpha)*p.Faraday/(p.Rg*T))+"  "+num2str((p.Rg*T))+"  "+num2str((1-p.alpha)*nu*p.Faraday/(p.Rg*T)))

						disp(num2str(i)+"  "+num2str(i0)+"   "+num2str(ce(i)^(1-p.alpha))+"   "+num2str(ce(i))+"   "+ num2str((p.csp_max-cseloc)^(1-p.alpha)*cseloc^p.alpha)+"  "+num2str(p.k0)+"  "+num2str(p.Faraday))

						disp(num2str(i)+"  "+num2str(nu)+"  "+num2str(Ueq(i))+"  "+num2str(ps(indloc))+"  "+num2str(pe(i)))

						disp(transpose(Ueq))
						disp(transpose(pe))
						disp(transpose(ps))
						disp(transpose(ce))
						disp((cse))

						if abs(j(i))>10^10
							j(i)=10^10;
						end
					end
				end
			end
		end
	end

end
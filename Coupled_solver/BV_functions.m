% Butler Volmer equation
classdef BV_functions
    methods

    	function [j,nu,i0]= butler_volmer_singlecell_standalone(obj,i,pe,ps,ce,cse,Ueq)
			global eq_build_fun
			global deb
			global sol
			global ini
			global p
			
			if sol.nb_cell_n< i && i<=sol.nb_cell_n+sol.nb_cell_s
				j=0;
				nu=0;
				i0=0;
				return
			end

			maxnu=10^6;
			k0_loc=1;
			csmax=p.csn_max;
			if i>sol.nb_cell_n
				k0_loc=2;
				csmax=p.csp_max;
			end

			ps_ind=eq_build_fun.pe2ps_index(i);
			
			if ce(i)<0
				ce_term= - (-ce(i))^(1-p.alpha);
			else
				ce_term= ce(i)^(1-p.alpha);
			end

			if cse(ps_ind)<0
				i0= p.Faraday*p.k0(k0_loc)* ce_term * (csmax-cse(ps_ind))^(1-p.alpha)		* (-(-cse(ps_ind))^p.alpha);				
			elseif cse(ps_ind)>csmax
				i0= p.Faraday*p.k0(k0_loc)* ce_term * (-(-(csmax-cse(ps_ind)))^(1-p.alpha))	* cse(ps_ind)^p.alpha;
			else
				i0= p.Faraday*p.k0(k0_loc)* ce_term * (csmax-cse(ps_ind))^(1-p.alpha)		* cse(ps_ind)^p.alpha;
			end

			nu=ps(ps_ind)-pe(i)-Ueq(i); %Could add the ionic resistance of the film layer (+F*Rfilm*j(i)) and iterate to find the combined vaues of nu and j
			nu=max(-maxnu,min(nu,maxnu));
			j=i0/p.Faraday* (exp((1-p.alpha)*p.Faraday*nu/(p.Rg*ini.T0)) - exp((-p.alpha)*p.Faraday*nu/(p.Rg*ini.T0)));
		end


		function [j,nu,i0]= butler_volmer_singlecell(obj,pe,ps,ce,cse,k0,alpha,F,R,T,Ueq,Rfilm,csmax)
			global eq_build_fun
			global deb
			
			maxnu=10^6;
			
			if ce<0
				ce_term= - (-ce)^(1-alpha);
			else
				ce_term= ce^(1-alpha);
			end

			if cse<0
				i0= F*k0* ce_term * (csmax-cse)^(1-alpha)		* (-(-cse)^alpha);				
			elseif cse>csmax
				i0= F*k0* ce_term * (-(-(csmax-cse))^(1-alpha))	* cse^alpha;
			else
				i0= F*k0* ce_term * (csmax-cse)^(1-alpha)		* cse^alpha;
			end

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

				if (i>sol.nb_cell_n+sol.nb_cell_s || i<=sol.nb_cell_n)
					if ((abs(j(i))>10^10 || isnan(j(i))==1 || deb.prints>=5 )&& not(contains(ID,"LHS_Jac_")))
						if deb.prints>=1
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
						end
						if abs(j(i))>10^10
							j(i)=10^10;
						end
					end
				end
			end
		end
	end

end
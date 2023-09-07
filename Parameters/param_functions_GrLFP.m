classdef param_functions_GrLFP
    %from Dandeliion library and Chen et al. 2020
    methods

        function kappa_e = electrolyte_conductivity(obj,c)
            x=c*0.001;
            kappa_e = 0.2667*x^3 - 1.2983*x^2 + 1.7919*x +0.1726 ;
        end

        function De = electrolyte_diffusivity(obj,c)
            x=c*0.001;
            De = (2.646*10^(-6)*(0.2667*x^3 - 1.2983*x^2 + 1.7919*x +0.1726)/(max(x,0.00000001))) * 0.0001 ;
        end

        function Dcs = neg_electrode_diffusivity(obj,c)
            global p
            x=c/p.csn_max;
            Dcs = (8.4e-9 * exp(-11.3 * x) + 8.2e-11)*0.0001 ;
        end

        function Dcs = pos_electrode_diffusivity(obj,c)
            global p
            x=c/p.csp_max;
            Dcs = (9.e-14)*0.0001;
        end

        function Ueqn = neg_electrode_Ueq(obj,c,i)
            global p
            x=c/p.csn_max;
            x=min(max(0,x),1);
            Ueqn = 0.7165 * exp(-369.03 * x) + 0.1219 * exp(-35.648 * (x - 0.053095)) - 0.018919 * tanh(21.197 * (x - 0.19618)) ...
                 - 0.016964 * tanh(27.137 * (x - 0.31283)) - 0.019931 * tanh(28.57 * (x - 0.61422)) - 0.93115 * exp(36.328 * (x - 1.1074)) + 0.14003 ;     
        end

        function Ueqp = pos_electrode_Ueq(obj,c,i)
            global p
            x=c/p.csp_max;
            x=min(max(0,x),1);
            Ueqp = (3.114559 + 4.438792 * atan(-71.7352 * x + 70.85337) - 4.240252 * atan(-68.5605 * x + 67.730082));
        end
    end
end

classdef param_functions_GrNMC
    %from Dandeliion library and Chen et al. 2020
    methods

        function kappa_e = electrolyte_conductivity(obj,c)
            x=c*0.001;
            kappa_e = 0.2667*x^3 - 1.2983*x^2 + 1.7919*x +0.1726 ;
        end

        function De = electrolyte_diffusivity(obj,c)
            x=c*0.001;
            De = (2.646*10^(-10)*(0.2667*x^3 - 1.2983*x^2 + 1.7919*x +0.1726)/(max(x,0.00000001)))*0.0001 ;
        end

        function Dcs = neg_electrode_diffusivity(obj,c)
            global p
            x=c/p.csn_max;
            Dcs = (8.4e-9 * exp(-11.3 * x) + 8.2e-11)*0.0001 ;
        end

        function Dcs = pos_electrode_diffusivity(obj,c)
            global p
            x=c/p.csp_max;
            Dcs = (3.7e-9 - 3.4e-9 * exp(-12 * (x - 0.62) * (x - 0.62)))*0.0001;
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
            Ueqp = -2.3521 * x - 0.074706 * tanh(31.886 * (x - 0.021992)) + 6.3498 * tanh(2.664 * (x - 0.17435)) ...
                 - 0.64024 * tanh(5.4862 * (x - 0.43925)) - 3.8238 * tanh(4.1217 * (x - 0.17619)) - 0.054212 * tanh(18.292 * (x - 0.76227)) + 4.2329;
        end
    end
end
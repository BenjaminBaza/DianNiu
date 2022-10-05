classdef param_functions_LGM50
    %from Dandeliion library and Chen et al. 2020
    methods

        function kappa_e = electrolyte_conductivity(obj,c)
            x=c*0.001;
            kappa_e = 0.1297*x^3 - 2.51*x^1.5 + 3.329*x +0.00000001 ;
        end

        function De = electrolyte_diffusivity(obj,c)
            x=c*0.001;
            De = 8.794*10^(-7)*x^2 - 3.972*10^(-6)*x + 4.862*10^(-6) ;
        end

        function Ueqn = neg_electrode_Ueq(obj,c,i)
            global p
            x=c/p.csn_max;
            Ueqn = 1.9793 * exp(-39.3631 * x) + 0.2482 - 0.0909 * tanh(29.8538 * (x - 0.1234)) - 0.04478 * tanh(14.9159 * (x - 0.2769)) - 0.0205 * tanh(30.4444 * (x - 0.6103));
        end

        function Ueqp = pos_electrode_Ueq(obj,c,i)
            global p
            x=c/p.csp_max;
            Ueqp = -0.8090 * x + 4.4875 - 0.0428 * tanh(18.5138 * (x - 0.5542)) - 17.7326 * tanh(15.7890 * (x - 0.3117)) + 17.5842 * tanh(15.9308 * (x - 0.3120));
        end
    end
end
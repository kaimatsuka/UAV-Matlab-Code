function [ true_false ] = check_sderivs( SDERIV )
%% calc_stability_derivatives.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:
%   This function checks the stability derivatives signs and makes sure
%   they are the correct signs
%   
% INPUTS:
%   SDERIV: stability derivative structure
% 
% OUTPUTS:
%   true_false: number indicating if passed all stability derivatives

if(SDERIV.CL_a    >  0 & ...
   SDERIV.CL_adot >= 1 & ...
   SDERIV.CL_q    >= 3 & ...
   SDERIV.CD_a    >  0 & ...
   SDERIV.Cm_a    <  0 & ...
   SDERIV.Cm_q    <  -5 & ...
   SDERIV.Cm_adot <= -3 & ...
   SDERIV.CY_beta <=  -0.125 & ...
   SDERIV.CY_r    >= 0 & ...
   SDERIV.Cl_beta <= 0 & ...   %look into this
   SDERIV.Cl_p    <  -0.1 & ...
   SDERIV.Cl_r    >= 0.05 & ...
   SDERIV.Cn_beta >=  0.025 & ...
   SDERIV.Cn_p    <  -0.025 )
%   SDERIV.CL_de   >= 0 &&
%   SDERIV.Cl_da   <= 0 && 
%   SDERIV.CY_dr   >= 0 && 
%   SDERIV.Cm_de   <= 0 &&
%   SDERIV.CD_de   >  0 && 
%   SDERIV.Cn_da   >= 0 && 
% SDERIV.Cn_r    <= -0.1
%   SDERIV.Cn_dr   <= 0 )
        true_false = 1;
else
        true_false = 0;
end

end


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
   SDERIV.CL_adot >  1 & ...
   SDERIV.CL_q    >= 0 & ...
   SDERIV.CL_de   >= 0 & ...
   SDERIV.Cm_a    <  0 & ...
   SDERIV.Cm_q    <  0 & ...
   SDERIV.Cm_adot <= 0 & ...
   SDERIV.Cm_de   <= -0.5 & ...
   SDERIV.CY_beta <= -0.225 & ... %original -0.25
   SDERIV.CY_r    >= 0 & ...
   SDERIV.CY_dr   >= 0.07 & ... %original 0.05
   SDERIV.Cl_a    >= 0 & ...
   SDERIV.Cl_beta <= 0 & ...   %look into this
   SDERIV.Cl_p    <  0 & ... doesn't fail much
   SDERIV.Cl_r    >= 0 & ...
   SDERIV.Cl_dr   >= 0.0001 &...
   SDERIV.Cn_beta >= 0.01 & ...
   SDERIV.Cn_p    <  0    & ...
   SDERIV.Cn_da   >= 0.001 & ...
   SDERIV.Cn_r    <= 0 & ... %original -0.03 
   SDERIV.Cn_dr   <= -0.0005)
%{
      
   
%}
   
%   SDERIV.Cl_da   <= 0 && 
%   SDERIV.CD_a    >  0 & ... % ignore this 
%   SDERIV.CD_de   >= 0 & ...

        true_false = 1;
else
        true_false = 0;
end

end
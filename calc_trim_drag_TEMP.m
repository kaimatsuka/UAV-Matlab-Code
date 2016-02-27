% calc_trim_drag.TEMP.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%   This function calculates the the elevator deflection given CL value.
%
% INPUTS:
%      SDERIV.CL_a
%
%
% OUTPUTS:
%         de  % elevator deflection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%i_w = 0; % incident angle of the wing

%i_t = -(SDERIV.CM_ac*SDERIV.CL_a + SDERIV.CM_a*CL)/...
        %(SDERIV.CL_a*Cm_i-SDERIV.Cm_a*SDERIV.CL_i); % trim incident angle

%SDERIV.CL_h = CL_a((alpha + i_w)*(1-eps_a)+ ...
              %(i_t-i_w)-alpha_0);
          
%DRAG.CDtrim = K*(SDERIV.CL_a*(alpha + i_w))^2 + ...
              %Eta_h*htail.S/wing.S*K_h*(SDERIV.CL_h)^2; % trim drag


              
de = (SDERIV.Cm0*SDERIV.CL_a + SDERIV.Cm_a*CL)/...
     (SDERIV.CL_a*SDERIV.Cm_de - SDERIV.Cm_a*SDERIV.CL_de);


 


function y = C_f(Re,M)

%skin friction Raymer Equation

y = 0.455./(log10(Re).^2.58.*(1+0.144*M.^2).^0.65);

end
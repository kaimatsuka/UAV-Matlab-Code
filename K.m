function y = K(x,c,t,M,swp_ang)

%Raymer Form Factor Equation

y = (1+0.6/(x/c)*(t/c)+100*(t/c)^4)*...
    (1.34*M.^0.18*cos(swp_ang)^0.28);

end
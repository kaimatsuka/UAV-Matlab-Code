function y = W_FS(F_G,int,N_T, N_E) 

%Nicolai fuel system weight estimate

y = 2.49*(((F_G^0.6)*((1/(1+int))^0.3)*(N_T^0.2)*(N_E^0.13))^1.21);

end
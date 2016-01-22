function y = W_vt(W_TO,N,S_V,b_V,t_VR)

%Nicolai vertical tail weight estimate

y = 98.5*(((((W_TO*N)/(10^5))^0.87)*...
    ((S_V/100)^1.2)*((b_V/t_VR)^0.5))^0.458);

end
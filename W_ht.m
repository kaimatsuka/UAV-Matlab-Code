function y = W_ht(W_TO,N,l_T,S_H,b_H,t_HR)

%Nicolai horizontal tail weight estimate

y = 127*(((((W_TO*N)/(10^5))^0.87)*((l_T/10)^0.483)*...
    ((S_H/100)^1.2)*((b_H/t_HR)^0.5))^0.458);

end

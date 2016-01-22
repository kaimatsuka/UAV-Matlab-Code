function y = W_w(lam, S_W, W_TO, N, A, lam_q, mtr, V_e)

%Nicolai wing weight estimate

y = 96.948*(((((W_TO*N)/(10^5))^0.65)*((A/cos(lam_q))^0.57)*...
((S_W/100)^0.61)*(((1+lam)/(2*mtr))^0.36)*((1+(V_e/500))^0.5))^0.993);

end
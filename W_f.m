function y = W_f(W_TO,N,L,W,D,V_e)

%Nicolai fuselage weight estimate

y = 200*(((((W_TO*N)/(10^5))^0.286)*((L/10)^0.857)*...
    ((W+D)/10)*((V_e/100)^0.338))^1.1);

end
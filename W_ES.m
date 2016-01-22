function y = W_ES(W_FS,W_TRON)

%Nicolai electrical system weight estimate

y = 426*(((W_FS+W_TRON)/10^3)^0.51);

end
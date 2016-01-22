function y = C_Dpi(K,Q,C_f,S_wet,S_ref,C_Dmisc,C_DLP)

for i=1:length(K)
    
    y(i) = (K(i)*Q(i)*C_f(i)*S_wet(i))/S_ref;
    
end

y = sum(y);

y = y+C_Dmisc+C_DLP;

end
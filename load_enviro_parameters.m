%Environmental Parameters 

gamma = 1.4; %specific heat ratio
rho_low  = 0.002309; %low altitude density (slugs/ft^3) 1000 ft
rho_high = 0.0018985; %low altitude density (slugs/ft^3) 7500 ft
rho_avg  = (rho_low+rho_high)/2; %average density (slugs/ft^3)
R        = 1716.59; %gas constant (ft^2/s^2*R)
a_s      = 1107.815608; %speed of sound (ft/s) 50F
temp     = 50; %temperature (fahrenheit)
mu       = 3.82*10^-7; %dynamic viscosity of air 

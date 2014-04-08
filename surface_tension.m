% This will calculate the surface tension of water at a given T. T in K,
% sigma will be in N/m

function sigma_out = surface_tension(T)

%sigma_out = 0.07275*(1-0.002*(T-291));

Tc  =  647.096;
tau =  1-T/Tc;
B   =  0.2358;
b   = -0.625;
mu  =  1.256;

sigma_out = B*tau.^mu.*(1+b.*tau);
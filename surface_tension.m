% This will calculate the surface tension of water at a given T. T in K,
% sigma will be in N/m

function sigma = surface_tension(T)

    Tc  =  647.096;
    tau =  1-T/Tc;
    B   =  0.2358;
    b   = -0.625;
    mu  =  1.256;

    sigma = B*tau.^mu.*(1+b.*tau);
    
    return
end

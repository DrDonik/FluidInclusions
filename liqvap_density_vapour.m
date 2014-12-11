function rho = liqvap_density_vapour(T)

%Calculate the density on the liquid-vapor equilibrium curve (saturated vapor)
%Refer to equation (2.7) of Wagner and Pruss, J. Phys. Chem. Ref. Data, 31, p 387 (2002)

rhoc = 322; %Density of water in MKS units at the critical point
Tc = 647.096; %Temperature in K at the critical point of water

%Coefficients in the equation for the pressure on the liquid-vapor
%equilibrium curve
c1 = -2.03150240;
c2 = -2.68302940;
c3 = -5.38626492;
c4 = -17.2991605;
c5 = -44.7586581;
c6 = -63.9201063;

x = 1 - T/Tc;

rho =( rhoc * exp(c1 * x.^(2/6) + c2 * x.^(4/6) + c3 * x.^(8/6) + c4 * x.^(18/6) + c5 * x.^(37/6) + c6 * x.^(71/6)))/1000;

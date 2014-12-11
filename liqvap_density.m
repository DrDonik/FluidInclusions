function rho = liqvap_density(T)

% T in K, rho will be in g/cm³. Multiply by 1000 if you want kg/m³

%Calculate the density on the liquid-vapor equilibrium curve
%Refer to equation (2.6) of Wagner and Pruss, J. Phys. Chem. Ref. Data, 31, p 387 (2002)

rhoc = 322; %Density of water in MKS units at the critical point
Tc = 647.096; %Temperature in K at the critical point of water

%Coefficients in the equation for the pressure on the liquid-vapor
%equilibrium curve
b1 = 1.99274064;
b2 = 1.09965342;
b3 = -0.510839303;
b4 = -1.75493479;
b5 = -45.5170352;
b6 = -6.74694450e5;

x = 1 - T/Tc;

rho = (rhoc * (1 + b1 * x.^(1/3) + b2 * x.^(2/3) + b3 * x.^(5/3) + b4 * x.^(16/3) + b5 * x.^(43/3) + b6 * x.^(110/3)))/1000;

%This sum is not valid close to 4°C. Test for the difference to the
%equilibrium:

%p_isochore = pwater(rho*1000,T);
%p_water = liqvap(T);

%p_diff = p_water - p_isochore

function [reftemp, alpha_V] = expansion_coeff(T)

%The formulas are for T in °C. I expect the argument T given in K
%go to °C
T = T - 273.15;

mineralNumber = get_fi_mineral();

if mineralNumber == 1;
    %disp('Calcite, reference T = 30°, Rao')
    
    reftemp = 30;
    alpha_c = 24.670e-6 + 1.742e-8*T - 5.141e-12*T.^2;

    alpha_a = -3.660e-6 - 7.112e-10*T - 3.339e-12*T.^2;

    alpha_V = 2*alpha_a + alpha_c;
end;

if mineralNumber == 2;
    %disp('Calcite, reference T = 20° Rao, Maximum Density at 5.15 °C')
    
    reftemp = 20;
    alpha_V = (17.229 + 0.016093*T - 1.1813e-5*T.^2)*1e-6;
end;

if mineralNumber == 3;
    %disp('Calcite, reference T = 20°, TPPM, Maximum Density at 4.65 °C')
    
    reftemp = 20;
    alpha_V = (10.996 + 0.030161*T - 2.3376e-5*T.^2)*1e-6;
end;

if mineralNumber == 4;
    %disp('Quartz, reference T = 20° TPPM, Maximum Density at 6.05 °C')
    
    reftemp = 20;
    alpha_V = (30.484 + 0.033128*T +3.9391e-5*T.^2)*1e-6;
end;

if mineralNumber == 5;
    %disp('Quartz, reference T = 20° HTMIAC, Maximum Density at 6.15 °C')
    
    reftemp = 20;
    alpha_V = (33.208 + 0.067459*T - 1.4725e-4*T.^2)*1e-6;
end;

if mineralNumber == 6;
    %disp('Gypsum, reference T = 19.4° Schofield, Maximum Density at 8.35 °C')
    
    reftemp = 19.4;
    alpha_V = (64.719 + 0.081573*T + 2.9088e-4*T.^2 + 1.2963e-6*T.^3)*1e-6;
end;

if mineralNumber == 7;
    %disp('No correction applied')
    
    reftemp = zeros(size(T));
    alpha_V = zeros(size(T));
end;

return;
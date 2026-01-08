% defines the Planck function
% param1: Ts is the surface temperature
% param2: r is the equivalent wavelength


function [B_Ts] = Function_B_Ts(Ts,r)
    C1=1.191*10^8;
    C2=1.439*10^4;
    a = r*r*r*r*r;
    b = exp(C2./(r.*Ts))-1;
    B_Ts = C1./(a.*b);
end
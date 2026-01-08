%% This function is used to inversely calculate Planck's formula for TES 
% param1: Ts is the surface temperature
% param2: r is the equivalent wavelength
% emiss: LSE, input the maximum emissivity in TES
% C1 and C2 are constants in Planck's formula

function [Ts] = Function_inverse_B_Ts_tes(B_Ts, r,emiss)
    C1=1.191*10^8;
    C2=1.439*10^4;
    Ts = (C2/r)./(log((C1.*emiss)./(B_Ts.*r.*r.*r.*r.*r)+1));
end
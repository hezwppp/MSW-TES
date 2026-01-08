
%% This function is used to calculate the Planck formula in reverse
% Ts: the surface temperature 
% r: the equivalent wavelength 
% e: LSE,
% C1, C2 are the constants in the Planck formula

function [Ts] = Function_inverse_B_Ts(B_Ts, r)
    C1=1.191*10^8;
    C2=1.439*10^4;
    Ts = (C2/r)./(log(C1./(B_Ts.*r.*r.*r.*r.*r)+1));
end
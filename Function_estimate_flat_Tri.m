
function [Tri] = Function_estimate_flat_Tri(Tgi, ewl, Ld, e_max)

% tb_Bi, the brightness temperature of Channel 1
% tb_Bj, the bright temperature of channel 2
% beta_BiBj, SW coefficient calculated using two channels
% ewl, effective wavelength
% Ld, atmospheric downflow radiation
% e_max, maximum emissivity
% author: HZW

    C1=1.191*10^8;
    C2=1.439*10^4;

    B_Tg = C1/(ewl*ewl*ewl*ewl*ewl*(exp(C2/(ewl*Tgi))-1));
    ri = B_Tg - (1 - e_max).*Ld;
    Tri = (C2/ewl)/(log((C1*e_max)/(ri*ewl*ewl*ewl*ewl*ewl)+1));
end


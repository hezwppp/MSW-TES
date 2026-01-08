
function [Tri] = Function_estimate_mountain_Tri(Tgi, ewl, Ld, SVF, e_max)
% Calculate the surface radiant temperature considering the T-A effect
% tb_Bi, the brightness temperature of Channel 1
% tb_Bj, the bright temperature of channel 2
% beta_BiBj, SW coefficient calculated using two channels
% ewl, effective wavelength
% Ld, atmospheric downflow radiation
% s, spherical albedo at the bottom of the atmosphere
% e_max, maximum emissivity
% author: HZW

    C1=1.191*10^8;
    C2=1.439*10^4;

    B_Tg = C1/(ewl*ewl*ewl*ewl*ewl*(exp(C2/(ewl*Tgi))-1));
    Rterrain = e_max.*B_Tg; % 假设邻近像元的离地辐射温度和目标像元的一致
    ri = B_Tg - (1-e_max).*Ld.*SVF - (1-e_max).*(1-SVF).*Rterrain;
    Tri = (C2/ewl)/(log((C1.*e_max)/(ri*ewl*ewl*ewl*ewl*ewl)+1));

end


%% Calculate the ground-level radiant temperature (Tgi)
% tb_Bi, the brightness temperature of Channel 1
% tb_Bj, the bright temperature of channel 2
% beta_BiBj, the SW coefficient calculated using two channels, has five coefficients
% author: HZW

function [Tgi] = Function_estimate_Tgi(tb_Bi, tb_Bj, beta_BiBj)
    co_Tg = beta_BiBj;

    Tgi = co_Tg(1) + co_Tg(2).*tb_Bi + co_Tg(3).*(tb_Bi - tb_Bj)...
            + co_Tg(4).*(tb_Bi - tb_Bj).*(tb_Bi - tb_Bj) + co_Tg(5).*(tb_Bi - tb_Bj)./(tb_Bi + tb_Bj);
end


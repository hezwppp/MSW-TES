
%% The MTES algorithm with considering the topographic and adjacent effects (T-A effect)

%%%%%%%%%%%%% input parameters %%%%%%%%%%%%%%%%%%%%
% SSP: Small-scale self-heat parameter
% SVF: Sky View factor
% tb_01~tb_05：The brightness temperature data of five TIR channels from ASTER sensor,respectively
% LST: LST retrieved by SW-TES algorithm without considering T-A effect
% LSE_01~LSE_05: LSE of five TIR bands retrieved by SW-TES algorithm without considering T-A effect
% ewl：The effective wavelength containing five TIR channels from ASTER sensor
% Ld_01~Ld_05: The atmospheric downward radiation containing five TIR channels from ASTER sensor
% LUT_3D：A SW coefficient lookup table constructed with considering T-A effect
% beta_ssp：The empirical relationship of εmin-MMD with considering topographic effect for different SSP
% e_max_list: The maximum LSE of each band retrieved by SW-TES algorithm without considering T-A effect
% band_number：Number of bands

%%%%%%%%%%%%% return parameters %%%%%%%%%%%%%%%%%%%%
% MLST: The MLST retrieved by the MSW-TES agorithm with considering T-A effect
% MLSE_01~MLSE_05: The MLSE of five TIR channels retrieved by the MSW-TES agorithm from ASTER sensor 
%                with considering T-A effect,respectively,

% date: 2025-07-17
% author: Zhi-Wei He

function [MLST,MLSE_01,MLSE_02,MLSE_03,MLSE_04,MLSE_05] = ...
    Function_TES_yes_TAeffect_satllite_data(SSP,SVF,tb_01,tb_02,tb_03,tb_04,tb_05,...
            LST,LSE_01,LSE_02,LSE_03,LSE_04,LSE_05,ewl,Ld_01,Ld_02,Ld_03,Ld_04,Ld_05,...
            LUT_3D,beta_ssp,e_max_list,bands_number)

% Obtain the rows and columns of the data
[Row,Col] = size(tb_01);

% Read in the effective wavelengths of five channels
ewl_01 = ewl(1);
ewl_02 = ewl(2);
ewl_03 = ewl(3);
ewl_04 = ewl(4);
ewl_05 = ewl(5);

% Create a zero matrix of LST and LSE to store the results
MLST = zeros(Row,Col);
MLSE_01 = zeros(Row,Col);
MLSE_02 = zeros(Row,Col);
MLSE_03 = zeros(Row,Col);
MLSE_04 = zeros(Row,Col);
MLSE_05 = zeros(Row,Col);

Tg1 = [];
Tr1 = [];
for row = 1:Row
    mrow = row
    for col = 1:Col
        col;
        if isnan(tb_01(row,col)) || LST(row,col)==0 % If the pixel is empty, exit this loop and proceed to the next one
            continue;
        end
        
        e_max = e_max_list(row,col);

        Ld_01_ =  Ld_01(row,col);
        Ld_02_ =  Ld_02(row,col);
        Ld_03_ =  Ld_03(row,col);
        Ld_04_ =  Ld_04(row,col);
        Ld_05_ =  Ld_05(row,col);

        % The sw fitting coefficients corresponding to ssp and svf are obtained through the interpolation function
        ssp = SSP(row,col);
        svf = SVF(row,col);
        beta_B1B2 = Function_interpolate_SW_coefficient(LUT_3D,ssp,svf, 1);
        beta_B2B1 = Function_interpolate_SW_coefficient(LUT_3D,ssp,svf, 2);
        beta_B3B1 = Function_interpolate_SW_coefficient(LUT_3D,ssp,svf, 3);
        beta_B4B5 = Function_interpolate_SW_coefficient(LUT_3D,ssp,svf, 4);
        beta_B5B4 = Function_interpolate_SW_coefficient(LUT_3D,ssp,svf, 5);

        %% execute NEM module 
        % Use the SW coefficient to calculate the above-ground radiation temperatures of the five channels
        Tg_01 = Function_estimate_Tgi(tb_01(row,col), tb_02(row,col), beta_B1B2);
        Tg_02 = Function_estimate_Tgi(tb_02(row,col), tb_01(row,col), beta_B2B1);
        Tg_03 = Function_estimate_Tgi(tb_03(row,col), tb_01(row,col), beta_B3B1);
        Tg_04 = Function_estimate_Tgi(tb_04(row,col), tb_05(row,col), beta_B4B5);
        Tg_05 = Function_estimate_Tgi(tb_05(row,col), tb_04(row,col), beta_B5B4);
        Tg1 =[Tg1;Tg_01,Tg_02,Tg_03,Tg_04,Tg_05];

        T_d = 1; % Assign an initial value of 1 K to the initial comparison temperature
        for num = 1:12
            Tr_01 = Function_estimate_mountain_Tri(Tg_01,ewl_01,Ld_01_,SVF(row,col),e_max);
            Tr_02 = Function_estimate_mountain_Tri(Tg_02,ewl_02,Ld_02_,SVF(row,col),e_max);
            Tr_03 = Function_estimate_mountain_Tri(Tg_03,ewl_03,Ld_03_,SVF(row,col),e_max);
            Tr_04 = Function_estimate_mountain_Tri(Tg_04,ewl_04,Ld_04_,SVF(row,col),e_max);
            Tr_05 = Function_estimate_mountain_Tri(Tg_05,ewl_05,Ld_05_,SVF(row,col),e_max);
            Tr1 =[Tr1;Tr_01,Tr_02,Tr_03,Tr_04,Tr_05];
            Tr_n1= [Tr_01;Tr_02;Tr_03;Tr_04;Tr_05];

            %  Calculate T_NEM
            T_NEM = max(Tr_n1);
            % Calculate the emissivity of five channels
            e_01_ = Function_B_Ts(Tr_01,ewl_01)/Function_B_Ts(T_NEM,ewl_01);
            e_02_ = Function_B_Ts(Tr_02,ewl_02)/Function_B_Ts(T_NEM,ewl_02);
            e_03_ = Function_B_Ts(Tr_03,ewl_03)/Function_B_Ts(T_NEM,ewl_03);
            e_04_ = Function_B_Ts(Tr_04,ewl_04)/Function_B_Ts(T_NEM,ewl_04);
            e_05_ = Function_B_Ts(Tr_05,ewl_05)/Function_B_Ts(T_NEM,ewl_05);
            e_n = [e_01_;e_02_;e_03_;e_04_;e_05_];
    
            %% execute RAT module 
            beta_01 = e_01_/(sum(e_n)/bands_number);
            beta_02 = e_02_/(sum(e_n)/bands_number);
            beta_03 = e_03_/(sum(e_n)/bands_number);
            beta_04 = e_04_/(sum(e_n)/bands_number);
            beta_05 = e_05_/(sum(e_n)/bands_number);
            beta_01_to_05 = [beta_01,beta_02,beta_03,beta_04,beta_05];
    
            %% execute MMD module   
            % The regression coefficients a,b, and c of the MMD module which needs to be re-determined based on the value of SSP

            if (ssp >= 0.0) &&  (ssp < 0.1)
                beta = beta_ssp(1,:);
            end
            if (ssp  >= 0.1) &&  (ssp  < 0.2)
                beta = beta_ssp(2,:);
            end
            if (ssp  >= 0.2) &&  (ssp  < 0.3)
                beta = beta_ssp(3,:);
            end
            if (ssp  >= 0.3) &&  (ssp  < 0.4)
                beta = beta_ssp(4,:);
            end
            if (ssp  >= 0.4) &&  (ssp  < 0.5)
                beta = beta_ssp(5,:);
            end
            if (ssp  >= 0.5) &&  (ssp  < 0.6)
                beta = beta_ssp(6,:);
            end
            if (ssp  >= 0.6) &&  (ssp  < 0.7)
                beta = beta_ssp(7,:);
            end
            if (ssp  >= 0.7) &&  (ssp  < 0.8)
                beta = beta_ssp(8,:);
            end
            if (ssp  >= 0.8) &&  (ssp  < 0.9)
                beta = beta_ssp(9,:);
            end
            if (ssp  >= 0.9) &&  (ssp  < 1.0)
                beta = beta_ssp(10,:);
            end
            if ssp  == 1.0
                beta = beta_ssp(11,:);
            end

            MMD_n = max(beta_01_to_05) - min(beta_01_to_05);
            e_min_n = Function_MMD(beta,MMD_n);           

            % Calculate the LSE of the five channels
            e_01 = (e_min_n/min(beta_01_to_05))*beta_01;
            e_02 = (e_min_n/min(beta_01_to_05))*beta_02;
            e_03 = (e_min_n/min(beta_01_to_05))*beta_03;
            e_04 = (e_min_n/min(beta_01_to_05))*beta_04;
            e_05 = (e_min_n/min(beta_01_to_05))*beta_05;

            % Obtain the maximum LSE among the five channels and update the maximum LSE in the iteration
            e_01_to_05 = [e_01,e_02,e_03,e_04,e_05];
            e_max_new = max(e_01_to_05);
            index = find(e_01_to_05 == e_max_new); % Find the channel corresponding to the maximum emissivity
            ch_name = num2str(index,'%02d');
    
            % obtain MLST
            eval(['Tg_b','=','Tg_',ch_name,';']);
            eval(['wl','=','ewl_',ch_name,';']);
            eval(['ld','=','Ld_',ch_name,'_',';']); 
            eval(['lse','=','LSE_',ch_name,'(row,col)',';']);

            % Calculate the mean LST and LSE of eight adjacent pixels to calculate the adjacent effect
            lst_temp = 0;
            lse_temp = 0;
            n = 0;
            
            % situation 1
            if (row-1 ~= 0) && (col-1 ~= 0)
                lst_temp = lst_temp + LST(row-1,col-1);
                eval(['lse_temp','=','lse_temp + ','LSE_',ch_name,'(row-1,col-1)',';']);
                n = n+1;
            end

            % situation 2
            if row-1 ~= 0
                lst_temp = lst_temp + LST(row-1,col);
                eval(['lse_temp','=','lse_temp + ','LSE_',ch_name,'(row-1,col)',';']);
                n = n+1;
            end

            % situation 3
            if (row-1 ~= 0) && (col+1 < Col)
                lst_temp = lst_temp + LST(row-1,col+1);
                eval(['lse_temp','=','lse_temp + ','LSE_',ch_name,'(row-1,col+1)',';']);
                n = n+1;
            end

            % situation 4
            if col-1 ~= 0
                lst_temp = lst_temp + LST(row,col-1);
                eval(['lse_temp','=','lse_temp + ','LSE_',ch_name,'(row,col-1)',';']);
                n = n+1;
            end

            % situation 5
            if (row ~= 0) && (col+1 < Col)
                lst_temp = lst_temp + LST(row,col+1);
                eval(['lse_temp','=','lse_temp + ','LSE_',ch_name,'(row,col+1)',';']);
                n = n+1;
            end

            % situation 6
            if (row+1 < Row) && (col-1 ~= 0)
                lst_temp = lst_temp + LST(row+1,col-1);
                eval(['lse_temp','=','lse_temp + ','LSE_',ch_name,'(row+1,col-1)',';']);
                n = n+1;
            end

            % situation 7
            if (row+1 < Row)
                lst_temp = lst_temp + LST(row+1,col);
                eval(['lse_temp','=','lse_temp + ','LSE_',ch_name,'(row+1,col)',';']);
                n = n+1;
            end

            % situation 8
            if (row+1 < Row) && (col+1 < Col)
                lst_temp = lst_temp + LST(row+1,col+1);
                eval(['lse_temp','=','lse_temp + ','LSE_',ch_name,'(row+1,col+1)',';']);
                n = n+1;
            end

            if n > 0
                ave_LST = lst_temp / n;
                ave_LSE = lse_temp / n;
            else
                ave_LST = LST(row, col);
                eval(['ave_LSE = LSE_', ch_name, '(row, col);']);
            end

            % leaving ground radiation
            B_Tg_b = Function_B_Ts(Tg_b,wl);

            % topographic effect
            B_T_terrain = Function_B_Ts(ave_LST,wl);
            Rterrain = ave_LSE*B_T_terrain; % The LST is used as the initial value for calculating the adjacent radiation

            % calculate Tr_
            Tr_ = (B_Tg_b...
                   -(1-e_max_new)*ld*SVF(row,col)...
                   -(1-e_max_new)*Rterrain*(1-SVF(row,col)))...
                   /e_max_new;

            LST_ = Function_inverse_B_Ts_tes(Tr_, wl, e_max_new);
            e_max_list(row,col) = e_max_new;

            % Separate out LST and LSE
            MLSE_01(row,col) = e_01;
            MLSE_02(row,col) = e_02;
            MLSE_03(row,col) = e_03;
            MLSE_04(row,col) = e_04;
            MLSE_05(row,col) = e_05;

            MLST(row,col) = LST_;
            delta_T = abs(LST_ - T_d);
            if delta_T <= 0.1
                delta_T;
                num;
                break;
            end
            e_max = e_max_new;
            T_d = LST_;
        end
    end
%     if row==100
%         break;
%     end
end
end
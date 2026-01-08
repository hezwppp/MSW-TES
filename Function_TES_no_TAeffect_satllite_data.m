
%% The TES algorithm without considering the topographic and adjacent effects (T-A effect)

%%%%%%%%%%%%% input parameters %%%%%%%%%%%%%%%%%%%%
% tb_01~tb_05：The brightness temperature data of five TIR channels from ASTER sensor,respectively
% ewl：The effective wavelength containing five TIR channels from ASTER sensor
% Ld_01~Ld_05: The atmospheric downward radiation containing five TIR channels from ASTER sensor
% LUT_2D：A SW coefficient lookup table constructed without considering T-A effect
% beta：The empirical relationship of εmin-MMD without considering topographic effect
% band_number：Number of bands

%%%%%%%%%%%%% return parameters %%%%%%%%%%%%%%%%%%%%
% LST0: The LST retrieved by the SW-TES agorithm without considering T-A effect,respectively
% LSE_01~LSE_05: The LSE of five TIR channels retrieved by the SW-TES agorithm from ASTER sensor 
%                without considering T-A effect,respectively,
% e_max_list: 

% date: 2025-07-17
% author: Zhi-Wei He


function [LST0,LSE_01,LSE_02,LSE_03,LSE_04,LSE_05,e_max_list] = ...
   Function_TES_no_TAeffect_satllite_data(tb_01,tb_02,tb_03,tb_04,tb_05,ewl,...
            Ld_01,Ld_02,Ld_03,Ld_04,Ld_05,LUT_2D,beta,band_number)


% Obtain the rows and columns of the data
[Row,Col] = size(tb_01);

% Read in the effective wavelengths of five channels
ewl_01 = ewl(1);
ewl_02 = ewl(2);
ewl_03 = ewl(3);
ewl_04 = ewl(4);
ewl_05 = ewl(5);

% Read in the SW coefficients of the five channels
beta_B1B3 = LUT_2D(1,2:6);
beta_B2B1 = LUT_2D(2,2:6);
beta_B3B1 = LUT_2D(3,2:6);
beta_B4B5 = LUT_2D(4,2:6);
beta_B5B4 = LUT_2D(5,2:6);

% Create a zero matrix of LST and LSE to store the results
LST0 = zeros(Row,Col);
LSE_01 = zeros(Row,Col);
LSE_02 = zeros(Row,Col);
LSE_03 = zeros(Row,Col);
LSE_04 = zeros(Row,Col);
LSE_05 = zeros(Row,Col);

e_max_list = zeros(Row,Col);
Tg0 = [];
Tr0 = [];

%% Start to execute the SW-TES algorithm
for row = 1:Row
    row    
    for col = 1:Col
        e_max = 0.99;
        if isnan(tb_01(row,col)) % If the pixel is empty, exit this loop and proceed to the next one
            continue;
        end
        Ld_01_ =  Ld_01(row,col);
        Ld_02_ =  Ld_02(row,col);
        Ld_03_ =  Ld_03(row,col);
        Ld_04_ =  Ld_04(row,col);
        Ld_05_ =  Ld_05(row,col);

        %% execute NEM module 
        % Use the SW coefficient to calculate the above-ground radiation temperatures of the five channels
        Tg_01 = Function_estimate_Tgi(tb_01(row,col), tb_03(row,col), beta_B1B3);
        Tg_02 = Function_estimate_Tgi(tb_02(row,col), tb_01(row,col), beta_B2B1);
        Tg_03 = Function_estimate_Tgi(tb_03(row,col), tb_01(row,col), beta_B3B1);
        Tg_04 = Function_estimate_Tgi(tb_04(row,col), tb_05(row,col), beta_B4B5);
        Tg_05 = Function_estimate_Tgi(tb_05(row,col), tb_04(row,col), beta_B5B4);
        Tg0 =[Tg0;Tg_01,Tg_02,Tg_03,Tg_04,Tg_05];

        T_d = 1; % Assign an initial value of 1 K to the initial comparison temperature
        for num = 1:12
            Tr_01 = Function_estimate_flat_Tri(Tg_01, ewl_01, Ld_01_, e_max);
            Tr_02 = Function_estimate_flat_Tri(Tg_02, ewl_02, Ld_02_, e_max);
            Tr_03 = Function_estimate_flat_Tri(Tg_03, ewl_03, Ld_03_, e_max);
            Tr_04 = Function_estimate_flat_Tri(Tg_04, ewl_04, Ld_04_, e_max);
            Tr_05 = Function_estimate_flat_Tri(Tg_05, ewl_05, Ld_05_, e_max);
            Tr0 =[Tr0;Tr_01,Tr_02,Tr_03,Tr_04,Tr_05];
            Tr_n0= [Tr_01;Tr_02;Tr_03;Tr_04;Tr_05];

            % Calculate T_NEM
            T_NEM = max(Tr_n0);
            % Calculate the emissivity of five channels
            e_01_ = Function_B_Ts(Tr_01,ewl_01)/Function_B_Ts(T_NEM,ewl_01);
            e_02_ = Function_B_Ts(Tr_02,ewl_02)/Function_B_Ts(T_NEM,ewl_02);
            e_03_ = Function_B_Ts(Tr_03,ewl_03)/Function_B_Ts(T_NEM,ewl_03);
            e_04_ = Function_B_Ts(Tr_04,ewl_04)/Function_B_Ts(T_NEM,ewl_04);
            e_05_ = Function_B_Ts(Tr_05,ewl_05)/Function_B_Ts(T_NEM,ewl_05);
            e_n0 = [e_01_;e_02_;e_03_;e_04_;e_05_];
    
            %% execute RAT module 
            beta_01 = e_01_/(sum(e_n0)/band_number);
            beta_02 = e_02_/(sum(e_n0)/band_number);
            beta_03 = e_03_/(sum(e_n0)/band_number);
            beta_04 = e_04_/(sum(e_n0)/band_number);
            beta_05 = e_05_/(sum(e_n0)/band_number);
            beta_01_to_05 = [beta_01,beta_02,beta_03,beta_04,beta_05];
    
            %% execute MMD module  
            % The regression coefficients a,b, and c of the MMD module can be obtained by fitting the pop library data
            MMD_n = max(beta_01_to_05) - min(beta_01_to_05);
            e_min_n0 = Function_MMD(beta,MMD_n);         

            % Calculate the LSE of the five channels
            e_01 = (e_min_n0/min(beta_01_to_05))*beta_01;
            e_02 = (e_min_n0/min(beta_01_to_05))*beta_02;
            e_03 = (e_min_n0/min(beta_01_to_05))*beta_03;
            e_04 = (e_min_n0/min(beta_01_to_05))*beta_04;
            e_05 = (e_min_n0/min(beta_01_to_05))*beta_05;

            % Obtain the maximum emissivity among the five channels and update the maximum emissivity in the iteration
            e_01_to_05 = [e_01,e_02,e_03,e_04,e_05];
            e_max_new = max(e_01_to_05);
            index = find(e_01_to_05 == e_max_new); % Find the channel corresponding to the maximum emissivity
    
            % obtain LST
            eval(['Tg_b','=','Tg_',num2str(index,'%02d'),';']);
            eval(['wl','=','ewl_',num2str(index,'%02d'),';']);
            eval(['ld','=','Ld_',num2str(index,'%02d'),'_',';']);
            B_Tg_b = Function_B_Ts(Tg_b,wl);
            Tr_ = (B_Tg_b - (1-e_max_new)*ld)/e_max_new;
            LST_ = Function_inverse_B_Ts_tes(Tr_,wl, e_max_new);
            e_max_list(row,col) = e_max_new;

            if ~isreal(e_01)
                row
                col
            end

            LSE_01(row,col) = e_01;
            LSE_02(row,col) = e_02;
            LSE_03(row,col) = e_03;
            LSE_04(row,col) = e_04;
            LSE_05(row,col) = e_05;

            LST0(row,col) = LST_;
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
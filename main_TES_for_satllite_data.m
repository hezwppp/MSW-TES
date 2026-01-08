% The implementation of the MSW-TES algorithm on satellite data
% author: He Z.-W

% clc;
% clear;

bands_number = 5; % Five thermal infrared bands
path_mat = '...\data\';
load([path_mat,'LUT_3D.mat']);
load([path_mat,'LUT_2D.mat']);
load([path_mat,'beta_flat.mat'])
load([path_mat,'beta_ssp.mat']); 

% Calculate the Tgi of each band using the SW coefficient obtained from the band combination
% The effective wavelength of 5 channels
ewl_01 = 8.2815;  
ewl_02 = 8.6330;  
ewl_03 = 9.0792; 
ewl_04 = 10.6621; 
ewl_05 = 11.2929;
ewl = [ewl_01,ewl_02,ewl_03,ewl_04,ewl_05];

%% Read in the brightness temperature data of the five channels of ASTER

% Radiation calibration coefficient of five ASTER TIR bands
co_b10 = 6.822*10^(-3); co_b11 = 6.780*10^(-3); co_b12 = 6.590*10^(-3);
co_b13 = 5.693*10^(-3); co_b14 = 5.225*10^(-3);


path_rad = 'tif\';
rad01_ = double(readgeoraster([path_rad,'B10.tif']));
rad02_ = double(readgeoraster([path_rad,'B11.tif']));
rad03_ = double(readgeoraster([path_rad,'B12.tif']));
rad04_ = double(readgeoraster([path_rad,'B13.tif']));
rad05_ = double(readgeoraster([path_rad,'B14.tif']));

rad01_(rad01_<=0)=nan;rad02_(rad02_<=0)=nan;rad03_(rad03_<=0)=nan;
rad04_(rad04_<=0)=nan;rad05_(rad05_<=0)=nan;
index = find(isnan(rad01_) | isnan(rad02_) | isnan(rad03_) | isnan(rad04_) | isnan(rad05_));
rad01_(index)=nan;rad02_(index)=nan;rad03_(index)=nan;rad04_(index)=nan;rad05_(index)=nan;
rad01 = (rad01_ - 1)*co_b10;
rad02 = (rad02_ - 1)*co_b11;
rad03 = (rad03_ - 1)*co_b12;
rad04 = (rad04_ - 1)*co_b13;
rad05 = (rad05_ - 1)*co_b14;

% The radiance is converted to brightness temperature
tb_01 = Function_inverse_B_Ts(rad01, ewl_01);
tb_02 = Function_inverse_B_Ts(rad02, ewl_02);
tb_03 = Function_inverse_B_Ts(rad03, ewl_03);
tb_04 = Function_inverse_B_Ts(rad04, ewl_04);
tb_05 = Function_inverse_B_Ts(rad05, ewl_05);

% Read in SSP and SVF
SSP = double(readgeoraster('SSP.tif'));
SVF = double(readgeoraster('SVF.tif'));

% Read the downwelling radiation
ld_path = '...\';
Ld_01 = double(readgeoraster([ld_path,'Ld_01.tif']));
Ld_02 = double(readgeoraster([ld_path,'Ld_02.tif']));
Ld_03 = double(readgeoraster([ld_path,'Ld_03.tif']));
Ld_04 = double(readgeoraster([ld_path,'Ld_04.tif']));
Ld_05 = double(readgeoraster([ld_path,'Ld_05.tif']));
Ld_01(index)=nan;Ld_02(index)=nan;Ld_03(index)=nan;Ld_04(index)=nan;Ld_05(index)=nan;
Ld = [Ld_01,Ld_02,Ld_03,Ld_04,Ld_05];

%% Step1: Implementation of the TES algorithm without considering the T-A effect
[LST,LSE_01,LSE_02,LSE_03,LSE_04,LSE_05,e_max_list] = ...
    Function_TES_no_TAeffect_satllite_data(tb_01,tb_02,tb_03,...
            tb_04,tb_05,ewl,Ld_01,Ld_02,Ld_03,Ld_04,Ld_05,...
            LUT_2D,beta_flat,bands_number);
AA = max(max(LST));
LST(LST==0)=nan;LSE_01(LSE_01==0)=nan;LSE_02(LSE_02==0)=nan;LSE_03(LSE_03==0)=nan;
LSE_04(LSE_04==0)=nan;LSE_05(LSE_05==0)=nan;e_max_list(e_max_list==0)=nan;

%% Step1: Implementation of the TES algorithm with considering the T-A effect
[MLST,MLSE_01,MLSE_02,MLSE_03,MLSE_04,MLSE_05] = ...
    Function_TES_yes_TAeffect_satllite_data(SSP,SVF,tb_01,tb_02,tb_03,tb_04,tb_05,...
            LST,LSE_01,LSE_02,LSE_03,LSE_04,LSE_05,ewl,Ld_01,Ld_02,Ld_03,Ld_04,Ld_05,...
            LUT_3D,beta_ssp,e_max_list,bands_number);







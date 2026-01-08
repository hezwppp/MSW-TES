% TES 算法在卫星数据上的实现
% author: HZW
% date: 2024.0306

% clc;
% clear;
% [Ts] = Function_inverse_B_Ts(4.5, 8.2891, 0.993)

bands_number = 5; % 5个热红外波段
path_mat = 'C:\Users\HZW\Desktop\Phd_hzw\Code\matlab\paper5e_SW-TES_ASTER_major_revision_1st_ISPRS\data\';
load([path_mat,'LUT_3D.mat']);
load([path_mat,'LUT_2D.mat']);
load([path_mat,'beta_flat.mat']) % 最小发射率和MMD的拟合系数
load([path_mat,'beta_ssp.mat']); % 分组最小发射率和MMD的拟合系数

% 利用波段组合得到的SW系数计算各波段的Tgi
% 5个通道的有效波长
ewl_01 = 8.2815;  
ewl_02 = 8.6330;  
ewl_03 = 9.0792; 
ewl_04 = 10.6621; 
ewl_05 = 11.2929;
ewl = [ewl_01,ewl_02,ewl_03,ewl_04,ewl_05];

%% 读入 ASTER 5个通道的亮温数据
date_num = 6; % 1-dali_0405_2008; 2-dali_1121_2016; 3-lasa_new_1125_2015; 4-lasa_new_0404_2017
             % 5-kailash_0314_2014；6-kailash_0320_2022
             % 最后的论文只留下了1256

% 辐射定标系数
co_b10 = 6.822*10^(-3); co_b11 = 6.780*10^(-3); co_b12 = 6.590*10^(-3);
co_b13 = 5.693*10^(-3); co_b14 = 5.225*10^(-3);

if date_num == 1 
    path_rad = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\ASTER_Image\tif_sa\Dali\';
    [~,geoinfo] = readgeoraster([path_rad,  'AST_L1B_0405_2008_035808_B10.tif']);
    info = geotiffinfo([path_rad,  'AST_L1B_0405_2008_035808_B10.tif']);
    rad01_ = double(readgeoraster([path_rad,'AST_L1B_0405_2008_035808_B10.tif']));
    rad02_ = double(readgeoraster([path_rad,'AST_L1B_0405_2008_035808_B11.tif']));
    rad03_ = double(readgeoraster([path_rad,'AST_L1B_0405_2008_035808_B12.tif']));
    rad04_ = double(readgeoraster([path_rad,'AST_L1B_0405_2008_035808_B13.tif']));
    rad05_ = double(readgeoraster([path_rad,'AST_L1B_0405_2008_035808_B14.tif']));
    rad01_(rad01_<=0)=nan;rad02_(rad02_<=0)=nan;rad03_(rad03_<=0)=nan;
    rad04_(rad04_<=0)=nan;rad05_(rad05_<=0)=nan;
    index = find(isnan(rad01_) | isnan(rad02_) | isnan(rad03_) | isnan(rad04_) | isnan(rad05_));
    rad01_(index)=nan;rad02_(index)=nan;rad03_(index)=nan;rad04_(index)=nan;rad05_(index)=nan;
    rad01 = (rad01_ - 1)*co_b10;
    rad02 = (rad02_ - 1)*co_b11;
    rad03 = (rad03_ - 1)*co_b12;
    rad04 = (rad04_ - 1)*co_b13;
    rad05 = (rad05_ - 1)*co_b14;
end

if date_num ==2
    path_rad = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\ASTER_Image\tif_sa\Dali\';
    [~,geoinfo] = readgeoraster([path_rad,  'AST_L1B_1121_2016_035831_B10.tif']);
    rad01_ = double(readgeoraster([path_rad,'AST_L1B_1121_2016_035831_B10.tif']));
    rad02_ = double(readgeoraster([path_rad,'AST_L1B_1121_2016_035831_B11.tif']));
    rad03_ = double(readgeoraster([path_rad,'AST_L1B_1121_2016_035831_B12.tif']));
    rad04_ = double(readgeoraster([path_rad,'AST_L1B_1121_2016_035831_B13.tif']));
    rad05_ = double(readgeoraster([path_rad,'AST_L1B_1121_2016_035831_B14.tif']));
    rad01_(rad01_<=0)=nan;rad02_(rad02_<=0)=nan;rad03_(rad03_<=0)=nan;
    rad04_(rad04_<=0)=nan;rad05_(rad05_<=0)=nan;
    index = find(isnan(rad01_) | isnan(rad02_) | isnan(rad03_) | isnan(rad04_) | isnan(rad05_));
    rad01_(index)=nan;rad02_(index)=nan;rad03_(index)=nan;rad04_(index)=nan;rad05_(index)=nan;
    rad01 = (rad01_ - 1)*co_b10;
    rad02 = (rad02_ - 1)*co_b11;
    rad03 = (rad03_ - 1)*co_b12;
    rad04 = (rad04_ - 1)*co_b13;
    rad05 = (rad05_ - 1)*co_b14;
end

if date_num == 3
    path_rad = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\ASTER_Image\tif_sa\Lasa_new\';
    [~,geoinfo] = readgeoraster([path_rad,  'AST_L1B_1125_2015_155900_B10.tif']);
    rad01_ = double(readgeoraster([path_rad,'AST_L1B_1125_2015_155900_B10.tif']));
    rad02_ = double(readgeoraster([path_rad,'AST_L1B_1125_2015_155900_B11.tif']));
    rad03_ = double(readgeoraster([path_rad,'AST_L1B_1125_2015_155900_B12.tif']));
    rad04_ = double(readgeoraster([path_rad,'AST_L1B_1125_2015_155900_B13.tif']));
    rad05_ = double(readgeoraster([path_rad,'AST_L1B_1125_2015_155900_B14.tif']));
    rad01_(rad01_<=0)=nan; rad02_(rad02_<=0)=nan; rad03_(rad03_<=0)=nan;
    rad04_(rad04_<=0)=nan; rad05_(rad05_<=0)=nan;
    index = find(isnan(rad01_) | isnan(rad02_) | isnan(rad03_) | isnan(rad04_) | isnan(rad05_));
    rad01_(index)=nan;rad02_(index)=nan;rad03_(index)=nan;rad04_(index)=nan;rad05_(index)=nan;
    rad01 = (rad01_ - 1)*co_b10;
    rad02 = (rad02_ - 1)*co_b11;
    rad03 = (rad03_ - 1)*co_b12;
    rad04 = (rad04_ - 1)*co_b13;
    rad05 = (rad05_ - 1)*co_b14;
end

if date_num == 4
    path_rad = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\ASTER_Image\tif_sa\Lasa_new\';
    [~,geoinfo] = readgeoraster([path_rad,  'AST_L1B_0404_2017_155827_B10.tif']);
    rad01_ = double(readgeoraster([path_rad,'AST_L1B_0404_2017_155827_B10.tif']));
    rad02_ = double(readgeoraster([path_rad,'AST_L1B_0404_2017_155827_B11.tif']));
    rad03_ = double(readgeoraster([path_rad,'AST_L1B_0404_2017_155827_B12.tif']));
    rad04_ = double(readgeoraster([path_rad,'AST_L1B_0404_2017_155827_B13.tif']));
    rad05_ = double(readgeoraster([path_rad,'AST_L1B_0404_2017_155827_B14.tif']));
    rad01_(rad01_<=0)=nan;rad02_(rad02_<=0)=nan;rad03_(rad03_<=0)=nan;
    rad04_(rad04_<=0)=nan;rad05_(rad05_<=0)=nan;
    index = find(isnan(rad01_) | isnan(rad02_) | isnan(rad03_) | isnan(rad04_) | isnan(rad05_));
    rad01_(index)=nan;rad02_(index)=nan;rad03_(index)=nan;rad04_(index)=nan;rad05_(index)=nan;
    rad01 = (rad01_ - 1)*co_b10;
    rad02 = (rad02_ - 1)*co_b11;
    rad03 = (rad03_ - 1)*co_b12;
    rad04 = (rad04_ - 1)*co_b13;
    rad05 = (rad05_ - 1)*co_b14;
end

if date_num == 5 % kailash_0314_2014
    path_rad = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\ASTER_Image\tif_sa\Kailash\';
    [~,geoinfo] = readgeoraster([path_rad,  'AST_L1B_0314_2014_052907_B10.tif']);
    rad01_ = double(readgeoraster([path_rad,'AST_L1B_0314_2014_052907_B10.tif']));
    rad02_ = double(readgeoraster([path_rad,'AST_L1B_0314_2014_052907_B11.tif']));
    rad03_ = double(readgeoraster([path_rad,'AST_L1B_0314_2014_052907_B12.tif']));
    rad04_ = double(readgeoraster([path_rad,'AST_L1B_0314_2014_052907_B13.tif']));
    rad05_ = double(readgeoraster([path_rad,'AST_L1B_0314_2014_052907_B14.tif']));
    rad01_(rad01_<=0)=nan; rad02_(rad02_<=0)=nan; rad03_(rad03_<=0)=nan;
    rad04_(rad04_<=0)=nan; rad05_(rad05_<=0)=nan;
    index = find(isnan(rad01_) | isnan(rad02_) | isnan(rad03_) | isnan(rad04_) | isnan(rad05_));
    rad01_(index)=nan;rad02_(index)=nan;rad03_(index)=nan;rad04_(index)=nan;rad05_(index)=nan;
    rad01 = (rad01_ - 1)*co_b10;
    rad02 = (rad02_ - 1)*co_b11;
    rad03 = (rad03_ - 1)*co_b12;
    rad04 = (rad04_ - 1)*co_b13;
    rad05 = (rad05_ - 1)*co_b14;
end

if date_num == 6 % kailash_0320_2022
    path_rad = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\ASTER_Image\tif_sa\Kailash\';
    [~,geoinfo] = readgeoraster([path_rad,  'AST_L1B_0320_2022_052042_B10.tif']);
    rad01_ = double(readgeoraster([path_rad,'AST_L1B_0320_2022_052042_B10.tif']));
    rad02_ = double(readgeoraster([path_rad,'AST_L1B_0320_2022_052042_B11.tif']));
    rad03_ = double(readgeoraster([path_rad,'AST_L1B_0320_2022_052042_B12.tif']));
    rad04_ = double(readgeoraster([path_rad,'AST_L1B_0320_2022_052042_B13.tif']));
    rad05_ = double(readgeoraster([path_rad,'AST_L1B_0320_2022_052042_B14.tif']));
    rad01_(rad01_<=0)=nan;rad02_(rad02_<=0)=nan;rad03_(rad03_<=0)=nan;
    rad04_(rad04_<=0)=nan;rad05_(rad05_<=0)=nan;
    index = find(isnan(rad01_) | isnan(rad02_) | isnan(rad03_) | isnan(rad04_) | isnan(rad05_));
    rad01_(index)=nan;rad02_(index)=nan;rad03_(index)=nan;rad04_(index)=nan;rad05_(index)=nan;
    rad01 = (rad01_ - 1)*co_b10;
    rad02 = (rad02_ - 1)*co_b11;
    rad03 = (rad03_ - 1)*co_b12;
    rad04 = (rad04_ - 1)*co_b13;
    rad05 = (rad05_ - 1)*co_b14;
end

% 再将读入的辐亮度转为亮温
tb_01 = Function_inverse_B_Ts(rad01, ewl_01);% 在该函数中需要输入地表发射率，但是此处为1即可
tb_02 = Function_inverse_B_Ts(rad02, ewl_02);
tb_03 = Function_inverse_B_Ts(rad03, ewl_03);
tb_04 = Function_inverse_B_Ts(rad04, ewl_04);
tb_05 = Function_inverse_B_Ts(rad05, ewl_05);

% 读入 SSP 和 SVF
if date_num == 1 || date_num == 2
    SSP = double(readgeoraster('D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\DEM\dem_sa\dali_SSP_sa.tif'));
    SSP(SSP<-300) = nan; SSP(index) = nan;
    SVF = double(readgeoraster('D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\DEM\dem_sa\Dali_SVF_sa.tif'));
    SVF(SVF<-300) = nan; SVF(index) = nan;
end

if date_num == 3 || date_num == 4
    SSP = double(readgeoraster('D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\DEM\dem_sa\lasa_new_SSP_sa.tif'));
    SSP(SSP<-300) = nan; SSP(index) = nan;
    SVF = double(readgeoraster('D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\DEM\dem_sa\Lasa_new_SVF_sa.tif'));
    SVF(SVF<-300) = nan; SVF(index) = nan;
end

if date_num == 5 || date_num == 6
    SSP = double(readgeoraster('D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\DEM\dem_sa\Kailash_SSP_sa.tif'));
    SSP(SSP<-300) = nan; SSP(index) = nan;
    SVF = double(readgeoraster('D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\DEM\dem_sa\Kailash_SVF_sa.tif'));
    SVF(SVF<-300) = nan; SVF(index) = nan;
end

% issp = find(~isnan(SSP));
% isvf = find(~isnan(SVF));
% 读入大气下行辐射 % 1-dali_0405_2008; 2-dali_1121_2016; 3-lasa_0610_2013; 4-lasa_1011_2017
if date_num == 1 
    ld_path = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\EAR5\pixel_by_pixel\dali_0405_2008\';
    Ld_01 = double(readgeoraster([ld_path,'Ld_dali_0405_2008_sa_01.tif']));
    Ld_02 = double(readgeoraster([ld_path,'Ld_dali_0405_2008_sa_02.tif']));
    Ld_03 = double(readgeoraster([ld_path,'Ld_dali_0405_2008_sa_03.tif']));
    Ld_04 = double(readgeoraster([ld_path,'Ld_dali_0405_2008_sa_04.tif']));
    Ld_05 = double(readgeoraster([ld_path,'Ld_dali_0405_2008_sa_05.tif']));
    Ld_01(index)=nan;Ld_02(index)=nan;Ld_03(index)=nan;Ld_04(index)=nan;Ld_05(index)=nan;
    Ld = [Ld_01,Ld_02,Ld_03,Ld_04,Ld_05];
end
if date_num == 2 
    ld_path = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\EAR5\pixel_by_pixel\dali_1121_2016\';
    Ld_01 = double(readgeoraster([ld_path,'Ld_dali_1121_2016_sa_01.tif']));
    Ld_02 = double(readgeoraster([ld_path,'Ld_dali_1121_2016_sa_02.tif']));
    Ld_03 = double(readgeoraster([ld_path,'Ld_dali_1121_2016_sa_03.tif']));
    Ld_04 = double(readgeoraster([ld_path,'Ld_dali_1121_2016_sa_04.tif']));
    Ld_05 = double(readgeoraster([ld_path,'Ld_dali_1121_2016_sa_05.tif']));
    Ld_01(index)=nan;Ld_02(index)=nan;Ld_03(index)=nan;Ld_04(index)=nan;Ld_05(index)=nan;
    Ld = [Ld_01,Ld_02,Ld_03,Ld_04,Ld_05];
end
if date_num == 3 
    ld_path = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\EAR5\pixel_by_pixel\lasa_new_1125_2015\';
    Ld_01 = double(readgeoraster([ld_path,'Ld_lasa_new_1125_2015_sa_01.tif']));
    Ld_02 = double(readgeoraster([ld_path,'Ld_lasa_new_1125_2015_sa_02.tif']));
    Ld_03 = double(readgeoraster([ld_path,'Ld_lasa_new_1125_2015_sa_03.tif']));
    Ld_04 = double(readgeoraster([ld_path,'Ld_lasa_new_1125_2015_sa_04.tif']));
    Ld_05 = double(readgeoraster([ld_path,'Ld_lasa_new_1125_2015_sa_05.tif']));
    Ld_01(index)=nan;Ld_02(index)=nan;Ld_03(index)=nan;Ld_04(index)=nan;Ld_05(index)=nan;
    Ld = [Ld_01,Ld_02,Ld_03,Ld_04,Ld_05];
end
if date_num == 4 
    ld_path = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\EAR5\pixel_by_pixel\lasa_new_0404_2017\';
    Ld_01 = double(readgeoraster([ld_path,'Ld_lasa_new_0404_2017_sa_01.tif']));
    Ld_02 = double(readgeoraster([ld_path,'Ld_lasa_new_0404_2017_sa_02.tif']));
    Ld_03 = double(readgeoraster([ld_path,'Ld_lasa_new_0404_2017_sa_03.tif']));
    Ld_04 = double(readgeoraster([ld_path,'Ld_lasa_new_0404_2017_sa_04.tif']));
    Ld_05 = double(readgeoraster([ld_path,'Ld_lasa_new_0404_2017_sa_05.tif']));
    Ld_01(index)=nan;Ld_02(index)=nan;Ld_03(index)=nan;Ld_04(index)=nan;Ld_05(index)=nan;
    Ld = [Ld_01,Ld_02,Ld_03,Ld_04,Ld_05];
end

if date_num == 5 
    ld_path = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\EAR5\pixel_by_pixel\kailash_0314_2014\';
    Ld_01 = double(readgeoraster([ld_path,'Ld_kailash_0314_2014_sa_01.tif']));
    Ld_02 = double(readgeoraster([ld_path,'Ld_kailash_0314_2014_sa_02.tif']));
    Ld_03 = double(readgeoraster([ld_path,'Ld_kailash_0314_2014_sa_03.tif']));
    Ld_04 = double(readgeoraster([ld_path,'Ld_kailash_0314_2014_sa_04.tif']));
    Ld_05 = double(readgeoraster([ld_path,'Ld_kailash_0314_2014_sa_05.tif']));
    Ld_01(index)=nan;Ld_02(index)=nan;Ld_03(index)=nan;Ld_04(index)=nan;Ld_05(index)=nan;
    Ld = [Ld_01,Ld_02,Ld_03,Ld_04,Ld_05];
end

if date_num == 6 
    ld_path = 'D:\0Mypaper\Paper5_SW-TES\ASTER\1.majorRevison_20240720\EAR5\pixel_by_pixel\kailash_0320_2022\';
    Ld_01 = double(readgeoraster([ld_path,'Ld_kailash_0320_2022_sa_01.tif']));
    Ld_02 = double(readgeoraster([ld_path,'Ld_kailash_0320_2022_sa_02.tif']));
    Ld_03 = double(readgeoraster([ld_path,'Ld_kailash_0320_2022_sa_03.tif']));
    Ld_04 = double(readgeoraster([ld_path,'Ld_kailash_0320_2022_sa_04.tif']));
    Ld_05 = double(readgeoraster([ld_path,'Ld_kailash_0320_2022_sa_05.tif']));
    Ld_01(index)=nan;Ld_02(index)=nan;Ld_03(index)=nan;Ld_04(index)=nan;Ld_05(index)=nan;
    Ld = [Ld_01,Ld_02,Ld_03,Ld_04,Ld_05];
end

%% Step1 未考虑T-A效应的TES算法实现 
[LST,LSE_01,LSE_02,LSE_03,LSE_04,LSE_05,e_max_list] = ...
    Function_TES_no_TAeffect_satllite_data(tb_01,tb_02,tb_03,...
            tb_04,tb_05,ewl,Ld_01,Ld_02,Ld_03,Ld_04,Ld_05,...
            LUT_2D,beta_flat,bands_number);
AA = max(max(LST));
LST(LST==0)=nan;LSE_01(LSE_01==0)=nan;LSE_02(LSE_02==0)=nan;LSE_03(LSE_03==0)=nan;
LSE_04(LSE_04==0)=nan;LSE_05(LSE_05==0)=nan;e_max_list(e_max_list==0)=nan;
%% Step2 考虑T-A效应的TES算法实现
[MLST,MLSE_01,MLSE_02,MLSE_03,MLSE_04,MLSE_05] = ...
    Function_TES_yes_TAeffect_satllite_data(SSP,SVF,tb_01,tb_02,tb_03,tb_04,tb_05,...
            LST,LSE_01,LSE_02,LSE_03,LSE_04,LSE_05,ewl,Ld_01,Ld_02,Ld_03,Ld_04,Ld_05,...
            LUT_3D,beta_ssp,e_max_list,bands_number);
MLST(MLST==0)=nan;MLSE_01(MLSE_01==0)=nan;MLSE_02(MLSE_02==0)=nan;MLSE_03(MLSE_03==0)=nan;
MLSE_04(MLSE_04==0)=nan;MLSE_05(MLSE_05==0)=nan;e_max_list(e_max_list==0)=nan;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: He Z.-W

%%Data1: Regarding the field description of simulation_TbTg.mat
Note: the simudata_SSP_00 to simudata_SSP_10 represent all the simulated data, and the rest are separated from them for the convenience of operation.
'SSP': small-scale self-heating parameter	
'SVF': sky-view factor	
'TS(K)' : the bottom layer temperature of the 98 atmosphere profiles
'WVC': the water vapor content of the 98 atmosphere profiles	
'emiss_01': LSE of ASTER TIR band 10 obtained using the selected 51 emissivity spectra
'emiss_02': LSE of ASTER TIR band 11 obtained using the selected 51 emissivity spectra
'emiss_03': LSE of ASTER TIR band 12 obtained using the selected 51 emissivity spectra
'emiss_04': LSE of ASTER TIR band 13 obtained using the selected 51 emissivity spectra
'emiss_05': LSE of ASTER TIR band 14 obtained using the selected 51 emissivity spectra
'BT_01': Simulated top-of-atmosphere brightness temperature of ASTER TIR band 10  using the proposed 3D TIR radiative transfer model
'BT_02': Simulated top-of-atmosphere brightness temperature of ASTER TIR band 11  using the proposed 3D TIR radiative transfer model	
'BT_03': Simulated top-of-atmosphere brightness temperature of ASTER TIR band 12  using the proposed 3D TIR radiative transfer model	
'BT_04': Simulated top-of-atmosphere brightness temperature of ASTER TIR band 13  using the proposed 3D TIR radiative transfer model	
'BT_05': Simulated top-of-atmosphere brightness temperature of ASTER TIR band 14  using the proposed 3D TIR radiative transfer model
'Tg_01': Simulated ground brightness temperature of ASTER TIR band 10  using the proposed 3D TIR radiative transfer model	
'Tg_02': Simulated ground brightness temperature of ASTER TIR band 11  using the proposed 3D TIR radiative transfer model	
'Tg_03': Simulated ground brightness temperature of ASTER TIR band 12  using the proposed 3D TIR radiative transfer model	
'Tg_04': Simulated ground brightness temperature of ASTER TIR band 13  using the proposed 3D TIR radiative transfer model	
'Tg_05': Simulated ground brightness temperature of ASTER TIR band 14  using the proposed 3D TIR radiative transfer model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Data2: Regarding the field description of LUT_3D.mat
The value in the first column divided by 10 indicates SSP
The value in the second column divided by 10 indicates SVF
The third column indicates that this coefficient is used to estimate the ground BT of the corresponding bands. 1,2,3,4, and 5 represent the ASTER TIR bands 10, 11, 12, 13, and 14 respectively
The fourth column represents band combinations, for example, 13 indicates the combination of band 10 and band 12
Columns 5 to 9 represent the five coefficients of ISW
Column 10 represents the RMSE when using simulated data for fitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Data3: Regarding the field description of LUT_2D.mat
The first column represents band combinations, for example, 13 indicates the combination of band 10 and band 12
Columns 2 to 6 represent the five coefficients of ISW
Column 7 represents the RMSE when using simulated data for fitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Data4: description of beta_ssp.mat
The three numbers inside are respectively Three regression coefficients and RMSE of the εmin-MMD relationship for different SSPs over mountainous
areas for five ASTER TIR bands.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Data5: description of beta_flat.mat
The three numbers inside are respectively three regression coefficients of the εmin-MMD relationship

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




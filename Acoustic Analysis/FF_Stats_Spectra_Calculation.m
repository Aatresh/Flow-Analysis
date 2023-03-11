%% CODE CALCULATES FAR FIELD ACOUSTIC STATISTICS & SPECTRA FROM PRESSURE FILES
clc;   clearvars;   set(0,'defaultfigurecolor',[1 1 1]);   code = 'FF_Comp'; 
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\'; 
cd([root1 'Jet_Analysis\Global_Functions\']); addpath([root1 'Jet_Analysis\Far_Field_Acoustic_codes\']);
% Test Conditions & nozzle configuration
tests = {'NPR_2p5_TR_1p0';'NPR_2p6_TR_1p0';'NPR_2p9_TR_1p0';'NPR_3p0_TR_1p0';'NPR_3p6_TR_1p0';'NPR_4p0_TR_1p0';'NPR_4p5_TR_1p0';...
         'NPR_2p5_TR_3p0';'NPR_2p9_TR_3p0';'NPR_3p0_TR_3p0';'NPR_3p6_TR_3p0';'NPR_4p0_TR_3p0';'NPR_4p5_TR_3p0'};

condition = [tests(1)];   config = 'TS';   nozzle = 'Major';     

% Input & Output Drive Selection 
[OutputStruct] = GF_DriveSelect(config,nozzle,code);   nozzle = OutputStruct.nozzle; 
% No. of mics;           Refernce pressure(Pa)
micNos  = 16;             Pref  = 20e-6;
% Acquistion frequency;  FFT Block Size
fs      = 204800;         N     = 4096;
% Frequency Spectrum based on Block size
freq = fs/N:fs/N:fs;        
%% Computation Loop
for n = 1:length(condition)
  driveIn  = [OutputStruct.in_root condition{n}(9:14) '\' condition{n}(1:7) '\'];
  driveOut = [OutputStruct.out_root];
% Pressure RMS OASPL;                    Time Varying Pressure Derivative RMS 
  pres_rms_mat = zeros(micNos);          dpdt_rms_mat = zeros(micNos);
% Pressure Skewness;                     Time Varying Pressure Derivative Skewness
  pres_skew_mat = zeros(micNos);         dpdt_skew_mat = zeros(micNos);
% Pressure Kurtosis;                     Time Varying Pressure Derivative Kurtosis
  pres_kurt_mat = zeros(micNos);         dpdt_kurt_mat = zeros(micNos);
% Pressure SPL Result Matrix;            Pressure Power Spectral Density Matrix           
  PHI_yy_amp_avg = zeros(micNos,N);      PSD = zeros(micNos,N);
% Input File name;                       Input Pressure Channels
  fileName = 'Pressure_FF_Pos001.dat';   chanVals = load([driveIn fileName]);     
 % SPL Computation Loop
  for ctr = 1:micNos
  % Raw pressure Signal;                 Pressure Signal Fluctuation Component
    signalRaw = chanVals(ctr,:);         signalFluc = signalRaw-mean(signalRaw);
  % Pressure Statictics Computation
    [pres_skew_mat(ctr),dpdt_skew_mat(ctr),pres_rms_mat(ctr),pres_kurt_mat(ctr),dpdt_rms_mat(ctr),dpdt_kurt_mat(ctr)] = stats_calculation(signalFluc,1/fs);
  % Pressure SPL & PSD Computation
    [PHI_yy_amp_avg(ctr,:),PSD(ctr,:)] = FFT_Code_Function_Nick(N,fs,signalFluc);
  end            
%% SAVING DATA
  pres_skew = pres_skew_mat(:,:);       dpdt_skew = dpdt_skew_mat(:,:);
  pres_rms  = pres_rms_mat(:,:);        pres_kurt = pres_kurt_mat(:,:);
  dpdt_rms  = dpdt_rms_mat(:,:);        dpdt_kurt = dpdt_kurt_mat(:,:);
  pres_spl  = 20*log10(pres_rms/Pref);
% File Names - Pressure Statictics               
  fileOutSpl      = ['spl_pres_N'  condition{n}(9:14) '_' condition{n}(1:7)];
  fileOutPresSkew = ['skew_pres_N' condition{n}(9:14) '_' condition{n}(1:7)];
  fileOutDpdtSkew = ['skew_dpdt_N' condition{n}(9:14) '_' condition{n}(1:7)];
  fileOutPresKurt = ['kurt_pres_N' condition{n}(9:14) '_' condition{n}(1:7)];
  fileOutDpdtKurt = ['kurt_dpdt_N' condition{n}(9:14) '_' condition{n}(1:7)];
% Saving to Disk - Pressure Statictics            
  save([driveOut fileOutSpl],'pres_spl');       save([driveOut fileOutPresSkew],'pres_skew');
  save([driveOut fileOutDpdtSkew],'dpdt_skew'); save([driveOut fileOutPresKurt],'pres_kurt');
  save([driveOut fileOutDpdtKurt],'dpdt_kurt');
% File Name - Pressure SPL & PSD               
  fileNameSPL=['FFT_Nick_nofilter_N' condition{n}(9:14) '_' condition{n}(1:7)];
  fileNamePSD=['PSD_Function_N' condition{n}(9:14) '_' condition{n}(1:7)];
% File Name - Pressure SPL & PSD   
  save([driveOut fileNameSPL],'freq','PHI_yy_amp_avg');
  save([driveOut fileNamePSD],'PSD');
% Plotting mics 16 & 1 for reference
  semilogx(freq(1:N/2),PHI_yy_amp_avg(16,1:N/2),'k','LineWidth',1.2,'DisplayName','Mic - 16'); grid on; hold on;      box off;
  semilogx(freq(1:N/2),PHI_yy_amp_avg(13,1:N/2),'b','LineWidth',1.2,'DisplayName','Mic - 13'); aX = gca;
  semilogx(freq(1:N/2),PHI_yy_amp_avg(1,1:N/2) ,'r','LineWidth',1.2,'DisplayName','Mic - 1');  xlim([300 100000]);
  aX.TickLabelInterpreter = 'latex';           xlabel('Hz'); ylabel('SPL(dB)');                cL = legend;           cL.Interpreter = 'latex';              
  aX.XLabel.Interpreter   = 'latex';           aX.YLabel.Interpreter   = 'latex';              cL.EdgeColor = [1 1 1];
end

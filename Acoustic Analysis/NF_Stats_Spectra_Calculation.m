%% CALCULATES THE OASPL, SPECTRA & HIGHER ORDER METRICS FOR NEAR FIELD DATA
clc; clearvars ; root = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\Near_Field_Data';
calculation = [1,2];   % 1=stats      2=spectra
tests = {'NPR_2p5_TR_1p0', 'NPR_2p9_TR_1p0','NPR_3p0_TR_1p0', 'NPR_3p6_TR_1p0','NPR_4p0_TR_1p0'};

condition = tests(4);   nozzle = 'Minor';    jet_typ = 'T2';    

fs = 204800;                N = 4096;           Pref = 20e-6;
freq = fs/N:fs/N:fs;        dt = 1/fs;          mics = 16;      
switch jet_typ
  case 'S'
   drive  = [root '\Single_Jet\' nozzle '\'];       postn = 145;
   output = [root 'Single_Jet\Results\' nozzle '\'];
 case 'T'   
   drive  = [root '\Twin_Jet\' nozzle '\'];         postn = 145;
   output = [root '\Twin_Jet\Results\' nozzle '\']; 
 case 'T2'   
   drive  = [root '\Twin_Jet_V2\' nozzle '\'];      postn = 63;
   output = [root '\Twin_Jet_V2\Results\' nozzle '\'];
end
tic  
for j = 1:length(condition)
  drive_in = [drive condition{j}(9:14) '\' condition{j}(1:7) '\'];
  drive_out= [output condition{j}(9:14) '\' condition{j}(1:7) '\'];
%  DATA METRICS
  pres_rms_mat = zeros(mics);           dpdt_rms_mat = zeros(mics);
  pres_skew_mat = zeros(mics);          dpdt_skew_mat = zeros(mics);
  pres_kurt_mat = zeros(mics);          dpdt_kurt_mat = zeros(mics);
  PHI_yy_amp_avg = zeros(mics,N);       PSD = zeros(mics,N);
  for x = 1:postn
    chanvals = load(sprintf('%sPressure_NF_Pos%0.3d.dat',drive_in,x));
    for ch = 1:size(chanvals,1)
      sig_analyze = chanvals(ch,:);
      sig_analyze2 = sig_analyze-mean(sig_analyze);
      for calc = calculation
        if calc == 1
           [pres_skew_mat(ch),dpdt_skew_mat(ch),pres_rms_mat(ch),pres_kurt_mat(ch),dpdt_rms_mat(ch),dpdt_kurt_mat(ch)]=stats_calculation(sig_analyze2,dt);
        elseif calc == 2
           [PHI_yy_amp_avg(ch,:),PSD(ch,:)] = FFT_Code_Function_Nick(N,fs,sig_analyze2);
%          [PHI_yy_amp_avg(ch,:)] = FFT_Code_Function_Nick(N,fs,sig_analyze2);
        end
      end
    end
   disp([newline '>> Saving data for Position - ',num2str(x)]);
   for calc=calculation       
     if calc == 1
        pres_skew=pres_skew_mat(:,:);           dpdt_skew=dpdt_skew_mat(:,:);
        pres_rms=pres_rms_mat(:,:);             pres_kurt=pres_kurt_mat(:,:);
        dpdt_rms=dpdt_rms_mat(:,:);             dpdt_kurt=dpdt_kurt_mat(:,:);
        pres_spl=20*log10(pres_rms/Pref);
%  ASSIGNING FILE NAMES           
        file_out_p =         ['spl_pres_' condition{j}(9:14) '_' condition{j}(1:7) '_' sprintf('%0.3d',x)];
        file_out_pres_skew = ['skew_pres_' condition{j}(9:14) '_' condition{j}(1:7) '_' sprintf('%0.3d',x)];
        file_out_dpdt_skew = ['skew_dpdt_' condition{j}(9:14) '_' condition{j}(1:7) '_' sprintf('%0.3d',x)];
        file_out_pres_kurt = ['kurt_pres_' condition{j}(9:14) '_' condition{j}(1:7) '_' sprintf('%0.3d',x)];
        file_out_dpdt_kurt = ['kurt_dpdt_' condition{j}(9:14) '_' condition{j}(1:7) '_' sprintf('%0.3d',x)];
%  SAVING OUTPUT            
        save([drive_out file_out_p],'pres_spl');
        save([drive_out file_out_pres_skew],'pres_skew');
        save([drive_out file_out_dpdt_skew],'dpdt_skew');
        save([drive_out file_out_pres_kurt],'pres_kurt');
        save([drive_out file_out_dpdt_kurt],'dpdt_kurt');
     elseif calc == 2
        file_name_out=['FFT_Nick_nofilter_' condition{j}(9:14) '_' condition{j}(1:7) '_' sprintf('%0.3d',x)];
        save([drive_out file_name_out],'freq','PHI_yy_amp_avg');
     end
   end
  end
end
toc
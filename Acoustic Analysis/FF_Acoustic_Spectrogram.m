%% CODE TO COMPUTE THE SPECTROGRAM FOR FAR FIELD ACOUSTIC RESULTS
clc; clearvars; fclose all; set(0,'defaultfigurecolor',[1 1 1]);
Pref = 20e-6;       fs = 204800;        N = 4096;           chc = 'y';
angles = [152 148 144 140 136 132 126 120 116 110 105 100 90 70 60 45];
root = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';  cd([root 'Jet_Analysis\Global_Functions\']);
addpath([root 'Jet_Analysis\Far_Field_Acouctic_Codes\']);  
tests = {'NPR_2p5_TR_1p0';'NPR_2p9_TR_1p0';'NPR_3p0_TR_1p0';'NPR_3p6_TR_1p0';'NPR_4p0_TR_1p0';'NPR_5p0_TR_1p0'};

condition = [tests(1)];    nozzle = 'Minor';      jet_typ = 'T';
switch jet_typ
   case 'S'
      drive  = [root 'Far_Field_Data\Single_Rectan\' nozzle '\'];  Deq = 0.0206502;
      output = [root 'Far_Field_Data\Single_Rectan\Results\' axis '\'];
   case 'T'
      drive  = [root 'Far_Field_Data\Twin_Rectan\' nozzle '\'];    Deq = 0.01847;
      output = [root 'Far_Field_Data\Twin_Rectan\Results\' nozzle '\'];
end

for n = 1:length(condition)
  [Mj,Uj,~,~] = GF_Velocity(condition{n});
  drive_in = [drive condition{n}(9:14) '\' condition{n}(1:7) '\'];
  chanvals = load([drive_in 'Pressure_FF_Pos001.dat']);
  while strcmp(chc,'y')
    prompt = [newline '>> SELECT MIC NO (1-16) - ']; mic = input(prompt);
    if mic < 1 || mic > 16
       disp('>> MIC NOT FOUND');  break;
    end
    sig = chanvals(mic,:);      sig = sig - mean(sig);
%  COMPUTING SPECTROGRAM USING SHORT TIME FFT
    [res_fft,freq,time,~] = spectrogram(sig,N,[],[],fs);
%  COMPUTING POWER SPECTRAL DENSITY & SPL FROM FFT RESULT
    PSD = abs(res_fft).^2/(fs*N);       PSD(2:end,:) = 2*PSD(2:end,:);
%  CONVERT TO SPL
    SPL = 10*log10(PSD*freq(2)/(20e-6)^2);      St = freq*Deq/Uj;
%  PLOTTING SPECTROGRAM - Hz
    figure; pcolor(time,freq,SPL); shading interp; colormap jet;  ylim([400 90000]);  c = colorbar;
    ax = gca; ax.YScale = 'log'; xlabel('\bfTime(s)'); ylabel('$Frequency(Hz)$','Interpreter','latex');
    title(c,'$dB/Hz$','Interpreter','latex');  set(gca,'FontSize',12,'FontName','times');
    title(['$Spectrogram:\psi_o = {',num2str(angles(mic)),'}^o$'],'Interpreter','latex');
%  PLOTTING SPECTROGRAM - St
    figure; pcolor(time,St,SPL); shading interp; colormap jet;  ylim([0.014 4.9]);  c = colorbar;
    ax = gca; ax.YScale = 'log'; xlabel('\bfTime(s)'); ylabel('$Frequency(St)$','Interpreter','latex');
    title(c,'$dB/St$','Interpreter','latex');  set(gca,'FontSize',12,'FontName','times');
    title(['$Spectrogram:\psi_o = {',num2str(angles(mic)),'}^o$'],'Interpreter','latex');
    prompt = [newline '>> PLOT ANOTHER MIC(y/n) - '];  chc = input(prompt,'s');
  end
end
    
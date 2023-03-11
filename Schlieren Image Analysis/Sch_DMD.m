%% CODE TO EXTRACT DMD MODE ENERGIES FROM SCHLIEREN IMAGES 
%  ACCEPTS IMAGE MATRCI & COMPUTES DMD ENERGY PROFILES
clear; clc; close all; fclose all; set(0,'defaultfigurecolor',[1 1 1]);  code = 'SchDMD'; 
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\'; 
cd([root1 'Jet_Analysis\Global_Functions\']); addpath([root1 'Jet_Analysis\DMD_Codes\']); load('dmd_Colormap.mat');      
tests = {'NPR_2p5_TR_1p0', 'NPR_2p6_TR_1p0', 'NPR_2p9_TR_1p0','NPR_3p0_TR_1p0',...
         'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_5p0_TR_1p0', 'NPR_6p0_TR_1p0'};

condition = tests(2);      config = 'S';     nozzle = 'Minor';    NF = 'H';   

%%  DRIVE SELECTION & DATA ORGANIZATION
[OutputStruct] = GF_DriveSelect(config,nozzle,code);
in_root = OutputStruct.in_root; out_root = OutputStruct.out_root;  nozzle = OutputStruct.nozzle; dt = OutputStruct.dt;
%  Jet parameters
[Mj,Uj,NPR,NTR] = GF_Velocity(condition{1});
disp([newline '>> ACQUISITION RATE - ',num2str(1/dt),' Hz']);
%  Loading Image matrix
driveIn  = [in_root condition{1}(9:14) '\' condition{1}(1:7)  '\'];
driveOut = [out_root condition{1}(9:14) '\' condition{1}(1:7) '\'];  
load([driveIn condition{1}(1:7) '_DAT']);   load([driveIn 'X']);  load([driveIn 'Y']);
%  Background subtraction
if isempty(OutputStruct.bckgrnd_root) ~= 1  
   disp([newline '->> Subtracting Background' newline]);
   Bckgrnd = load([OutputStruct.bckgrnd_root 'Bckgrnd_DAT']);     cLim1 = [0 1.5];
   Bckgrnd = Bckgrnd.Master_U;                                    %avgTwo = mean(Master_U,3);  
   M2 = Master_U./mean(Bckgrnd,3);   %avgOne = mean(M2,3);         M2Fluc = M2 - avgOne;
else
   Bckgrnd = [];  cLim1 = [0 500];   avgOne = mean(Master_U,3);   avgTwo = avgOne;                               %#ok
end;   clear Bckgrnd Master_U  avgTwo;
%  Axis Normalization & labelling
[Xn,Yn,~,~,xName,yName,lenScales,figSiz] = GF_AxisDefnSch(config,nozzle,NF,X,Y);   Xn = Xn';   Yn = Yn';
%  Cropping images based on normalization factor & computing fluctuating component
[M2,Xn,Yn,limX,limY] = imgCrop(M2,Xn,Yn,NF);   avgOne = mean(M2,3);   M2Fluc = M2 - avgOne;  
%  DEFINING NO. OF MODES & RESHAPING MATRIX
matDims = size(M2Fluc);   M2Fluc = reshape(M2Fluc,matDims(1)*matDims(2),matDims(3));      disp([newline '----> Done <----']);
%% DMD COMPUTATION FUNCTION
tic;    nModes = size(M2Fluc,2)-1;    [dmdOutput] = SchFn_DMDCalc(M2Fluc,nModes,dt,matDims);   toc;   
dmdModes = dmdOutput.dmdModes;   dmdFreq = dmdOutput.dmdFreq;
disp([newline '----> DMD Computation Complete <----']);    textSize = 12;
%  Find dominant oscillatory mode & define St number
domModeNo = find(dmdOutput.dmdModeAmp==max(dmdOutput.dmdModeAmp));   dmdSt = dmdFreq*lenScales(1)/Uj;   stLim = (1/(2*dt))*lenScales(1)/Uj; 
%  Display dominant mode number & Frewuency(Hz,St)
disp([newline '->> Frequency Resolution  - ',num2str(dmdOutput.dmdFreq(2)-dmdOutput.dmdFreq(1))]);
disp([newline '>> Dominant Mode Number   - ', num2str(domModeNo)]);
disp([newline '>> Dominant Frequency(Hz) - ', num2str(dmdOutput.dmdFreq(domModeNo))]);
disp([newline '>> Dominant Strouhal No.  - ', num2str(dmdSt(domModeNo))]);
%% PLOTTING MODE AMPLITUDES, GROWTH RATES & ENERGY DISTRIBUTION
%  Mode Energy amplitudes: Frequency
figure(1);   stem(dmdFreq,dmdOutput.dmdModeAmp,'k:o', 'filled','MarkerSize',6,'MarkerFaceColor','k');   hold on;   
plot(dmdFreq(domModeNo),dmdOutput.dmdModeAmp(domModeNo),'co','MarkerSize',8,'MarkerFaceColor','c');     ax = gca;
xlabel('$Frequency(Hz)$','Interpreter','latex');   xlim([300 (1/(2*dt))]);   grid on;   box off;   ax.FontSize = textSize;      
ylabel('$Amplitude(\xi)$','Interpreter','latex');  ylim([0 1.5]);  ax.XScale = 'log';   ax.YScale = 'log';  
title('Mode amplitude vs Frequency','Interpreter','latex');   ax.TickLabelInterpreter = 'latex';   set(gcf,'Position',[20 300 660 410]);
%  Mode Energy amplitudes: Strouhal number
figure(2);   stem(dmdSt,dmdOutput.dmdModeAmp,'k:o', 'filled','MarkerSize',6,'MarkerFaceColor','k');   hold on;   
plot(dmdSt(domModeNo),dmdOutput.dmdModeAmp(domModeNo),'co','MarkerSize',8,'MarkerFaceColor','c');     ax = gca;
xlabel('$St$','Interpreter','latex');    xlim([0.01 stLim]);   grid on;   box off;   ax.TickLabelInterpreter = 'latex';   
ylabel('$Amplitude(\xi)$','Interpreter','latex');  ylim([0 1.5]);  ax.XScale = 'log';   ax.YScale = 'log';   
title('Mode amplitude vs St','Interpreter','latex');   ax.FontSize = textSize;   set(gcf,'Position',[685 300 660 410]);
%  Growth Rates
figure(3);   stem(dmdFreq,dmdOutput.dmdGrwthRate,'k:o', 'filled','MarkerSize',6,'MarkerFaceColor','k');   hold on; 
plot(dmdFreq(domModeNo),dmdOutput.dmdGrwthRate(domModeNo),'co','MarkerSize',8,'MarkerFaceColor','c');     ax = gca;
xlabel('$Frequency(Hz)$','Interpreter','latex');        xlim([100 (1/(2*dt))]);   grid on;   box off;     ax.XScale = 'log';  
ylabel('$Growth Rate(\Omega)$','Interpreter','latex');  ylim([-1.05 1.05]);     ax.TickLabelInterpreter = 'latex';
title('Growth Rate vs Frequency','Interpreter','latex');   ax.FontSize = textSize;   set(gcf,'Position',[520 90 850 425]);
%  Eigen Value Distribution
figure(4);   plot(dmdOutput.allLambda,'ko','LineWidth',1.2,'MarkerFaceColor','k');   grid on;   box off;   ylim([-1.02 1.02]);  ax = gca;   axis equal;   
xlabel('$Modes(n)$','Interpreter','latex');   ylabel('$Time \thinspace Eigens(\lambda_t)$','Interpreter','latex');   ax.TickLabelInterpreter = 'latex';   
 ax.FontSize = textSize;   set(gcf,'Position',[20 90 500 425]);
%  Energy Distribution vs Frequency
figure(5);   stem(dmdFreq,dmdOutput.dmdModeAmp,'r:o','LineWidth',1.2,'MarkerFaceColor','k','MarkerSize',6);   grid on;   box off;   ax = gca;   ax.FontSize = textSize;      
ylim([0 1.1]);   xlabel('$Frequency(Hz)$','Interpreter','latex');   ylabel('$Energy(\epsilon_n)$','Interpreter','latex');   ax.XScale = 'log';   ax.YScale = 'log';    
xlim([100 (1/(2*dt))]);   title('Energy vs Frequency','Interpreter','latex');   set(gcf,'Position',[20 200 660 410]);   ax.TickLabelInterpreter = 'latex';

%% PLOTTING MODE SHAPES OF DOMINANT FREQUENCY
%  Plotting energy distribution for the dominant mode number - Real part
nCtr = 4;   T = 1/dmdFreq(domModeNo);   nt = 10;   tSteps = linspace(0,T,nt);   domMode = dmdModes(:,:,domModeNo)*exp(2i*pi*dmdFreq(domModeNo)*tSteps(nCtr));   
figure(6);   pcolor(Xn,Yn,real(domMode));   shading interp;   colormap(dmd_Colormap);   caxis([-0.004 0.004]);   axis equal;   xlim(limX);   ylim(limY);   
xlabel(xName,'Interpreter','latex');   ylabel(yName,'Interpreter','latex');    ax = gca;   c = colorbar;   set(gcf,"Position",[120 485 690 330]);
ax.FontSize = textSize;   ax.TickLabelInterpreter = 'latex';   c.TickLabelInterpreter = 'latex';   cT = title(c,'$R_e(\epsilon_f)$','Interpreter','latex');    
title(['Dominant Mode \it (f = ',num2str(dmdOutput.dmdFreq(domModeNo)),'Hz)'],'Interpreter','latex');  cT.Position = [50 82 0];    cT.FontSize = 13;
%  Plotting energy distribution for the dominant mode number - absolute part(nomalized with 80% of peak magnitude for better visibility of modes)
peakMag = max(max(abs(domMode)));   peakMag = peakMag*0.9;
figure(7);  pcolor(Xn,Yn,abs(domMode)/peakMag);   shading interp;   colormap(dmd_Colormap);    c = colorbar;   caxis([0 0.8]);   axis equal;   ax = gca;   
xlim(limX);    ylim(limY);   xlabel(xName,'Interpreter','latex');   ylabel(yName,'Interpreter','latex');     set(gcf,"Position",[120 75 690 330]);   %ax.ColorScale = 'log';
ax.FontSize = textSize;   ax.TickLabelInterpreter = 'latex';   c.TickLabelInterpreter = 'latex';   cT = title(c,'$|\epsilon_f|$','Interpreter','latex');    
title(['Dominant Mode \it (f = ',num2str(dmdOutput.dmdFreq(domModeNo)),'Hz)'],'Interpreter','latex');   cT.Position = [50 82 0];   cT.FontSize = 13;
%% SPATIAL FOURIER ANALYSIS
inputStruct.Xn = Xn;        inputStruct.inpMat = (dmdModes(:,:,domModeNo));              
inputStruct.Yn = Yn;        inputStruct.normFac = lenScales(1)*1000;      
[outputStruct] = GF_SpatialFourier(inputStruct);
%% PLOTTING SPECIFIC FREQUENCIES
prompt = [newline '>> Enter Frequency - '];   freqVal = input(prompt);   srtdFreq = sort(dmdFreq);
[a1,~] = find(srtdFreq<freqVal);  freqVal = srtdFreq(a1(end)+1);   [modeNo,~] = find(dmdFreq==freqVal);   theta = linspace(0,4*pi,100);
%  Mode Energy (amplitude fields)
figure;   pcolor(Xn,Yn,squeeze(abs(dmdModes(:,:,modeNo))));   shading interp;   colormap(dmd_Colormap);    c = colorbar;   caxis([0.0005 0.008]);  axis equal;
set(gca,'ColorScale','log');   c.Ticks = [0.0005 0.0008 0.002 0.005 0.008];    xlim(limX);    ylim(limY);  ax = gca;   set(gcf,"Position",[120 485 690 330]);  
xlabel(xName,'Interpreter','latex');   ylabel(yName,'Interpreter','latex');    ax.FontSize = textSize;   ax.TickLabelInterpreter = 'latex';
c.TickLabelInterpreter = 'latex';      title(c,'$\epsilon_f$','Interpreter','latex');    title(['Mode Amplitude \it (f=',num2str(dmdFreq(modeNo)),'Hz)'],'Interpreter','latex');
%  Mode Phase
figure;   pcolor(Xn,Yn,squeeze(angle(dmdModes(:,:,modeNo))));   shading interp;   c = colorbar;   caxis([-pi pi]);   colormap bluewhitered;
xlabel(xName,'Interpreter','latex');   axis equal;   xlim(limX);    ylim(limY);   ax = gca;   set(gcf,"Position",[120 75 690 330]);   c.Ticks = [-pi -pi/2 0 pi/2 pi];    
ylabel(yName,'Interpreter','latex');   ax.FontSize = textSize;   ax.TickLabelInterpreter = 'latex';   c.TickLabelInterpreter = 'latex';
title(c,'$\Delta \phi$','Interpreter','latex');   title(['Mode Phase \it (f=',num2str(dmdFreq(modeNo)),'Hz)'],'Interpreter','latex');    c.TickLabels = {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'};  
%% .gif Generation
gifName  = ['DMD-ModeFreq - ', num2str(abs(dmdFreq(modeNo)))];   
gifTitle = ['Mode Energy-Phase Transition(\it f=',num2str(dmdFreq(modeNo)),'Hz)'];
%  ADDING PARAMETERS TO INPUT STRUCT
DMD_INP.mode = modeNo;    DMD_INP.theta = theta;   DMD_INP.Xn = Xn;        DMD_INP.dmd_Colormap = dmd_Colormap;         
DMD_INP.lim_x = limX;     DMD_INP.lim_y = limY;    DMD_INP.xname = xName;  DMD_INP.yname = yName;
DMD_INP.nozzle = nozzle;  DMD_INP.drive_out = driveOut;                    DMD_INP.colormap = dmd_Colormap; 
DMD_INP.gifName = gifName;DMD_INP.NF = NF;         DMD_INP.TYP = 'DMD';    DMD_INP.Yn = Yn;
DMD_INP.gifTitle = gifTitle;                       dmdInput = permute(dmdModes,[3 1 2]);
%  GIF DISPLAY & SAVE 
GF_GIF_DisplayWrite(dmdInput,DMD_INP);
%% ------------------------------------ Functions ---------------------------------------------------------------%%
%%      Function 1: CROP IMAGE BASED ON NORMALIZATION FACTOR
function [M2,Xn,Yn,limX,limY] = imgCrop(M2,Xn,Yn,NF)
%  Defining radia limits based on normalization factor
    if strcmp(NF,'D')
       yLim = 2.5;   limX = [0 6];   limY = [-1.5 1.5];
    else
       yLim = 3.5;   limX = [0 8];   limY = [-2.5 2.5];
    end;                            [~,zLoc] = find(Xn==0);       
% Re-sizing axes & image matrix according to limits
    [upLimit,~] = find(Yn<yLim);   [lowLimit,~] = find(Yn>(-1*yLim));   
    upLimit = upLimit(end)+1;      lowLimit = lowLimit(1)-1;
    Yn = Yn(lowLimit:upLimit);      Xn = Xn(zLoc:end);
    M2 = M2(lowLimit:upLimit,zLoc:end,:);
end

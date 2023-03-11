%% Code to visualize the staging phenomenon occuring at specific operating conditions
clc;    clearvars;  fclose all; set(0,'defaultfigurecolor',[1 1 1]); tic;
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';
addpath([root1 'Jet_Analysis\Exit_Test_Codes\']); cd([root1 'Jet_Analysis\Global_Functions\']);   
%% Condition selection & Input path definition 
root = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\Exit_Tests';
tests = {'NPR_2p5_TR_1p0';'NPR_2p6_TR_1p0';'NPR_2p7_TR_1p0';'NPR_2p8_TR_1p0';'NPR_2p9_TR_1p0';'NPR_3p0_TR_1p0';...
         'NPR_3p1_TR_1p0';'NPR_3p2_TR_1p0';'NPR_3p3_TR_1p0';'NPR_3p4_TR_1p0';'NPR_3p5_TR_1p0';'NPR_3p6_TR_1p0'; ...
         'NPR_4p0_TR_1p0';'NPR_5p0_TR_1p0'; 'NPR_6p0_TR_1p0'};        load('blckToWhite.mat');

condition = [tests(1:12)];   nozzle = 'S';   config = 'CF-3';   run = 'R2';   Peaks = 'on';  txtSiz = 13;

%%  Condition parameters
if strcmp(config,'CF-3')
   mics = 4;
else
   mics = 3;
end
%  Input path definition
if strcmp(nozzle,'S')
   driveIn  = [root '\Single_Rectan\' config '\' run '\'];  driveImg  = [root '\Single_Rectan\' config '\']; 
   driveRes = [root '\Single_Rectan\' config '\Results\' run '\'];                           cLim = [81 141];
   driveOut = [root '\Single_Rectan\Images\' config '\'];          nozName = 'Rectangular';  cLim2 = [85 130];
else
   driveIn  = [root '\Circular_Med\' config '\' run '\'];   driveImg  = [root '\Circular_Med\' config '\']; 
   driveRes = [root '\Circular_Med\' config '\Results\' run '\'];                            cLim = [85 120];
   driveOut = [root '\Circular_Med\Images\' config '\'];           nozName = 'Circular';     cLim2 = [85 130];    
end
%  Acquisition frequency, Reference pressure, Equivalent/exit diameter, recording time 
fs = 204800;   Pref = 20e-6;   dt = 1/fs;   Deq = 0.0206502;   time = 0:dt:2;   
%  SPL & Strouhal number
splMat = zeros(mics,4096,length(condition));   StMat  = zeros(4096,length(condition));
%  Jet Mach number, velocity & NPR
jetMach = zeros(length(condition),1);     jetVel = zeros(length(condition),1);     jetNPR  = zeros(length(condition),1);
%  Nozzle & mic layout
figure(1);   imshow([driveImg config '.jpg']);   set(gcf,'Position',[999 250 475 440]); 
%  Color pallete
lineColor = GF_CustomColormap([0 0.5 1],[1 0 0;0.5 0.1 0.5;0.2 0.2 1],length(condition));
%  Data loading loop
for ctr = 1:length(condition)
    file_name = ['FFT_Nick_nofilter_' condition{ctr}(9:14) '_' condition{ctr}(1:7)];
    load([driveRes  '50Hz\' condition{ctr}(9:14) '\'  file_name]);
    splMat(:,:,ctr) = PHI_yy_amp_avg;    
%  Scaling SPL of Rectangular jet to same distance as Circular jet dInit = 45mm to dFinal = 100mm
    if strcmp(nozzle,'S')
       splMat(:,:,ctr) = splMat(:,:,ctr) - 20*log10(100/45);          % Mic distance correction 
       splMat(:,:,ctr) = splMat(:,:,ctr) + 10*log10(334.9179/335.4);  % Area/Thrust Correction
       %lineColor = GF_CustomColormap([1 0.5 0],[0.1 1 0.1;0.1 0.5 0.5;0.2 0.2 1],length(condition)); 
    end
%  Flow conditions & Strouhal number matrix
    [jetMach(ctr),jetVel(ctr),jetNPR(ctr),~] = GF_Velocity(condition{ctr});   StMat(:,ctr) = (freq*Deq)/jetVel(ctr);
end;  clear ctr PHI_yy_amp_avg;   
%  Peak Frequencies, corresponding SPL & OASPL values
[peaksMat,oasplMat] = PeaksOaspl(splMat,mics,freq,nozzle);
%  NPR Legend
nprNames = num2str(jetNPR); 
%  NPR vs Mach Number (looping to plot custom color palette
for ctr = 1:length(condition)
    figure(2);   plot(jetNPR(ctr),jetMach(ctr),'h','MarkerSize',7,'MarkerFaceColor',lineColor(ctr,:),'MarkerEdgeColor',lineColor(ctr,:));   hold on;
end
title('$Test \thinspace Conditions$','Interpreter','latex');   aX = gca;   aX.XTick = [2.5 2.7 2.9 3.1 3.3 3.5 3.67];   ylim([1.2 1.51]);   box off;   grid on;   set(gcf,'Position',[330 340 670 350]); 
xlabel('$NPR$','Interpreter','latex');   ylabel('$Mach \thinspace No.$','Interpreter','latex');   set(gca,'FontSize',txtSiz);   aX.TickLabelInterpreter = 'latex';
%  Saving figure - NPR vs Mach
figName = ['NPR-vs-Mach-',nozName];   commandwindow;   GF_FigureSave(figName,driveOut,2);
%% Microphone selection & Related plot
commandwindow;   prompt = [newline,'Enter Mic Number(1 to 4) - '];   micNo = input(prompt);  
micMat = splMat(micNo,:,:);   micMat = reshape(micMat,size(splMat,2),length(condition));
%  micMrkrs = {'o','^','s','h'};     mrkr = micMrkrs{micNo};  -> Changing markers based on nozzle config
if strcmp(nozzle,'S')
   sideName = {'$Rectangular:C-D \thinspace Side';'$Rectangular:Flat \thinspace Side';'$Rectangular:C-D \thinspace Side';'$Rectangular:Flat \thinspace Side';};
   orientName = sideName{micNo};    mrkr = '^';  
else
   orientName = '$Circular \thinspace Nozzle';   mrkr = 'o';
end
%  SPL Contour (Frequency) vs  Mach No.
figure(3);   subplot(121);   pcolor(jetMach,freq(1:2048),micMat(1:2048,:));   shading interp;   figAx = gca;     figAx.YScale = 'log';
xlabel('$Mach \thinspace no.(M_j)$','Interpreter','latex');  xlim([1.223 1.5]);   c1 = colorbar;   caxis(cLim);  title(c1,'$SPL(dB)$','Interpreter','latex');  
ylabel('$Frequency(Hz)$','Interpreter','latex');   ylim([500 100000]);    figAx.TickLabelInterpreter = 'latex';  c1.TickLabelInterpreter = 'latex';
title([orientName,';Mic.no:',num2str(micNo),'$'],'Interpreter','latex');  set(gca,'FontSize',txtSiz,'FontName','times');
%  SPL Contour (Strouhal)  vs  Mach No.
subplot(122);   pcolor(jetMach,StMat(1:2048,:),micMat(1:2048,:));   shading interp;   colormap hot;   figAx = gca;   figAx.YScale = 'log';
xlabel('$Mach \thinspace no.(M_j)$','Interpreter','latex');   xlim([1.223 1.5]);      c2 = colorbar;  caxis(cLim);   figAx.TickLabelInterpreter = 'latex';
ylabel('$St$','Interpreter','latex');   ylim([0.02 4.9]);     title(c2,'$SPL(dB)$','Interpreter','latex');           c2.TickLabelInterpreter = 'latex';   
title([orientName,';Mic.no:',num2str(micNo),'$'],'Interpreter','latex');   set(gca,'FontSize',txtSiz);
set(gcf,'Position',[85 270 1260 400]);
%  Peak SPL levels vs Mach Number & NPR for instability peaks
for ctr = 1:size(peaksMat,2)
%  Mach Number
    figure(4); plot(jetMach(ctr),peaksMat(micNo,ctr,1),mrkr,'MarkerSize',7,'MarkerFaceColor',lineColor(ctr,:),'MarkerEdgeColor',lineColor(ctr,:),'DisplayName',['$NPR-',nprNames(ctr,:),'$']); hold on; grid on;
%  NPR
    figure(5); plot(jetNPR(ctr), peaksMat(micNo,ctr,1),mrkr,'MarkerSize',7,'MarkerFaceColor',lineColor(ctr,:),'MarkerEdgeColor',lineColor(ctr,:)); hold on; grid on;
end
%  Formatting
figure(4);   box off;   c3 = legend;   c3.Location = 'northeastoutside';    c3.EdgeColor = [1 1 1];   c3.Interpreter = 'latex';   set(gcf,'Position',[75 200 700 420]);   
xlabel('$Mach \thinspace no.$','Interpreter','latex');   xlim([1.2 1.51]);  set(gca,'FontSize',txtSiz);   aX = gca;   aX.TickLabelInterpreter = 'latex';
ylabel('$SPL(dB/Hz)$','Interpreter','latex');            ylim([95 140]);    title([orientName,':Mic - ',num2str(micNo),'$'],'Interpreter','latex');
figure(5);   box off; ylim([90 140]);    aX = gca;   aX.XTick = [2.5 2.7 2.9 3.1 3.3 3.5 3.67];   title([orientName,':Mic - ',num2str(micNo),'$'],'Interpreter','latex');   aX.TickLabelInterpreter = 'latex';  
xlabel('$NPR$','Interpreter','latex');   xlim([2.45 3.7]);   ylabel('$SPL(dB/Hz)$','Interpreter','latex');   set(gca,'FontSize',txtSiz);   set(gcf,'Position',[780 200 700 420]);
%  Peak frequency(St) vs Mach Number & NPR
for ctr = 1:size(peaksMat,2)
    stVal = (peaksMat(micNo,ctr,2)*Deq)/jetVel(ctr);
%  Mach number
    figure(6);   plot(jetMach(ctr),stVal,mrkr,'MarkerSize',7,'MarkerFaceColor',lineColor(ctr,:),'MarkerEdgeColor',lineColor(ctr,:),'DisplayName',['$NPR-',nprNames(ctr,:),'$']);   hold on;   grid on;
%  NPR
    figure(7);   plot(jetNPR(ctr), stVal,mrkr,'MarkerSize',7,'MarkerFaceColor',lineColor(ctr,:),'MarkerEdgeColor',lineColor(ctr,:));   hold on;   grid on;
end
figure(6);   box off;   set(gcf,'Position',[75 200 700 420]);    %c3 = legend; c3.Location = 'northeastoutside';   c3.EdgeColor = [1 1 1];   c3.Interpreter = 'latex';
xlabel('$Mach \thinspace no.$','Interpreter','latex');   xlim([1.2 1.51]);   set(gca,'FontSize',txtSiz);   textLabels(nozzle,6);   aX = gca; aX.TickLabelInterpreter = 'latex';
ylabel('$Strouhal \thinspace No.(St)$','Interpreter','latex');   ylim([0.2 0.9]);   title([orientName,':Mic - ',num2str(micNo),'$'],'Interpreter','latex');
figure(7);   box off;   xlabel('$NPR$','Interpreter','latex');   set(gca,'FontSize',txtSiz);    title([orientName,':Mic - ',num2str(micNo),'$'],'Interpreter','latex');    
ylabel('$Strouhal \thinspace No.(St)$','Interpreter','latex');   ylim([0.2 0.9]);   aX = gca;   aX.XTick = [2.5 2.7 2.9 3.1 3.3 3.5 3.67]; 
textLabels(nozzle,7);   set(gcf,'Position',[780 200 560 420]);   aX.TickLabelInterpreter = 'latex';  
%% Plotting peaks on top of SPL Contour
%  Contour
figure(8);   pcolor(jetMach,StMat(1:2048,:),micMat(1:2048,:));   shading interp;   colormap(blckToWhite);   c2 = colorbar;   caxis(cLim2);
xlabel('$Mach \thinspace no.(M_j)$','Interpreter','latex');      xlim([1.22 1.5]);   title(c2,'$SPL(dB)$','Interpreter','latex');   hold on;    
ylabel('$St$','Interpreter','latex');   ylim([0.2 1.2]);   box off;   set(gcf,'Position',[350 260 510 530]); 
aX = gca;   aX.TickLabelInterpreter = 'latex';   aX.TickLength = [0.001 0.025];   c2.TickLabelInterpreter = 'latex'; 
title([orientName,';Mic.no:',num2str(micNo),'$'],'Interpreter','latex');   set(gca,'FontSize',12);   set(gca,'Layer','top');
%  Peak values
if strcmp(Peaks,'on')
%  Max vlues for each mic spectrum
   for ctr = 1:size(peaksMat,2)
       stVal = (peaksMat(micNo,ctr,2)*Deq)/jetVel(ctr);
       plot(jetMach(ctr),stVal,mrkr,'LineWidth',1.2,'MarkerSize',6,'MarkerEdgeColor','w','MarkerFaceColor','k'); 
   end
end;  textLabels(nozzle,8);   aX.Color = [0 0 0.5];
%% Frequency spectra & Peaks
lineColor2 = GF_CustomColormap([0 0.5 1],[0 0 1;0.8 0.7 0;1 0.1 0.1],length(condition));   figure(9);   set(gcf,'Position',[175 285 825 525]);   
%  Spectra
for ctr = 1:size(splMat,3)
    imgOne(ctr) = semilogx(freq(1:2048),splMat(micNo,1:2048,ctr),'LineWidth',1.2,'Color',lineColor2(ctr,:),'DisplayName',['$NPR-',nprNames(ctr,:),'$']);   hold on;   %#ok
end;   box off;   grid on;   xlim([300 100000]);  set(gca,'FontSize',txtSiz);   aX = gca;   aX.TickLabelInterpreter = 'latex';  
title([orientName,';Mic.no:',num2str(micNo),'$'],'Interpreter','latex');   xlabel('$Frequency(Hz)$','Interpreter','latex');   ylabel('$SPL(dB/Hz)$','Interpreter','latex');
%  Peaks
for ctr = 1:size(peaksMat,2)
    plot(peaksMat(micNo,ctr,2),peaksMat(micNo,ctr,1),'h','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','c');   hold on;
end;   c3 = legend(imgOne);   c3.EdgeColor = [1 1 1];   c3.Interpreter = 'latex';   c3.Location = 'northwest';
%% Saving figures
%  Staging - SPL vs Mach No., NPR
disp([newline,'->> Staging vs Mach no.']);   figure(3); figName = ['Staging-',nozName,'-Mic-',num2str(micNo)];            commandwindow; GF_FigureSave(figName,driveOut,3);
%  SPL Peaks vs Mach No., NPR
disp([newline,'->> SPL Peaks vs Mach no.']); figure(4); figName = ['Peaks-SPL-vs-Mach-',nozName,'-Mic-',num2str(micNo)];  commandwindow; GF_FigureSave(figName,driveOut,4)
disp([newline,'->> SPL Peaks vs NPR']);      figure(5); figName = ['Peaks-SPL-vs-NPR-',nozName,'-Mic-',num2str(micNo)];   commandwindow; GF_FigureSave(figName,driveOut,5);
%  Peak Frquency(St) vs Mach No., NPR
disp([newline,'->> Freq(St) vs Mach no.']);  figure(6); figName = ['Mach-vs-St-',nozName,'-Mic-',num2str(micNo)];         commandwindow; GF_FigureSave(figName,driveOut,6);
disp([newline,'->> Freq(St) vs NPR']);       figure(7); figName = ['Mach-vs-NPR-',nozName,'-Mic-',num2str(micNo)];        commandwindow; GF_FigureSave(figName,driveOut,7);
if strcmp(Peaks,'on')
%  Staging with Peaks
   disp([newline,'->> Staging with peaks']); figure(8); figName = ['Staging-with-Peaks-',nozName,'-Mic-',num2str(micNo)]; 
else
   disp([newline,'->> Staging no peaks']);   figure(8); figName = ['Staging-no-Peaks-',nozName,'-Mic-',num2str(micNo)];
end;                                                                                                                      commandwindow; GF_FigureSave(figName,driveOut,8);
%  Frequency Spectra
disp([newline,'->>Staging Spectra']);        figure(9); figName = ['Staging-Spectra-',nozName,'-Mic-',num2str(micNo)];    commandwindow; GF_FigureSave(figName,driveOut,9);

%% ------------------------------------ Functions ---------------------------------------------------------------%%
    %% Function 1: Find the peak frequencies & OASPL
function [peaksMat,oasplMat] = PeaksOaspl(splMat,mics,freq,nozzle)
% Peak SPL & Frequency Matrix: (no. of Mics,Conditions,SPL/Freq)
peaksMat = zeros(mics,size(splMat,3),2);
% OASPL : (No. of Mics,Conditions)
oasplMat = zeros(mics,size(splMat,3));
% Instablity frequencies derived from Schlieren Fourier analysis or
% spectral peaks(Refer to Jet_Parameters for further details)
% Circular jet
instabFreqCirc   = [10900,10300,8000,7700,7350,7100,6900,6700,6550,6400,6250,6050];
% Rectangular - Minor or C-D side
instabFreqRectCD = [8100,9900,8400,8000,7650,7300,7000,6800,6550,6300,6150,5850];
% Rectangualr - Major or Flat Side
instabFreqRectFS = [8200,9900,9550,8000,7650,7300,7000,6800,6450,6300,6250,6050];
% Loop to find peak frequencies/SPL values & compute OASPL for each mic
for ctr = 1:mics 
  if strcmp(nozzle,'C')
     % Assigning peak frequencies
     peaksMat(ctr,:,2) = instabFreqCirc;
    % Finding the peak SPL levels for screech/instability inducing freq.
     for ctr2 = 1:size(peaksMat,2)
       [~,locTn] = find(freq==instabFreqCirc(ctr2));
       peaksMat(ctr,ctr2,1) = splMat(ctr,locTn,ctr2);
     end
  elseif strcmp(nozzle,'S')
      if ctr == 1 || ctr == 3
         % Assigning peak frequencies
         peaksMat(ctr,:,2) = instabFreqRectCD;
         % Finding the peak SPL levels for screech/instability inducing freq.
         for ctr2 = 1:size(peaksMat,2)
           [~,locTn] = find(freq==instabFreqRectCD(ctr2));
           peaksMat(ctr,ctr2,1) = splMat(ctr,locTn,ctr2);
         end
      else
         % Assigning peak frequencies
         peaksMat(ctr,:,2) = instabFreqRectFS;
         % Finding the peak SPL levels for screech/instability inducing freq.
         for ctr2 = 1:size(peaksMat,2)
           [~,locTn] = find(freq==instabFreqRectFS(ctr2));
           peaksMat(ctr,ctr2,1) = splMat(ctr,locTn,ctr2);
         end
      end
  end
% Computing PSD & OASPL
  PSD = reshape(splMat(ctr,:,:),size(splMat,2),size(splMat,3))';
  PSD = 10.^(PSD/10);         PSD =(PSD.*(20e-6)^2./50)./2;
  for ctr2 = 1:size(PSD,1)                                       
    oasplMat(ctr,ctr2) = trapz(freq,PSD(ctr2,:));   
    oasplMat(ctr,ctr2) = 10*log10(oasplMat(ctr,ctr2)/(20e-6)^2);
  end
end
end
    %% Function 2: Label modes on staging plot
function textLabels(nozzle,figNo)
if strcmp(nozzle,'C')
   if figNo == 8
   % Mode Type identification text - Circular
     figure(figNo);
     text(1.2,0.65,'$\bf A_{mode}$','Color','white','FontSize',15,'Interpreter','latex');
     text(1.35,0.3,'$\bf B_{mode}$','Color','white','FontSize',15,'Interpreter','latex');
   elseif figNo == 6
     text(1.215,0.65,'$\bf A_{mode}/$','Color','k','FontSize',15,'Interpreter','latex');
     text(1.32,0.3,'$\bf B_{mode}$','Color','k','FontSize',15,'Interpreter','latex');
   elseif figNo == 7
     text(2.6,0.65,'$\bf A_{mode}$','Color','k','FontSize',15,'Interpreter','latex');
     text(3.1,0.28,'$\bf B_{mode}$','Color','k','FontSize',15,'Interpreter','latex');
   end
else
   if figNo == 8
      % Mode Type identification text - Rectangular - Staging Spectra
      text(1.25,0.6,'$\bf {\hat A}_{mode}$','Color','w','FontSize',15,'Interpreter','latex');
      text(1.35,0.3, '$\bf {\hat B}_{mode}$','Color','w','FontSize',15,'Interpreter','latex');
   else
      % Mode Type identification text - Rectangular - Other plots
      text(1.25,0.65,'$\bf {\hat A}_{mode}$','Color','k','FontSize',15,'Interpreter','latex');
      text(1.35,0.3, '$\bf /{\hat B}_{mode}$','Color','k','FontSize',15,'Interpreter','latex');
   end
end
end

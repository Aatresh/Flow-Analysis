%% CODE TO CONVERT SCHLIEREN IMAGES TO FFT ACOUSTIC MAPS 
%  Accepts Image Matrices & computes FFT Energy profiles
% Configuration Legend: Active Use
% C      - CircularMedium --------------------------> Acquisition Rate: 45,000  Hz
% S      - SingleRectangular -----------------------> Acquisition Rate: 41,000  Hz
% S2/SS2 - SingleRectangular(Schlieren/Shadowgraph)-> Acquisition Rate: 204,800 Hz
% TR     - TwinRectagular --------------------------> Acquisition Rate: 41,000  Hz
% TR2    - TwinRectangular -------------------------> Acquisition Rate: 204,800 Hz
% TS     - TwinSquare ------------------------------> Acquisition Rate: 41,000  Hz
% TS1    - TwinSquare ------------------------------> Acquisition Rate: 112,000 Hz
% TS2    - TwinSquare ------------------------------> Acquisition Rate: 204,800 Hz
%---------------------------------------------------------------------------------------
tic;   fclose all;    clc;    clearvars;    set(0,'defaultfigurecolor',[1 1 1]);    code = 'SchFourier';    
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';         textSize = 13;
driveManu = [root1(1:57),'Publicaitons & Conferences\Journal Drafts\Influence of nozzle geometry' ...
    ' on screech instability closure\Physics of Fluids Draft\Submission Material\Images_Media\Media\'];
addpath(driveManu);
driveCode = [root1 'Jet_Analysis\Global_Functions\']; cd(driveCode); addpath([root1 'Jet_Analysis\Schlieren & SPOD_Codes\']);       
load custom_map3.mat;      load schl_sc.mat;     load blckRedGold.mat;   load schJet.mat;
tests = {'NPR_2p5_TR_1p0', 'NPR_2p6_TR_1p0', 'NPR_2p9_TR_1p0', 'NPR_3p0_TR_1p0',...
         'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_4p5_TR_1p0', 'NPR_5p0_TR_1p0', 'NPR_6p0_TR_1p0'};

condition = tests(1);    config = 'C';    nozzle = 'Minor';    NF = 'D';    aAmb = 345.501931; 

%  Selecting drive based on frame rate
[OutputStruct] = GF_DriveSelect(config,nozzle,code);   nozzle = OutputStruct.nozzle;   dt = OutputStruct.dt;
%  Jet Parameters
[Mj,Uj,NPR,NTR] = GF_Velocity(condition{1});           disp([newline '==> Acquisition Rate - ',num2str(1/dt),' Hz']);
%  Loading Video Matrix
driveIn  = [OutputStruct.in_root condition{1}(9:14) '\' condition{1}(1:7)  '\'];
driveOut = [OutputStruct.out_root condition{1}(9:14) '\' condition{1}(1:7) '\'];
driveCheck(driveIn);     load([driveIn condition{1}(1:7) '_DAT']);    load([driveIn 'X']);   load([driveIn 'Y']);
%% 1. * LOADING BACKGROUND FILE
if isempty(OutputStruct.bckgrnd_root) ~= 1 
   disp([newline '==>> Subtracting Background' newline]);
   Bckgrnd = load([OutputStruct.bckgrnd_root 'Bckgrnd_DAT']);    cLim1 = [0 1.8];
   Bckgrnd = Bckgrnd.Master_U;   M2 = Master_U./mean(Bckgrnd,3);    
   avgOne = mean(M2,3);          avgTwo = mean(Master_U,3);      M2Fluc = M2 - avgOne;
else
   Bckgrnd = [];  cLim1 = [0 500];    avgOne = mean(Master_U,3); avgTwo = avgOne; 
end; disp([newline,'    <====   Done!   ====>']);  load schJet.mat;   bckMean = mean(Bckgrnd,3);   clear Bckgrnd;   
%  Axis Normalization & Labelling
[Xn,Yn,limX,limY,xName,yName,lenScales,figSize] = GF_AxisDefnSch(config,nozzle,NF,X,Y);   Xn = Xn';   Yn = Yn';
maxX = round(max(Xn),2);   minX = round(min(Xn),2);   maxY = round(max(Yn),2);   minY = round(min(Yn),2);   
if strcmp(config,'TS') && strcmp(nozzle,'Minor')
   [avgOne,M2Fluc] = bckgCorctnOne(Master_U,M2,Bckgrnd);
end
%  Original fluctuating matrix; dimension variables; data clearing
Master_U2 = Master_U - avgTwo;   rowNos = size(Master_U2,1);     colNos = size(Master_U2,2);    imgNos = size(Master_U2,3);  
clear M2 Master_U;
%  Function to generate tag, lipline location & nozzle config for filename
[tag,lipLine,fNameCon] = tagNLine(NF,config,nozzle,lenScales);
%% 2. * PLOTTING AVERAGE SCHLIEREN
fig1 = figure;   pcolor(Xn,Yn,avgOne);    shading interp,    colormap gray,    axis equal;    caxis(cLim1);  set(gca,'FontSize',13);  
xlim(limX);      xlabel(xName,'Interpreter','latex');        aX = gca;         aX.TickLabelInterpreter = 'latex';
ylim(limY);      ylabel(yName,'Interpreter','latex');        set(gcf,'Position',figSize);  
title('$Average \thinspace Schlieren \thinspace Image$','Interpreter','latex');   hold on;       
%  Plotting Shock Cell Locations
shckLocPlot(driveOut,NF,lenScales,Xn,Yn,fig1,'k');
%%    => 2.1 FIGURE SAVE: AVERAGE SCHLIEREN 
figName = ['Avg_Schlieren',tag,'-',condition{1}];             GF_FigureSave(figName,driveOut,fig1.Number);
%% 3.    PLOT CENTERLINE 
[Shcks,Expnsn,shockCellWidth] = SchFn_ShockCellLoc(Xn,Yn,avgOne,config,nozzle,xName,yName,NPR,limY);
avgShockWidth = mean(shockCellWidth(2:end),2);    disp([newline,'->> Average shock cell width = ',num2str(avgShockWidth)]);
% disp([newline,'-> Shock & epxansion locations computed with shorter X axis(starting from 0)'])
% save([drive_out 'ShckExpansn'],'Shcks','Expnsn');
% figName = ['Shock_Expansion_Contour',tag,'-',condition{1}];   GF_Figure_Save(figName,drive_out,figX.Number);
% figName = ['Shock_Expansion_Locations',tag,'-',condition{1}]; GF_Figure_Save(figName,drive_out,figY.Number);
%  SAVE AVG SCHLIEREN MAT FILE
% save([driveOut 'SC - Average' tag],'avgOne','Xn','Yn');
%% 4. * SPATIAL FFT: AVERAGE IMAGE CENTERLINE
inpStruct.Xn = Xn;     inpStruct.avgOne = avgTwo;   inpStruct.lenScales = lenScales;   
inpStruct.Yn = Yn;     inpStruct.NPR    = NPR;      inpStruct.leakAge   = 1;
[pltFig,pSpec,kPspecH,kPosPeak] = SchFn_CenLineSPFFT(inpStruct);
%% 5.    STANDARD DEVIATION
stdMap = std(Master_U2,0,3);    stdMap = log10(stdMap);      
fig2 = figure;    contourf(Xn,Yn,stdMap,30);    colormap(schJet),    axis equal,  ylim(limY),    xlim(limX);    
c = colorbar;     set(gcf,'Position',figSize);  caxis([1.3 2.68]);   aX = gca;    aX.TickLabelInterpreter = 'latex';    
xlabel(xName,'Interpreter','latex');    ylabel(yName,'Interpreter','latex');      set(gca,'FontSize',13); 
title(c,'$log(\sigma_{\partial \rho})$','Interpreter','latex');                   c.TickLabelInterpreter = 'latex';
title('$Standard \thinspace Deviation$','Interpreter','latex');                                                                   
% Plotting Shock Cell Locations
shckLocPlot(driveOut,NF,lenScales,Xn,Yn,fig2,'w')
%%    => 4.1 FIGURE SAVE: INTENSITY STANDARD DEVIATION
figName = ['Log(Std_Dev)',tag,'-',condition{1}];    GF_FigureSave(figName,driveOut,fig2.Number);    
%% 6.   PLOTTING THE STD DEV DISTRIBUTION TO COMPUTE STANDING WAVE WAVELENGTH
strtPos = 0.7;    endPos = 1;       [radPos,~]  = find((strtPos<Yn) & (Yn<endPos));     stdVals = stdMap(radPos,:);   
meanStdVals = mean(stdVals);        maxVals = max(stdVals);     figure;     
plot(Xn,stdVals,'.c');    hold on;       plot(Xn,meanStdVals,'ok','LineWidth',1.2);    grid on;    ax = gca;    ax.FontSize = 13;           
xlabel(xName,'Interpreter','latex');     ylabel('$\sigma_I$','Interpreter','latex');   xlim([0.5 5.5]);  box off;
plot(Xn,maxVals,'sr','LineWidth',1.2);   ax.TickLabelInterpreter = 'latex';     set(gcf,'Position',[30 250 1250 570]);   
title(['$\bf Standard \thinspace Deviation(yRange = ',num2str(strtPos),NF,':',num2str(endPos),NF,')$'],'Interpreter','latex');    clear meanVals maxVals   
%% 7. * RESOLVED INTENSITY & INTEGRAL ENERGY CALCULATION
[blockSize,resIntnstySclr] = blkSizer(config);                          
% Method-1: Using background image as basline value & mulitplying by factor to enhance variations
% pRefBase = mean(mean(bckMean))*6E-10;                 % without background image subtraction
% New Method: Reference pressure value is now computed based on the acoustic OASPL
% distribution obtained from Near field Data
pRefSub  = mean(mean(bckMean))*6*10^(-resIntnstySclr);   % with background image subtraction
nOvlp    = 0.75;                                         % percentage of overlap
nBlks    = (size(Master_U2,3)- (blockSize*nOvlp))/(blockSize - (blockSize*nOvlp)); 
schFreq  = 0:1/dt/blockSize:1/dt-1;
%  Without background subtraction for intgEnergy
% [~,~,~,~,intgEng] = Sch_ImgFFT(Master_U2,blockSize,nOvlp,dt,pRefBase);
%  With background subtraction for frequency SPL
[imgPSD,imgFFT,imgPhs,imgSpec,intgEng] = SchFn_ImgFFT(M2Fluc,blockSize,nOvlp,dt,pRefSub);
disp([newline '-------------------------------------------']);
disp(['=> Number of blocks     :- ',num2str(floor(nBlks))]);
disp(['=> Frequency resolution :- ',num2str(schFreq(2))]);
disp('-------------------------------------------');
%% 8.   PLOTTING RESOLVED INTENSITY SPECTRUM AT (X,Y)
prompt = [newline '-> Enter Axial  point for SPL(between:',num2str(minX),' to ',num2str(maxX),'): '];  inpX = input(prompt);
prompt = [newline '-> Enter Radial point for SPL(between:',num2str(minY),' to ',num2str(maxY),'): '];  inpY = input(prompt);
[xCol,xVal,yRow,yVal] = precisLoc(inpX,inpY,Xn,Yn);         ptSpec = reshape(imgSpec(yRow,xCol,:),[1,size(imgSpec,3)]);  
fig10 = figure;    semilogx(schFreq,ptSpec,'LineWidth',1);   grid on;   hold on;   set(gcf,'Position',[60 110 620 405]);
xlabel('$\bf Hz$','Interpreter','latex');        xlim([300 (1/dt)/2]);   box off;             
ylabel('$\bf SPL(dB)$','Interpreter','latex');   ylim([95 135]);         set(gca,'FontSize',13,'TickLabelInterpreter','latex');              
title(['$SPL \thinspace at',tag,':x=',num2str(xVal),',y=',num2str(yVal),'$'],'Interpreter','latex');  
%%    => 8.1 FIGURE SAVE: SPATIAL RESOLVED INTENSITY SPECTRA
figName = ['SPL-Spectra-X_',num2str(xVal),'-Y_',num2str(yVal),tag,'-',condition{1}];       GF_FigureSave(figName,driveOut,fig10.Number); clear xLoc yLoc;
%% 9. * PLOT FREQUENCY SPECIFIC RESOLVED INTENSITY & PHASE
prompt = [newline '-> Enter Frequency to plot(range:300-',num2str((1/dt)/2),'): '];        St = schFreq*lenScales(1)/(Uj);
inpFreq = input(prompt); % [~,p] = find(freq>inp_freq); pos = p(1)-1;
freqPos = round(inpFreq/schFreq(2))+1;    cLim2 = [100 150];    cLim3 = [135 170];  figSizPhs = [790 380 620 405];
%  Plotting SPL
fig3 = figure;    pcolor(Xn,Yn,imgSpec(:,:,freqPos)),    shading interp, axis equal,    colormap(schJet);
caxis(cLim2);     xlim(limX),        ylim(limY);         xlabel(xName,'Interpreter','latex');          ylabel(yName,'Interpreter','latex');
title(['$Resolved \thinspace Intensity(\zeta), St \approx ',num2str(round(St(freqPos),3)),'$'],'Interpreter','latex');
c = colorbar;     set(gcf,'Position',figSize);           title(c,'$\zeta_{i}$','interpreter','latex');             
set(gca,'FontSize',13);              aX = gca;           aX.TickLabelInterpreter = 'latex';            c.TickLabelInterpreter = 'latex';
%  Plotting Shock Cell Locations
shckLocPlot(driveOut,NF,lenScales,Xn,Yn,fig3,'k');   box off;
%  Plotting Phase 
fig4 = figure;    pcolor(Xn,Yn,-angle(imgPhs(:,:,freqPos))),     shading interp, axis equal,    caxis([-pi pi]); colormap bluewhitered; 
xlim(limX),       ylim(limY);        xlabel(xName,'Interpreter','latex');  ylabel(yName,'Interpreter','latex'); 
c = colorbar;     c.Ticks = [-pi,-pi/2,0,pi/2,pi];       set(gcf,'Position',figSizPhs); c.TickLabelInterpreter = 'latex';
title(['$St \approx',num2str(round(St(freqPos),2)),'$'],'Interpreter','latex');
title(c,'$\Delta\phi$','Interpreter','latex');           c.TickLabels = {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'};  
set(gca,'FontSize',13);              aX = gca;           aX.TickLabelInterpreter = 'latex';        %aX.XTick = [0,2,4,6,8];
%  Plotting Shock Cell Locations
shckLocPlot(driveOut,NF,lenScales,Xn,Yn,fig4,'w');   box off;
%%    => 9.1 FIGURE SAVE: RESOLVED INTENSITY & PHASE
figName = ['SPL-Energy-Dist',tag,'-',condition{1},'-Freq-',num2str(schFreq(freqPos))];  commandwindow;    GF_FigureSave(figName,driveOut,fig3.Number);
figName = ['Phase-Dist',tag,'-',condition{1},'-Freq-',num2str(schFreq(freqPos))];       commandwindow;    GF_FigureSave(figName,driveOut,fig4.Number); 

%% 10.  PHASE PROPAGATION .gif
nt = 80;    T = 1/schFreq(freqPos);   time = linspace(0,T,nt);       commandwindow;
%  Adding above parameters to input Struct file 
phaseLoopInput.time = time;    phaseLoopInput.nt = nt;           phaseLoopInput.TYP = code;        phaseLoopInput.config = config;
phaseLoopInput.freq = schFreq; phaseLoopInput.freqNo = freqPos;  phaseLoopInput.St = St;           phaseLoopInput.Sz = figSizPhs;
phaseLoopInput.Xn = Xn;        phaseLoopInput.nozzle = nozzle;   phaseLoopInput.Yn = Yn;           phaseLoopInput.NF = NF;  
phaseLoopInput.lim_x = limX;   phaseLoopInput.lim_y = limY;      phaseLoopInput.clim = [-pi pi];   phaseLoopInput.cMap = schJet;
phaseLoopInput.gifName = ['PhaseShift-',num2str(schFreq(freqPos)),'Hz']; phaseLoopInput.driveOut = driveOut;
GF_GIF_DisplayWrite(imgPhs,phaseLoopInput); 
%% 11.  PLOTTING PHASE TRANSITION
%  Time constant & Phase position counter
nt = 29;    T = 1/schFreq(freqPos);   time = linspace(0,T,nt);   time2 = time/T;   lipLineLoc = zeros(size(Xn));   lipLineLoc(:) = 0.5;
for ctr = 1:nt
%  Plotting phase
    cFig = figure;   pcolor(Xn,Yn,-angle(squeeze(imgPhs(:,:,freqPos)*exp(2i*pi*schFreq(freqPos)*time(ctr)))));   shading interp;   axis equal;   box off;
    xlabel(xName,'Interpreter','latex');   xlim([0 6]),    c = colorbar;   colormap bluewhitered;   caxis([-pi pi]);   set(gcf,'Position',[55 480 700 240]);            
    ylabel(yName,'Interpreter','latex');   ylim([0 2]);    c.Ticks = [-pi,-pi/2,0,pi/2,pi];   c.TickLabelInterpreter = 'latex';   aX = gca;
    title(['$\tau/T =',num2str(time2(ctr)),'$'],'Interpreter','latex');    aX.FontSize = textSize;   aX.TickLabelInterpreter = 'latex';   %aX.XTick([0 2 4 6]);        
    title(c,'$\Delta\phi$','Interpreter','latex');    c.TickLabels = {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'};  
    disp(['Step = ',num2str(ctr-1)]);   hold on;   set(gca,'Layer','top');   plot(Xn,lipLineLoc,':w','LineWidth',1.2);   aX.TickLength = [0 0]; 

%  Defautl figure save command
%     saveDrive = 'X:\OneDrive - University of Cincinnati\Working_Directory\Publicaitons & Conferences\Journal Drafts\Effect of Flow Instabilities on Circular & Rectangular Supersonic Jets\Manuscript Images\PhaseShift-Circular-NPR2p5-10900Hz\V2\';
%     figName = ['Interval-',num2str(time2(ctr)),'T'];   commandwindow;   GF_FigureSave(figName,saveDrive,cFig.Number);  
%  Command set to save high resolution .tif files
%     cd(saveDrive);   figName = ['Interval-',num2str(time2(ctr)),'T.tif'];   print(cFig, '-dtiffn', '-r300', figName);   cd(driveCode);
%  Shock cell locations
shckLocPlot(driveOut,NF,lenScales,Xn,Yn,cFig,'w')
end
%% 12.  PLOT PHASE & INTENSITY VARIATION AT 3 RADIAL DISTANCES(0.5,0.55,0.6)
%  Execute Sec. 10 first
[d1,d2,d3,d4,d5] = radDistAlloc(config,NF);
rD1 = find(Yn>d1);     rD2 = find(Yn>d2);     rD3 = find(Yn>d3);   rD4 = find(Yn>d4);   rD5 = find(Yn>d5);       
%  Recompiling the phase to get the best distribution from module above
ctr = 7;         newImgPhs = -angle(imgFFT(:,:,freqPos)*exp(2i*pi*schFreq(freqPos)*time(ctr)));
phsRD1 = newImgPhs(rD1(1),:);         phsRD2 = newImgPhs(rD2(1),:);         phsD3 = newImgPhs(rD3(1),:);          phsD4 = newImgPhs(rD4(1),:);          phsD5 = newImgPhs(rD5(1),:);
intRD1 = imgSpec(rD1(1),:,freqPos);   intRD2 = imgSpec(rD2(1),:,freqPos);   intRD3 = imgSpec(rD3(1),:,freqPos);   intRD4 = imgSpec(rD4(1),:,freqPos);   intRD5 = imgSpec(rD5(1),:,freqPos);  
%  Plotting phase at three specific radial locations
fig5 = figure;      plot(Xn,phsRD1,'r','LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d1),'$']);   grid on;   hold on;   box off;  set(gcf,'Position',[30 250 850 400]);
plot(Xn,phsRD2,'-.b','LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d2),'$']);   xlabel(xName,'Interpreter','latex');   ax1 = gca; 
plot(Xn,phsD3 ,'--k','LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d3),'$']);   ylabel('$Phase(\Delta \phi)$','Interpreter','latex');
plot(Xn,phsD4 ,':r', 'LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d4),'$']);   plot(Xn,phsD3 ,':b','LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d5),'$']);
xlim([0.5 5.5]);    ylim([-pi pi]);    ax1.YTick = [-pi,-pi/2,0,pi/2,pi];   ax1.YTickLabel = {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'};
ax1.FontSize = 13;   ax1.TickLabelInterpreter = 'latex';   lG = legend;   lG.Interpreter = 'latex';   lG.EdgeColor = [1 1 1];  lG.Location = 'northwest';
%  Plotting intensity at three specific radial locations
fig6 = figure;      plot(Xn,intRD1,'r','LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d1),'$']);   grid on;   hold on;   box off;  set(gcf,'Position',[30 250 850 400]);   %#ok
plot(Xn,intRD2,'-.b', 'LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d2),'$']);   xlabel(xName,'Interpreter','latex');   ax1 = gca; 
plot(Xn,intRD3 ,'--k','LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d3),'$']);   ylabel('$Intensity(\zeta_{i})$','Interpreter','latex');
plot(Xn,intRD4 ,':r', 'LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d4),'$']);   plot(Xn,intRD5 ,':b','LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d5),'$']);
xlim([0.5 5.5]);   ylim([100 150]);   ax1.FontSize = 13;   ax1.TickLabelInterpreter = 'latex';   lG = legend;   lG.Interpreter = 'latex';   lG.EdgeColor = [1 1 1];  lG.Location = 'southeast';

%% 13.  FIND PEAKS IN RESOLVED INTENSITY PLOTS AT SPECIFIC FREQUENCY
% Compute the average of resolved intensity distribution at various radial locations & fit a 
% smoothing spline to find the peak locations
strtPos = 0.6;    endPos = 1;     [radPos,~] = find((strtPos<Yn) & (Yn<endPos));   
splVals = imgSpec(radPos,:,freqPos);    meanVals = mean(splVals);    maxVals = max(splVals);   figure;    plot(Xn,splVals,'.c');    hold on;   ax = gca;  
plot(Xn,meanVals,'ok','LineWidth',1.2);   grid on;  hold on;     set(gcf,'Position',[30 250 1250 570]);    ax.FontSize = 13;     ax.TickLabelInterpreter = 'latex';
xlabel(xName,'Interpreter','latex');    ylabel('$Resolved \thinspace Intensity(\zeta)$','Interpreter','latex');  xlim([0.5 5.5]);    box off;
plot(Xn,maxVals,'sr','LineWidth',1.2);
%  Plotting Enerygy at three specific radial locations
fig6 = figure;      plot(Xn,imgSpec(rD1(1),:,freqPos),'r','LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d1),'$']);   grid on;   hold on;   box off;
plot(Xn,imgSpec(rD2(1),:,freqPos),'-.b','LineWidth',1.2,'DisplayName' ,['$\bf ',yName(2:end-1),' = ',num2str(d2),'$']);   xlabel(xName,'Interpreter','latex');   ax2 = gca;
plot(Xn,imgSpec(rD3(1),:,freqPos), '--k','LineWidth',1.2,'DisplayName',['$\bf ',yName(2:end-1),' = ',num2str(d3),'$']);   ylabel('$Resolved \thinspace Intensity(\zeta)$','Interpreter','latex');
plot(Xn,imgSpec(rD4(1),:,freqPos), '.m','LineWidth',1.2,'DisplayName' ,['$\bf ',yName(2:end-1),' = ',num2str(d4),'$']);   xlim([0.5 4.5]);    ylim(cLim2);    ax2.FontSize = 13;   
ax2.TickLabelInterpreter = 'latex';   lG = legend;  lG.Interpreter = 'latex';   lG.EdgeColor = [1 1 1];      lG.Location = 'southeast';      set(gcf,'Position',[30 250 1250 570]);
%  Computing & plotting mean
meanInts = [];    meanInts(1,:) = imgSpec(rD1(1),:,freqPos);    meanInts(2,:) = imgSpec(rD2(1),:,freqPos);    meanInts(3,:) = imgSpec(rD3(1),:,freqPos);    meanInts(4,:) = imgSpec(rD4(1),:,freqPos);
meanInts = mean(meanInts);   fig6.Number;   plot(Xn,imgSpec(rD2(1),:,freqPos),'-ok','LineWidth',1.2,'MarkerFaceColor','g','MarkerEdgeColor','k','MarkerSize',4,'DisplayName','$\bf Mean \thinspace Value$');
%  Fitting a smoothing spline to the centerline data to get peak locations
[xPeak, intsPeak] = prepareCurveData(Xn,meanInts);       smoothingParam = 0.9997;
%  Assigning Fit Type & Fitting options
ft = fittype( 'smoothingspline' );    opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';                opts.SmoothingParam = smoothingParam;
%  Fitting to the intensity & extracting peak locations
[fitIntensty, ~] = fit(xPeak,intsPeak,ft,opts );     cofVal = coeffvalues(fitIntensty);    intsVals = cofVal.coefs(:,4);    clear cofVal;   [~,p1] = findpeaks(intsVals);    
plot(Xn(1:end-1),intsVals,'.c','DisplayName','$\bf Spline \thinspace Fit$');   plot(Xn(p1),intsVals(p1),'d','MarkerFaceColor','c','MarkerEdgeColor','k','LineWidth',1.2,'HandleVisibility','off');
%% 14.  PHASE PROPAGATION .mp4
prompt = [newline,'->> Generate .mp4?(y/n) - '];                     chc = input(prompt,'s');
movieInput.fTitle = ['$\bf Movie - 7,','Flat \thinspace Side; NPR \thinspace',num2str(NPR),';\it f_s:',num2str(schFreq(freqPos-1)),' Hz$']; 
movieInput.videoName = 'Fig-15(d)-C-D Side-Movie';
% imgPhs = imgPhs*exp(2i*pi*schFreq(freqPos)*time(4));    % use for NPR 2.6 Major 9950Hz; nt = 29
figSiz = [80 220 640 500];
if strcmp(chc,'y')
   nt = 80;    T = 1/schFreq(freqPos);   time = linspace(0,T,nt);     commandwindow;
   movieInput.Xn     = Xn;        movieInput.Yn        = Yn;          movieInput.St       = St;
   movieInput.code   = code;      movieInput.time      = time;        movieInput.driveOut = driveManu;
   movieInput.cLim   = [-pi pi];  movieInput.StNo      = freqPos;     movieInput.freq     = schFreq;
   movieInput.figSiz = figSiz;    movieInput.schJet    = schJet;      movieInput.config   = config;
   movieInput.nt     = nt;        movieInput.xName     = xName;       movieInput.yName    = yName; 
   movieInput.nozzle = nozzle;    movieInput.driveCode = driveCode;   movieInput.style    = 'hf1';                                      
   movieInput.descrpText1 = 'Karnam, Mohammad, Gutmark';
   movieInput.descrpText2 = '$Physics \thinspace of \thinspace Fluids$';
%  Expected source Location
%    sX = [1.65,1.656,2.309,2.309];   sY = [0.425,-0.425,0.425,-0.425];   sourceLocX = movieInput.sX;   sourceLocY = movieInput.sY;
%    save([driveOut,['sourceLocX-',num2str(schFreq(freqPos))]],'sourceLocX');   save([driveOut,['sourceLocY-',num2str(schFreq(freqPos))]],'sourceLocY');
%    movieInput.videoName = ['PhaseShift',tag,'-',nozzle,'-NPR-',num2str(NPR),'-Freq-',num2str(schFreq(freqPos)),'Hz.mp4'];
%    movieInput.fTitle = ['$St \approx',num2str(round(St(freqPos),2)),'$'];
   GF_MovieMaker(movieInput,imgPhs);     disp([newline,'->> Video file saved']);  
end
%% 15.  EXPORTING TO TECPLOT FORMAT
prompt = [newline,'->> Choose parameter to export: Phase(1); Energy distribution(2); Cancel(3)',newline,'->> Input: '];     chc = input(prompt);
if chc < 3
   writeInp.X = Xn;                           writeInp.Y = Yn;                      writeInp.chc = chc;
   writeInp.codeRoot = root1;                 writeInp.driveOut  = driveOut;
   writeInp.nozzle   = nozzle;                writeInp.condition = condition{1};
   if chc == 1
      writeInp.contour  = imgPhs(:,:,freqPos);    writeInp.figName  = ['Phase-Dist',tag,'-',condition{1},'-Freq-',num2str(schFreq(freqPos)),'.plt'];
   else
      writeInp.contour  = imgSpec(:,:,freqPos);   writeInp.figName  = ['SPL-Energy-Dist',tag,'-',condition{1},'-Freq-',num2str(schFreq(freqPos)),'.plt']; 
   end
   tecWrite(writeInp);
end
%% 16.  INTEGRAL ENERGY
fig7 = figure;    pcolor(Xn,Yn,intgEng),  shading interp,    axis equal,    colormap(schJet);    caxis(cLim3);
xlim(limX),       xlabel(xName,'Interpreter','latex');       ylim(limY);    ylabel(yName,'Interpreter','latex');
c = colorbar;     title(c,'\Sigma\zeta','FontWeight','bold');               set(gcf,'Position',figSize); 
ax = gca;         ax.TickLabelInterpreter = 'latex';         set(gca,'FontSize',13);             c.TickLabelInterpreter = 'latex'; 
title('$Integral \thinspace Energy$','Interpreter','latex');             
%%    => 16.1 FIGURE SAVE: INTEGRAL ENERGY
figName = ['OASPL-Dist',tag,'-',condition{1}];      GF_FigureSave(figName,driveOut,fig7.Number);
%% 17.  PLOT INTEGRAL OF FREQUENCY RANGE
%  Starting & Ending Frequencies
prompt = [newline '-> Enter Starting Frequency (range:300-',num2str((1/dt)/2),'): '];  strtFreq = input(prompt);
prompt = [newline '-> Enter Ending Frequency (range:300-',num2str((1/dt)/2),'): '];    endFreq  = input(prompt);
%  Frequency resolution                     Rounding to remove any smallfractions
delFreq = schFreq(3) - schFreq(2);         intgrlFreq = round(schFreq); 
%  Location of frequencies in Frequency Matrix & selecting those range of frequencies 
[~,strtLoc] = find(intgrlFreq == strtFreq);    [~,endLoc] = find(intgrlFreq == endFreq);
freqRange   = schFreq(strtLoc:endLoc);
%  Range of energy levles to be integrated
energyRange = imgSpec(:,:,strtLoc:endLoc);
%  Reversing log operation to get back real part of FFT 
fftRealVals = 10.^(energyRange/10);         fftRealVals   = (fftRealVals.*(20e-6)^2./50)/2;
fftRealVals = reshape(fftRealVals,[rowNos*colNos,size(fftRealVals,3)]);
%  Defining intergal matrix 
energyIntgrl = zeros(size(fftRealVals,1),1);
%  Performing numerical integration
wB = waitbar(0,'Computing Energy Integral'); 
for ctr = 1:size(fftRealVals,1)
    energyIntgrl(ctr) = trapz(freqRange,fftRealVals(ctr,:)); 
    energyIntgrl(ctr) = 10*log10(energyIntgrl(ctr)/(20e-6)^2);  waitbar(ctr/size(fftRealVals,1),wB);
end;   close(wB);
%  Reshaping back to image size
energyIntgrl = reshape(energyIntgrl,rowNos,colNos);
%  Plotting Selective Integral Energy 
fig8 = figure;     pcolor(Xn,Yn,energyIntgrl),  shading interp,      axis equal,      colormap(schJet);
caxis([150 180]);  xlim(limX), ylim(limY);      xlabel(xName,'Interpreter','latex');  ylabel(yName,'Interpreter','latex');
title(['$\Sigma\epsilon, St_{range} \approx ',num2str(round(St(strtLoc),3)),':',num2str(round(St(endLoc),3)),'$'],'Interpreter','latex');
c = colorbar;      title(c,'$\Sigma\epsilon_{range}$','interpreter','latex');         set(gcf,'Position',figSize); 
set(gca,'FontSize',13);    ax = gca;            ax.TickLabelInterpreter = 'latex';    c.TickLabelInterpreter = 'latex';   
%%    => 17.1 FIGURE SAVE: SELECTIVE INTEGRAL ENERGY
figName = ['Enrgy-Intgrl',tag,'-FreqRng-',num2str(schFreq(strtLoc)),'-',num2str(schFreq(endLoc))];   GF_FigureSave(figName,driveOut,fig8.Number);
%% 18.  ABSOLUTE OF MAGNITUDE OF SELECTED FREQUENCY
cohrFlucMag = abs(imgFFT(:,:,freqPos)*5);   GF_FigurePlot(cohrFlucMag,Xn,Yn);      fig9 = gcf;   xlim(limX);   ylim(limY);   colormap(schJet);     ax = gca;
xlabel(xName,'Interpreter','latex');      ylabel(yName,'Interpreter','latex');   ax.FontSize = textSize;     ax.TickLabelInterpreter = 'latex';  
title('$\bf Absolute \thinspace Fluctuation \thinspace Magnitude$','Interpreter','latex');   c = colorbar;   c.TickLabelInterpreter = 'latex';   caxis([0 50])
%%    => 18.1 FIGURE SAVE: ABSOLUTE FLUCTUATION MAGNITUDE
figName = ['Fluc-Magnitdue',tag,'-Freq-',num2str(schFreq(freqPos))];   GF_FigureSave(figName,driveOut,fig9.Number);
%% 19.  TWO POINT SPATIAL CORRELATION 
prompt = [newline '-> Enter Reference Axial location(between:0 to ',num2str(maxX),'): ']; inpX = input(prompt);
prompt = ['-> Enter Reference Radial location(between:',num2str(minY),' to ',num2str(maxY),'): ']; inpY = input(prompt);
%  Getting the precise location from input co-ordinates            Covariance Matrix
[colRef,xVal,rowRef,yVal] = precisLoc(inpX,inpY,Xn,Yn);           R_cov = zeros(size(avgTwo));
%  Dashed line to indicate Reference Location
ref_X = zeros(size(Yn));   ref_Y = zeros(size(Xn));   colr1 = [0.3 0.3 0.3];
ref_X(:) = xVal;           ref_Y(:) = yVal;    
%  Standard Deviation at Reference Location
stdRef = std(Master_U2(rowRef,colRef,:));
%  Computing Two-Point Spatial Correlation at all locations
wB = waitbar(0,'Computing Two-Point Spatial Co-relation');    
for row = 1:rowNos
  for col = 1:colNos
    R_cov(row,col) = (sum(Master_U2(row,col,:).*Master_U2(rowRef,colRef,:),3))/(stdRef * std(Master_U2(row,col,:)));
    waitbar(row/rowNos,wB);
  end
end; close(wB);
%  Normalizing with peak value & display Reference Location
R_cov = R_cov/max(max(R_cov));   disp([newline,'>> - X Ref: ',num2str(Xn(colRef))]);   disp([newline,'>> - Y Ref: ',num2str(Yn(rowRef))])  
%  Legend 
fig10 = figure;   pcolor(Xn,Yn,R_cov),   axis equal,   shading interp;   caxis([-1 1]);   hold on;   colormap bluewhitered;
plot(ref_X,Yn,'--','Color',colr1,'LineWidth',0.5);     plot(Xn,ref_Y,'--','Color',colr1,'LineWidth',0.5);
xlabel(xName,'Interpreter','latex');   xlim(limX);     legCond = ['NPR - ',num2str(NPR)];            box off;  
ylabel(yName,'Interpreter','latex');   ylim(limY);     legXloc = ['$X_{ref}: ' num2str(xVal),'$'];   legYloc = ['$Y_{ref}: ' num2str(yVal),'$'];
title(['$Reference \thinspace Loc. - Row:',num2str(rowRef),';Col.:',num2str(colRef),'$'],'Interpreter','latex'); 
c = legend(legCond,legXloc,legYloc);   set(gcf,'Position',[105 310 735 400]);   title('$Two \thinspace Point \thinspace Spatial \thinspace Correlation (R_{\Delta_x})$','Interpreter','latex');
c.EdgeColor = [1 1 1];   c.Location = 'northeastoutside';   c.Interpreter = 'latex';   ax = gca;   ax.TickLabelInterpreter = 'latex';   ax.FontSize = 12;
%%    => 19.1  FIGURE SAVE: TWO-POINT SPATIAL CORRELATION
figName = ['TwoPointSCorrel-X-',num2str(xVal),'Y-',num2str(yVal),tag];      GF_FigureSave(figName,driveOut,fig10.Number);
%% 20.* COMPUTING SHOCK-COHERENT INTERACTION WAVENUMBERS
%  Picking point on plot for wave number of first shock - Ks1
disp([newline,'=> Pick first three shock wave numbers']);   nozHt = lenScales(2);   figure(pltFig.Number);   [x,~] = ginput(3);   [~,zLoc] = find(kPspecH==0);                  
loc1 = find(kPspecH(kPosPeak+zLoc)<x(1));     loc2 = find(kPspecH(kPosPeak+zLoc)<x(2));     loc3 = find(kPspecH(kPosPeak+zLoc)<x(3));
kS1  = kPspecH(kPosPeak(loc1(end)+1)+zLoc);   kS2  = kPspecH(kPosPeak(loc2(end)+1)+zLoc);   kS3  = kPspecH(kPosPeak(loc3(end)+1)+zLoc);
%  Wavenumbers: kKH - KH waves; kIn1 - KH Shock-1 Interaction; kIn2 - KH Shock-2 Interaction; kIn3 - KH Shock-3 Interaction   
uConv = 0.71;   kKh = nozHt*(2*pi*schFreq(freqPos))/(uConv*Uj);   
kIn1 = kKh - kS1;   kIn2 = kKh - kS2;   kIn3 = kKh - kS3;  
disp([newline,'==> First Shock Wave Number   : ',num2str(kS1)]);
disp([newline,'==> Second Shock Wave Number  : ',num2str(kS2)]);
disp([newline,'==> Third  Shock Wave Number  : ',num2str(kS3)]);
disp([newline,'==> KH Instability Wave Number: ',num2str(kKh)]);
disp([newline,'==> KH Shock-1 Interaction    : ',num2str(kIn1)]);
disp([newline,'==> KH Shock-2 Interaction    : ',num2str(kIn2)]);
disp([newline,'==> KH Shock-3 Interaction    : ',num2str(kIn3)]);
%  Re-scaling X-axis to start from nozzle exit & using equivalent diameter to normalize length scales 
[~,xStrt] = find(Xn==0);    xNew = Xn(xStrt:end);   nozHt = lenScales(2);
%% 21.  SPATIO-TEMPORAL FOURIER TRANSFORM
%  Run Sec. 19 - SSFT on avg Schlieren first
%  Taking radial location as input 
prompt = [newline,'=>> Enter start radial location to extract data(',num2str(round(min(Yn),2)),' to ',num2str(round(max(Yn),2)),')- '];  
radLoc = input(prompt);   radRange = radLoc:0.2:radLoc+2;   fig11 = figure;   mainPSD = [];   fig12 = figure;
%  Average Schlieren Plot
figure(fig11.Number);  pcolor(Xn,Yn,avgOne);    shading interp,    colormap gray,    axis equal;    caxis(cLim1);  set(gca,'FontSize',13);  
xlim([0 6.5]);   xlabel(xName,'Interpreter','latex');        aX = gca;         aX.TickLabelInterpreter = 'latex';
ylim([-2 2]);    ylabel(yName,'Interpreter','latex');        set(gcf,'Position',figSize);   set(gca,'Layer','top');     
hold on;         box off;
for ctr = 1:size(radRange,2)
    if radRange(ctr) < 0
       [locTn,~] = find(Yn > radRange(ctr));    locTn = locTn(1) - 1;     
    else
       [locTn,~] = find(Yn < radRange(ctr));    locTn = locTn(end) + 1;
    end
%  Plotting data line on average figure
    yLoctn = zeros(size(Xn));   yLoctn(:) = Yn(locTn);   figure(fig11.Number); 
    plot(Xn,yLoctn,':w','LineWidth',1.2) ;   ax = gca;   ax.TickLength = [0 0];   clear ax; 
%  Extracting values at specified radial location for fftn input
    inpSignal = Master_U2(locTn,xStrt:end,:);    inpSignal = reshape(inpSignal,size(xNew,2),size(Master_U2,3),[]);
%  Checking if number of values is odd for centering
    if rem(size(inpSignal,1),2) == 0
       inpSignal = inpSignal(1:end-1,:);
    end
% FFT Block size, Frequency resolution & Frequency Range(using global function)
    [blockSize] = blkSizer(config);      delF = (1/dt)/blockSize;         freqRange = GF_FourierFreq(blockSize,1/dt); %freqRange = 0:delF:(1/dt)-delF;
% Nyquist Cut-off & Positive Frequency values
    nCut = blockSize/2;                 postFreq = freqRange(1:blockSize/2);       
%  Strouhal number normalized by height of the nozzle, convective velocity factor & axial resolution in mm
    stH = (postFreq*nozHt)/Uj;          dx = (xNew(3) - xNew(2))*lenScales(1);
%  Wavenumber definition
    kAct = GF_FourierFreq(size(xNew,2),2*pi/dx);   kNos = fftshift(kAct);    kRangeH = kNos*nozHt;   
%  Angular frequency, k = -omega/aO line & k = omega/uC line
    kNegAcous = -1*stH*Uj*2*pi/aAmb;     kConvc = stH*Uj*2*pi/(uConv*Uj);    kPosAcous = kNegAcous*-1;          
%  Overlap definition & PSD initialization
    nOvlp = 0.5;        schPSD = zeros(size(inpSignal,1),blockSize);
    nBlks = (size(Master_U2,3)- (blockSize*nOvlp))/(blockSize - (blockSize*nOvlp)); 
    nOvlpblcks = blockSize*nOvlp;
%  FFTN computation Loop
    for iBlk = 1:nBlks
        offset = min((iBlk-1)*(blockSize - nOvlpblcks)+blockSize,size(inpSignal,2))-blockSize;
        sampleNo = (1:blockSize) + offset;                                      % adding offset to get correct overlap
        blk = inpSignal(:,sampleNo(1):sampleNo(end));           
        resFFT = fftn(blk);                resFFT = fftshift(resFFT,1);    % Centering wavenumbers to zero in middle
        resFFT = fliplr(resFFT);
        schPSD  = abs(resFFT).^2+schPSD;   
    end
%  Normalizing & re-scaling PSD
    schPSD = schPSD/((1/dt)*blockSize*nBlks);    schPSD = 10*log10(schPSD*delF/(20e-10)^2); 
%  Saving all PSDs in a main variable
    mainPSD(:,:,ctr) = schPSD;   dataLine = smoothdata(schPSD(:,freqPos),'gaussian',2);              %#ok
% Plotting screech frequency radial trend
    figure(fig12.Number);   plot(kRangeH,dataLine,'-o','LineWidth',1.5,'DisplayName',['$',num2str(radRange(ctr)),xName(4:6),'$'],'MarkerFaceColor','auto');
    hold on;   grid on;  
end
psdVals = linspace(200,260,size(schPSD(:,freqPos),1));   kIn1Line = zeros(size(kRangeH));   kIn1Line(:) = kIn1; 
figure(fig12.Number);   plot(kIn1Line,psdVals,'-k','LineWidth',1.5,'DisplayName','$k_{KH}-k_{S1}$');   lG = legend;   lG.Interpreter = 'latex';   lG.EdgeColor = [1 1 1];
box off;   ylim([200 260]);   xlim([-10 10]);   ax = gca;   ax.TickLabelInterpreter = 'latex';   ax.FontSize = 12;     xlabel('$kh_j$','Interpreter','latex'); 
ylabel('$\Sigma |\hat{\psi_f}|^2$','Interpreter','latex');   set(gcf,'Position',[100 340 880 450]);
%% 22.  PLOTTING CONTOUR & SELECTING CANDIDATE G-JM WAVE NO.  
%  Plotting contour
GF_FigurePlot(mainPSD(:,1:nCut,1)',kRangeH,stH); fig12 = gcf;      hold on;   ax = gca;   axis normal;    caxis([200 240]);   
xlabel('$kh_j$','Interpreter','latex');       xlim([-15 15]);   box off;   ax.FontSize = 13;           colormap(jet);   colorbar off; 
ylabel('$St_{h_j}$','Interpreter','latex');   ylim([0 0.7]);    set(gcf,'Position',[370 350 670 415]); ax.TickLabelInterpreter = 'latex';  
title(['$Y_{\omega,k} = ',num2str(radRange(1)),xName(4:6),'$'],'Interpreter','latex');  grid on;  set(gca,'Layer','top');          
%  Limit Strouhal number lines & Convective velocity lines
plot(kNegAcous,stH,'--k','LineWidth',1.2);  plot(kConvc,stH,'-.k','LineWidth',1.2);   plot(kPosAcous,stH,':k','LineWidth',1.2);
text(-9,0.4,'$k= -\omega/ a_0$','Color','k','FontSize',15,'Interpreter','latex');     text(6,0.4 ,'$k= \omega/ u_c$','Color','k','FontSize',15,'Interpreter','latex');
%  Selecting candidate G-JM frequency
figure(fig12.Number);   [x,y] = ginput(1);   [~,kLoc] = find(kRangeH < x);   kLoc = kLoc(end)+1;   kG = kRangeH(kLoc);
disp([newline,'==> Candidate G-JM Wave no.: ',num2str(kG)]);
%  Placing Marker at the loaction of G-JM
figure(fig12.Number);   hold on;   set(gca,'Layer','top');   plot(kG,stH(freqPos-8),'^','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','c');
%% 23.  PLOTTING DISPERSION RELATION PLOTS
% List of dispersion lines - Symmetric
dLines = {'S4.DAT'}; 
fid = fopen([driveOut,'Dispersion-Symmetric\',dLines{1}]);   dispData = textscan(fid,'%n%n','Delimiter',' ');
kVals = dispData{1};   stVals = dispData{2};  figure(fig12.Number);   plot(kVals,stVals,'b','LineWidth',1.2);
%%     * => 23.1 FIGURE SAVE: SPATIO-TEMPORAL FFT
% figName = ['SpatialFFT',tag,'-','RadialLoc-Y=',num2str((round(Yn(locTn),2)))];             GF_FigureSave(figName,driveOut,fig11.Number);
%  Saving to Journal draft folder
figName = ['STFFT-',fNameCon,'-',nozzle,'-NPR_',strrep(num2str(NPR),'.','p')];              GF_FigureSave(figName,driveManu,fig12.Number);
%%     * => 23.2 FIGURE SAVE: SCHLIEREN WITH LOCATION MARKER
% figName = ['RadialLocatns',num2str(radLoc),tag,];             GF_FigureSave(figName,driveOut,fig11.Number);
%  Saving to Journal draft folder
figName = ['STFFT-Location-',fNameCon,'-',nozzle,'-NPR_',strrep(num2str(NPR),'.','p')];           GF_FigureSave(figName,driveManu,fig11.Number);
%% 24.* STREAM WISE SPATIAL FOURIER TRANSFORM AT SELECTED FREQUENCY
%  Section computes the stream wise Fourier transform to convert spatially resolved images to spatial wavenumber domain 
%  Execute Sec. 19 first, Note: Factor for imgFFT*val for Rectan:Minor NPR 3.0(5),Min-2.6(3),Maj-2.6(3);2.5(1.2)
sfftIn.Xn = xNew;       sfftIn.inpMat = imgFFT(:,xStrt:end,freqPos)  ;  sfftIn.scrFreq = schFreq(freqPos); 
sfftIn.Yn = Yn;         sfftIn.lenScales = lenScales*1000;              sfftIn.NF = NF;
sfftIn.xName = xName;   sfftIn.yName = yName;                           sfftIn.config = config;
[sfftOut] = GF_SFFTCalc(sfftIn);
% Plot Lip Lines
figNo = sfftOut.figNo;   colr = 'w';   xPoints = -20:20;         lipLineLoc = zeros(2,size(xPoints,2));   lipLineLoc(1,:) = lipLine(1);   lipLineLoc(2,:) = lipLine(2);
figure(figNo);           ax = gca;     ax.ColorScale = 'log';    caxis([0.04 1.05]);   hold on;   set(gca,'Layer','top');    
plot(xPoints,lipLineLoc(1,:),':','LineWidth',1.2,'Color',colr,'HandleVisibility','off');   
plot(xPoints,lipLineLoc(2,:),':','LineWidth',1.2,'Color',colr,'HandleVisibility','off');   ax.TickLength = [0 0];   
%  Plotting kKH, kS1, kS2, kIn1, kIn2 & kG line
kKhLines = zeros(size(Yn));   kS1Line = zeros(size(Yn));   kS2Line = zeros(size(Yn));   
kKhLines(:) = kKh;            kS1Line(:) = kS1;            kS2Line(:) = kS2;        
kIn1Line = zeros(size(Yn));   kGLine = zeros(size(Yn));    kIn2Line = zeros(size(Yn));
kIn1Line(:) = kIn1;           kIn2Line(:) = kIn2; 
plot(kKhLines,Yn,':y','LineWidth',1.2,'DisplayName','$k_{KH}$');   plot(kS1Line,Yn,':m','LineWidth',1.2,'DisplayName','$k_{S1}$'); 
plot(kS2Line,Yn,':c','LineWidth',1.2,'DisplayName','$k_{S2}$');    plot(kIn1Line,Yn,'r','LineWidth',1.2,'DisplayName','$k_{KH}-k_{S1}$');     
plot(kIn2Line,Yn,'m','LineWidth',1.2,'DisplayName','$k_{KH}-k_{S2}$');    
lG = legend;   lG.Interpreter = 'latex';   lG.TextColor = 'w';   lG.FontSize = 13;   lG.Color = 'k';   lG.String(3) = {'$k_{\pm a}$'};        
%  Defining graphics object & disabling HandleVisibility for pcolor contour
A = gca;   depObj = A.Children;   depObj(end).HandleVisibility = 'off';   depObj(end-1).HandleVisibility = 'off';   lG.Position = [0.13 0.65 0.14 0.27];
%%    * => 24.1 FIGURE SAVE: STREAMWISE FFT
figName = ['SFFT-',fNameCon,'-',nozzle,'-NPR_',strrep(num2str(NPR),'.','p'),'-',num2str(schFreq(freqPos)),'Hz'];   
GF_FigureSave(figName,driveManu,figNo);   commandwindow;   GF_FigureSave(figName,driveOut,figNo);   
%% 25.* STREAMWISE FFT USING pspectrum
sFFTpS.Xn = Xn;   sFFTpS.xName = xName;   sFFTpS.inpMat = imgFFT(:,:,freqPos);   sFFTpS.NF = NF;          
sFFTpS.Yn = Yn;   sFFTpS.yName = yName;   sFFTpS.scrFreq = schFreq(freqPos);     sFFTpS.lenScales = lenScales;                                   
[sTpSOut] = GF_SFFTCalcPspec(sFFTpS);   
% Plot Lip Lines
figNo1 = sTpSOut.figNo;   colr = 'w';   xPoints = -20:20;         lipLineLoc = zeros(2,size(xPoints,2));   lipLineLoc(1,:) = lipLine(1);   lipLineLoc(2,:) = lipLine(2);
figure(figNo1);           ax = gca;     ax.ColorScale = 'log';    caxis([0.001 1]);   hold on;   set(gca,'Layer','top');    
plot(xPoints,lipLineLoc(1,:),':','LineWidth',1.2,'Color',colr,'HandleVisibility','off');   
plot(xPoints,lipLineLoc(2,:),':','LineWidth',1.2,'Color',colr,'HandleVisibility','off');   ax.TickLength = [0 0]; 
%  Plotting kKH, kIn1, kIn2 & kG line
kKhLines = zeros(size(Yn));   kS1Line = zeros(size(Yn));   kS2Line = zeros(size(Yn));
kKhLines(:) = kKh;            kS1Line(:) = kS1;            kS2Line(:) = kS2;
kIn1Line = zeros(size(Yn));   kIn2Line = zeros(size(Yn));   
kIn1Line(:) = kIn1;           kIn2Line(:) = kIn2;           
plot(kKhLines,Yn,':y','LineWidth',1.2,'DisplayName','$k_{KH}$');   plot(kS1Line,Yn,':m','LineWidth',1.2,'DisplayName','$k_{S1}$');         
plot(kS2Line,Yn,':c','LineWidth',1.2,'DisplayName','$k_{S2}$');    plot(kIn1Line,Yn,'r','LineWidth',1.2,'DisplayName','$k_{KH}-k_{S1}$');   
plot(kIn2Line,Yn,'m','LineWidth',1.2,'DisplayName','$k_{KH}-k_{S2}$');      
lG = legend;   lG.Interpreter = 'latex';   lG.TextColor = 'w';   lG.FontSize = 13;   lG.Color = 'k';   lG.String(3) = {'$k_{\pm a}$'};   lG.Position = [0.72 0.72 0.1 0.2];     
%  Defining graphics object & disabling HandleVisibility for pcolor contour
A = gca;    depObj = A.Children;   depObj(end).HandleVisibility = 'off';   depObj(end-1).HandleVisibility = 'off';   lG.Position = [0.13 0.65 0.14 0.27]; 
%%    * => 25.1 FIGURE SAVE: STREAMWISE FFT USING pspectrum
figName = ['SFFTpS-',fNameCon,'-',nozzle,'-NPR_',strrep(num2str(NPR),'.','p'),'-',num2str(schFreq(freqPos)),'Hz'];    
GF_FigureSave(figName,driveManu,figNo1);   commandwindow;   GF_FigureSave(figName,driveOut,figNo1);   
%% 25.* MODE PLOT: STREAM WISE SPATIAL FOURIER TRANSFORM AT SELECTED WAVENUMBER
sfftPlot.X = sfftOut.xNew;    sfftPlot.dwStrmRes = sfftOut.dwStrmRes;   sfftPlot.upStrmRes = sfftOut.upStrmRes;
sfftPlot.Y = sfftOut.yNew;    sfftPlot.baseFFT = sfftOut.baseFFT;       sfftPlot.lenScales = lenScales;
sfftPlot.kD = sfftOut.kD;     sfftPlot.figNo = sfftOut.figNo;           sfftPlot.NF = NF;
sfftPlot.xName = xName;       sfftPlot.yName = yName;                   sfftPlot.fPeak = schFreq(freqPos); 
sfftPlot.Uj    = Uj;          [kVal] = GF_SFFTModePlot(sfftPlot);       curntFig = gcf;
disp([newline,'=> Central Wave no.: ',num2str(kVal)]);
%% 26.* PLOT LIP-LINES ON SFFT MODES & ADJUST AXIS TO HALF WIDTH
if kVal < 0
   cmp = 'u';
else
   cmp = 'd';
end;   prompt = [newline,'=> Mode(KH/a_c/G): '];   modTyp = input(prompt,'s');
if strcmp(modTyp,'KH') 
   lipStyl = '-.w';
else
   lipStyl = '-.k';
end
xPoints = -20:20;   lipLineLoc = zeros(2,size(xPoints,2));   lipLineLoc(1,:) = lipLine(1);   lipLineLoc(2,:) = lipLine(2);   
figure(curntFig.Number);   caxis([-40 40]);   hold on;    set(gca,'Layer','top');   plot(xPoints,lipLineLoc,lipStyl,'LineWidth',1.2);  
colormap bluewhitered;     ylim([-2 2]);      ax = gca;   ax.TickLength = [0 0];    cB = colorbar;   cB.TickLabelInterpreter = 'latex';   
cL = title(cB,['$R_e[I_{fft}(\hat{\psi}_',cmp,')]$'],'Interpreter','latex');   set(gcf,'Position',[100 470 720 400]);  %set(gcf,'Position',[100 470 735 270]);      
cL.Rotation = 90;   cL.Position = [55 120 0];    
%%    * => 26.1 FIGURE SAVE: STREAMWISE FFT MODE PLOT
figName = ['SMP-',fNameCon,'-',nozzle,'-NPR_',strrep(num2str(NPR),'.','p'),'-',num2str(schFreq(freqPos)),'Hz-',modTyp];    
GF_FigureSave(figName,driveManu,curntFig.Number);   commandwindow;   GF_FigureSave(figName,driveOut,curntFig.Number);   

%% -------------------------------------------- Functions ---------------------------------------------------------------%%
%%      Function 1: COMPUTE BLOCKSIZE BASED ON ACQUISITION RATE
function [blkSize,resIntnstySclr] = blkSizer(config)
    nozlConfig     = {'C'; 'S'; 'S2';'SS2'; 'TR';'TR2';'TS';'TS1';'TS2';'TRV0'};
    frameRate      = [45E3;41E3;204E3;204E3;41E3;204E3;41E3; 112E3;204E3; 45E3];
    fftBlks        = [900; 820; 2048; 2048; 820; 2048; 820; 2240; 2048;  900];
    intstyScalr    = [11.2;11.5;11.5;11.5;  11.5;11.5; 11.5;11.5; 11.5;  11.5]; 
    schTestConfig  = table(nozlConfig,frameRate,fftBlks,intstyScalr);
    inDx           = strcmp(schTestConfig.nozlConfig,config);
    blkSize        = schTestConfig.fftBlks(inDx);
    resIntnstySclr = schTestConfig.intstyScalr(inDx);
end
%%      Function 2: FIND PRECISE LOCATION
function [xCol,xVal,yRow,yVal] = precisLoc(inpX,inpY,Xn,Yn)
%  Finding the precise X location match for input value
    [~,x] = find(Xn == inpX);          
    if isempty(x)
       [~,preX] = find(Xn<inpX);       [~,postX] = find(Xn>inpX);
       preDif = inpX - Xn(preX(end));  postDif = Xn(postX(1)) - inpX;
       if preDif < postDif
          xCol = preX(end);   xVal = round(Xn(preX(end)),2);
       else
          xCol = postX(1);    xVal = round(Xn(postX(1)),2);
       end
    else
       xCol = x;   xVal = Xn(x);
    end
%  Finding the precise X location match for input value
    [y,~] = find(Yn == inpY);
    if isempty(y)
       [preY,~] = find(Yn<inpY);       [postY] = find(Yn>inpY);
       preDif = inpY - Yn(preY(end));  postDif = Yn(postY(1)) - inpY;
       if preDif < postDif
          yRow = preY(end);   yVal = round(Yn(preY(end)),2);
       else
          yRow = postY(1);    yVal = round(Yn(postY(1)),2);
       end
    else
       yRow = y;   yVal = Yn(y);
    end
end
%%      Function 3: DRIVE CHECK FOR VALID INPUT DRIVE
function driveCheck(driveIn)
    if isfolder(driveIn) == 0
       disp([newline,'->> Check nozzle & condition']); return;
    end
end
%%      Function 4: FUNCTION FOR CORRECTING BACKGROUND SUBTRACTION FOR TS-MINOR
function [avgOne,M2Fluc] = bckgCorctnOne(Master_U,M2,Bckgrnd)
% Compute true mean & background subtracted mean
    trueMean = mean(Master_U,3);             bkSubMean1 = mean(M2,3);
% Shifting masked region in background to avoid nozzle exit overlap in correction
    shftdBckgrnd = circshift(mean(Bckgrnd,3),-12,2); 
% Subtracting shifted background & new average image
    bkSubMean2 = trueMean./shftdBckgrnd;      avgOne = zeros(size(bkSubMean1));
% Taking parts of non-overlaped nozzle exit from shifted correction(colLocation) & 
% rest of the flow from original background subtracted image
    colLoctn = 62;                            avgOne(:,1:colLoctn) = bkSubMean2(:,1:colLoctn);     
    avgOne(:,colLoctn:end) = bkSubMean1(:,colLoctn:end);
% Recomputing background subtraction with shifted background
    M2 = Master_U./shftdBckgrnd;              M2Fluc = M2 - mean(M2,3);
end
%%      Function 5: CONVERT IMAGES TO TECPLOT FORMAT FOR PLOTTING
function tecWrite(writeInp)
% Adding code drive to path & changing drive to output folder
    addpath([writeInp.codeRoot 'Jet_Analysis\Global_Functions\']);      cd(writeInp.driveOut); 
% Defining meshed axes
    [X1,Y1] = meshgrid(writeInp.X,writeInp.Y);
% Input parameters based on nozzle config
    nozType = {'Major';'Minor'};                 zoneName = {'CDSide(MajorV)';'FlatSide(MinorV)'};
    axisOrder = [2,3];                           inDx1 = strcmp(writeInp.nozzle,nozType);
% Input parameters based on contour choice
    varNames = {'PhaseDistr'; 'EnergyDistr'};    inDx2 = writeInp.chc; 
% Holding Struct to write to Tecplot format
    holdMat = [];    holdMat.Nvar = 4;                              % Number of variables to store
    holdMat.varnames = {'x','y','z',cell2mat(varNames(inDx2))};     % Variable Names
    holdMat.surfaces(1).zonename = cell2mat(zoneName(inDx1));       % Contour Zone
    holdMat.surfaces(1).x = X1;                                     % Axial Co-ordinate
    if strcmp(writeInp.nozzle,'Major')                              % Radial Co-ordinate
       holdMat.surfaces(1).z = Y1;
    else
       holdMat.surfaces(1).y = Y1;
    end
    holdMat.surfaces(1).order = axisOrder(inDx1);                   % Axis order: 3(default) - xy; 2 - xz;
    holdMat.surfaces(1).v(1,:,:) = writeInp.contour;                % Contour Values(uBar,vBar,etc.)
    contourName = writeInp.figName;
% Function to export Tecplot format
    mat2tecplot(holdMat,contourName);
% Changing back to code root
    cd([writeInp.codeRoot 'Jet_Analysis\Global_Functions\']);      clc;
end
%%      Function 6: RADIAL DISTANCES TO FIND PEAKS IN RESOLVED INTENSITY & PHASE DISTRIBUTION
function [d1,d2,d3,d4,d5] = radDistAlloc(config,NF)
    if strcmp(config,'C')
       d1 = 0.5;    d2 = 0.6;    d3 = 0.65;   d4 = 0.7;   d5 = 1;
    elseif strcmp(config,'S') && strcmp(NF,'D') 
       d1 = 0.34;   d2 = 0.5;    d3 = 0.6;    d4 = 0.7;   d5 = 1;
    end
end
%%      Function 7: PLOTTING SHOCK LOCATIONS FROM SCHLIEREN
function shckLocPlot(driveOut,NF,lenScales,Xn,Yn,fig1,C)
%  Loading the location file from drive
    fid = fopen([driveOut,'ShockCellLoc.DAT']);   data = textscan(fid,'%n');   shockLoc = data{1};
%  Normalizing based on specified length scale(NF)
    if strcmp(NF,'D')
       shockLoc = shockLoc/(lenScales(1)*1000);
    else
       shockLoc = shockLoc/(lenScales(2)*1000);
    end
%  Plotting vertical lines to show shock cell locations
    vertLines = zeros(size(Yn,1),1); 
    for ctr = 1:size(shockLoc,1)
        [~,pos] = find(Xn<shockLoc(ctr));   vertLines(:) = Xn(pos(end)+1);   figure(fig1.Number);
        hold on;   set(gca,'Layer','top');  plot(vertLines,Yn,':','Color',C,'LineWidth',1.05);   ax = gca;   ax.TickLength = [0 0];
    end
end
%%      Function 8: GENERATING IMAGE SAVE TAGS & LIP LINE LOCATION
function [tag,lipLine,fNameCon] = tagNLine(NF,config,nozzle,lenScales)
    if strcmp(NF,'D') && strcmp(config,'C')
       tag = '(D)';   lipLine = ((lenScales(1)/2)/lenScales(1));       lipLine = [-lipLine lipLine];   fNameCon = 'Circ';
    elseif strcmp(NF,'H') && strcmp(config,'C')
           tag = '(h)';   lipLine = ((lenScales(1)/2)/lenScales(1));   lipLine = [-lipLine lipLine];   fNameCon = 'Circ';
    elseif strcmp(NF,'D') && strcmp(config,'C')~=1 
           tag = '(D)';
           if strcmp(nozzle,'Minor')
              lipLine = ((lenScales(2)/2)/lenScales(1));
           else
              lipLine = ((lenScales(2))/lenScales(1));
           end;   lipLine = [-lipLine lipLine];
    elseif strcmp(NF,'H') && strcmp(config,'C')~=1 
           tag = '(h)';   
           if strcmp(nozzle,'Minor')
              lipLine = ((lenScales(2)/2)/lenScales(2));
           else
              lipLine = ((lenScales(2))/lenScales(2));
           end;   lipLine = [-lipLine lipLine]; 
    end 
    if strcmp(config,'S')
       fNameCon = 'Rectan';
    elseif strcmp(config,'TR')
       fNameCon = 'TwinRect';
    elseif strcmp(config,'TS')
       fNameCon = 'TwinSqur';
    end
end

%%      Function 9: PLOTTING SHOCK LOCATIONS FROM PIV
% function plotShockLoc(figNo,config,NPR,NF,Yn)
%     if strcmp(config,'C') && NPR == 2.5 && strcmp(NF,'D')
%        prompt = [newline 'Plot Shock Locations(y/n) - '];    chc = input(prompt,'s');
%        shckLoc = [0.23278,1.00891,1.7188,2.36698,2.98428,3.57072,4.12629,4.62014,5.08312];
%        if strcmp(chc,'y')
%           figure(figNo);   aX = gca;   set(aX,'Layer','top');   aX.TickLength = [0.001,0.001]; 
%           for ctr = 1:length(shckLoc)
%               sL = ones(size(Yn))*shckLoc(ctr);
%               plot(sL,Yn,'-.k','LineWidth',1.2);   hold on;
%           end
%        end
%     end
% end
%% Appendix: Deprecated 
% STREAM WISE SPATIAL FOURIER TRANSFORM ALONG AVERAGE SCHLIEREN CENTERLINE
%  Section - 1: Compute the streamwise FFT along the centerline of the avg Schlieren image
% Defining Extent of axial domain
% Code Start
    % prompt = [newline,'Final Axial Location(0-',num2str(Xn(end)),'): '];   xEnd = input(prompt);
    % [~,xEnd] = find(Xn<xEnd);   endPos = xEnd(end)+1;   Deq = lenScales(1)*1000;
    % %  Making the number of axial locations odd numbered 
    % if rem(size(xStrt:endPos,2),2)
    %    axSpan = xStrt:endPos;       % Odd numbered
    % else
    %    axSpan = xStrt:endPos - 1;  % Changing Even numbered span to Odd
    % end
    % %  Defining the minimum & maximum wave numbers based on the resolution &
    % %  size of the axial domain(X in mm), number of points in the domain, domain resolution(m)
    % axDom = Xn(axSpan)*Deq;   sigPts = size(axDom,2);   dx = (axDom(3) - axDom(2))/1000;  
    % %  Wavenumber Range computation
    % kNos = GF_FourierFreq(sigPts,2*pi/dx);   [~,neg] = find(kNos<0);   [~,pos] = find(kNos>0);   kRange = zeros(size(kNos));
    % negVals = kNos(neg);   posVals = kNos(pos);  kRange(neg) = posVals;   kRange(pos-1) = negVals;   kRangeD = kRange*Deq/1000;
    % %  Spatial Fourier transform along the centerline of average Schlieren image
    % %  Applying fftshift & normalizing with peak absolute value
    % cYLoc = find(Yn==0);  cSig = avgOne(cYLoc,axSpan);   cSig   = cSig - mean(cSig);   
    % cFT  = fft(cSig);     cFT  = fftshift(cFT);         cFTMax = max(abs(cFT));     cFT = abs(cFT)/cFTMax; 
    % %  Plotting the spatial FFT result along the centerline vs axial wavenumber
    % fig12 = figure(12);   plot(kRangeD,cFT,'-.r','LineWidth',1.2);   grid on;   box off;   ax = gca;   ax.TickLabelInterpreter = 'latex';   hold on;   ax.FontSize = 12;
    % xlim([0 20]);   ylim([0 1.05]);   xlabel('$k_xD$','Interpreter','latex');   ylabel('$|\hat{I}|$','Interpreter','latex');  set(gcf,'Position',[70 400 640 350]);
% Code End
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%  SPL - CONTOUR
% Code Start
    % fig3 = figure; contourf(Xn,Yn,imgSpec(:,:,pos),40), shading interp, axis equal, colormap(schJet);
    % caxis(cLim); xlim(lim_x), ylim([0 2.5]);  xlabel(xname,'Interpreter','latex'); ylabel(yname,'Interpreter','latex');
    % title(['$Acoustic \thinspace Map: St = ',num2str(St(pos)),'$'],'Interpreter','latex');set(gcf,'Position',[180 300 740 360]);
    % c = colorbar; title(c,'\it\epsilon_{f_s}');   set(gca,'FontSize',12,'FontName','times'); 
%-------------------------------------- Spectrogram trial ------------------------------------------
    % [blockSize] = blkSizer(config);     pixlCtr = 1;    freqRes = (1/dt)/blockSize;    freqRng = 0:freqRes:(1/dt)/2;
    % fftMain = zeros(size(freqRng,2),round(nBlks),rowNos*colNos);                       wB = waitbar(0,'Computing Spectrogram'); 
    % for row = 1:rowNos
    %   for col = 1:colNos
    %       inp = reshape(Master_U2(row,col,:),[1,imgNos]);       
    %       [fftRes,~,tR,~] = spectrogram(inp,blockSize,blockSize/2,blockSize,1/dt);  waitbar(row/rowNos,wB);
    %       fftMain(:,:,pixlCtr) = fftRes;                pixlCtr = pixlCtr + 1;     
    %   end
    % end;  close(wB); 
    % for ctr = 1:3
    %   imgFFT = fftMain(pos,ctr,:);          imgFFT = reshape(imgFFT,colNos,rowNos);         imgPSD = abs(imgFFT).^2;   
    %   imgPSD = imgPSD*2;                  imgSPL = 10*log10(imgPSD*50/(pRef^2));          GF_FigurePlot(imgSPL');     GF_FigurePlot(-angle(imgFFT)');
    % end
%------------------- CALCULATING SPL AT EACH LOCATION FOR THE WHOLE IMAGE SET -----------------------
% REFERENCE PRESSURE & RESULTING FREQUENCY MATRIX 
    % [blockSize] = blkSizer(config);                         bckFluc = Bckgrnd - bckMean;                              
    % pRef    = mean(mean(bckMean))*6E-10;                    % Using background image as basline value & mulitplying by factor to enhance variations   
    % nOvlp = 0.5;                                            % percentage of overlap
    % nBlks   = (size(Master_U2,3)- (blockSize*nOvlp))/(blockSize - (blockSize*nOvlp)); 
    % schFreq = 0:1/dt/blockSize:1/dt-1;                      % frequency spectrum
    % % nFreq    = size(freq,2)/2+1;                          % Nyquist cut-off limit for frequency
    % nFreq   = size(schFreq,2);                              % Use if not computing single sided spectrum
    % schFreq    = schFreq(1,1:nFreq);                        % truncating frequencies to Nyquist limit
    % imgSpec = zeros(size(Master_U2,1),size(Master_U2,2),nFreq);    
    % imgPhs  = zeros(size(Master_U2,1),size(Master_U2,2),nFreq);
    % imgFFT  = zeros(size(Master_U2,1),size(Master_U2,2),nFreq);
    % schPSD     = zeros(size(Master_U2,1),size(Master_U2,2),nFreq);
    % OASPL   = zeros(size(Master_U2,1),size(Master_U2,2));        wB = waitbar(0,'Computing Image FFT');                     
    % for row = 1:rowNos
    %   for col = 1:colNos
    %     inp = reshape(Master_U2(row,col,:),[1,imgNos]);         
    %     [imgSpec(row,col,:),imgPhs(row,col,:),schPSD(row,col,:),OASPL(row,col),imgFFT(row,col,:)] = Sch_Img_FFT_v1(inp,blockSize,nOvlp,schFreq,1/dt,pRef); waitbar(row/rowNos,wB);   
    %   end
    % end;    close(wB);    maxX = round(max(Xn),2);     minY = round(min(Yn),2);     maxY = round(max(Yn),2);
    % disp([newline '-------------------------------------------']);
    % disp(['-> Number of blocks     :- ',num2str(floor(nBlks))]);
    % disp(['-> Frequency resolution :- ',num2str(schFreq(2))]);
    % disp('-------------------------------------------');
% Code End
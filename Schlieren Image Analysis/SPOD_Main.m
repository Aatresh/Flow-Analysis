%% SPOD ANALYSIS FROM SCHLIEREN & SHADOWGRAPH IMAGES 
%% WINDOW SIZE - JASA
% CONTOUR1 - set(gcf,'Position',[70 296 574 332]); MODE -  set(gcf,'Position',[328 395 441 324]);
% CONTOUR WITH COLORBAR - set(gcf,'Position',[80  310  630  330]);
% CONTOUR2 - ylim([-6 6]); xlim([0 12]); set(gcf,'Position',[137 290 413 407]);
% xlim([0 10]);   xticks([0 0.5 2 4 6 8 10])
% POD MODE - ylim([-1.5 1.5]); xticks([0.5 1 2 3 4 5]); set(gcf,'Position',[149 332 617 611]);
% RECONSTRUCTED VELOCITY (CONVEC, PHASED, STREAMLINES) - ylim([-1.2 1.2]); set(gcf,'Position',[164 395 643 308]);
% caxis([0 0.01]); xticks([0.5 1 2 3 4 5]);
% STREAMLINES - set(gcf,'Position',[79 372 467 399]);
% NEW WINDOW SIZES: 
% STD DEV - MAJOR(Deq): ylim([-2.5 2.5]); set(gcf,'Position',[655 340 639 417]); set(gcf,'Position',[70 296 574 332]);
% set(gcf,'Position',[105 310 560 310]);                                

clearvars;   clc;   set(0,'defaultfigurecolor',[1 1 1]);  code = 'SPOD';
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';
addpath([root1 'Jet_Analysis\Schlieren & SPOD_Codes\']);   cd([root1 'Jet_Analysis\Global_Functions\']);   
tests = {'NPR_2p5_TR_1p0', 'NPR_2p6_TR_1p0', 'NPR_2p9_TR_1p0','NPR_3p0_TR_1p0',...
         'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_4p5_TR_1p0','NPR_5p0_TR_1p0', 'NPR_6p0_TR_1p0'};

condition = tests(4);   config = 'TS';    nozzle = 'Major';    NF = 'D';  

%% 1.  * DRIVE PATH SELECTION 
[OutputStruct] = GF_DriveSelect(config,nozzle,code);
inRoot = OutputStruct.in_root; outRoot = OutputStruct.out_root;  nozzle = OutputStruct.nozzle; dt = OutputStruct.dt; 
disp([newline '->> Acquisition Rate - ',num2str(1/dt),' Hz']);
%  Figure tags
if strcmp(NF,'D')
   tag = '(D)';
else
   tag = '(h)';
end
%  Colormaps 
load custom_map3; 	 load spod_custom2.mat;       
%  Jet flow parameters
[Mj,Uj,NPR,NTR] = GF_Velocity(condition{1});
%  Loading case .mat files
driveIn  = [inRoot condition{1}(9:14) '\' condition{1}(1:7) '\'];   driveOut = [outRoot condition{1}(9:14) '\' condition{1}(1:7) '\'];  
load([driveIn condition{1}(1:7) '_DAT']);   load([driveIn 'X']);    load([driveIn 'Y']);
%  Loading Background file
if isempty(OutputStruct.bckgrnd_root) ~= 1 
   Bckgrnd = load([OutputStruct.bckgrnd_root 'Bckgrnd_DAT']);    bckLim = [0.1 1.8];
   Bckgrnd = Bckgrnd.Master_U;  M2 = Master_U./mean(Bckgrnd,3);  avgOne = mean(M2,3);   avgTwo = mean(Master_U,3);     
else
   Bckgrnd = [];  bckLim = [20 400];  M2 = Master_U;
end;   disp([newline '-> Done!']);   beep;
%% 2.  * AXIS NORMALIZATION & STANDARD DEVIATION
[Xn,Yn,limX,limY,xName,yName,lenScales,fig_siz] = GF_AxisDefnSch(config,nozzle,NF,X,Y);
%  Standard deviation calculation & plot
rho_std = std(M2,[],3);   rho_std = log10(rho_std);   fig = figure; pcolor(Xn,Yn,rho_std);     shading interp;   colormap jet;   axis equal;
caxis([-1.45 -0.2]);      yticks([-2 -1 0 1 2]);      set(gcf,'Position',[490 440 560 420]);   xlim(limX),       ylim(limY);     xlabel(xName,'Interpreter','latex'); 
ylabel(yName,'Interpreter','latex');   ax = gca;      ax.TickLabelInterpreter = 'latex';       ax.FontSize = 12; 
%  Matrix dimensions
rowNos = size(M2,1);   colNos = size(M2,2);   imgNos = size(M2,3);   Yn = Yn';   Xn = Xn';
%%  =>   2.1 FIGURE SAVE: STANDARD DEVIATION
figName = ['Std_Devtn_hot' tag '-',condition{1}];   GF_FigureSave(figName,driveOut,fig.Number);   clear figName;
%% 3.    PLOTTING INSTANTANEUOUS IMAGES
[figNo,imgNo] = plotInstantaneous(M2,Xn,Yn,limX,limY,xName,yName,bckLim); 
%%  =>   3.2 FIGURE SAVE: INSTANTANEOUS IMAGES
figName = ['Instant-imgNo-',num2str(imgNo),tag,'-',condition{1}];   GF_FigureSave(figName,driveOut,figNo);
%% 4.  * CHANGING IMAGE MATRIX ORDER TO TIME AS FIRST DIMENSION(NO. OF SNAPSHOTS)
Master_U = permute(Master_U,[3 1 2]);       M2 = permute(M2,[3 1 2]);
%% 5.  * STRUCT FORMATTING FOR CONTINUOUS IMAGE LOOPING
SPOD_INP.nozzle = nozzle;   SPOD_INP.lim_x = limX;                SPOD_INP.lim_y    = limY;
SPOD_INP.TYP = 'SC';        SPOD_INP.Sz    = [490 440 560 420];   SPOD_INP.driveOut = driveOut;
SPOD_INP.Xn  = Xn;           SPOD_INP.Yn   = Yn;                  SPOD_INP.NF       = NF;
SPOD_INP.gifName = [condition{1} tag '.gif'];                     SPOD_INP.config   = config;
SPOD_INP.bckLim = bckLim;   % SPOD_INP.Sz = [676 88 1116 886];
%%  =>   5.1 .gif GEN/SAVE: INSTANTANEIOUS IMAGES 
GF_GIF_DisplayWrite(M2,SPOD_INP);
%% 6.  * CLEARING FIELDS TO PREP FOR SPOD RESULTS
SPOD_INP = rmfield(SPOD_INP,'TYP');   SPOD_INP = rmfield(SPOD_INP,'gifName');   clear M2;     
%% 7.  * SPOD CALCULATION
% INPUT :
% Master_U  - Image set with Time as 1st dimension
% ones(1,w) - No. of images to group - Defines the block size & no. of modes
% OVLP      - No. of images overlapping between blocks
% dt        - 1/fs
% n_fft     - No. of realizatoins per block
% OUTPUTS :
%    L      - Modal Energy Spectra of all Frequencies(Frequencies X Modes)
%    P      - SPatial SPOD Modes arranged in descending order
%    freq   - Frequencies 
[nFFT] = blkSizer(config);    delF = 1/(dt*nFFT);       nOvlp = nFFT*0.5;   tic;  
[L,P,spodFreq] = spod(Master_U,ones(1,nFFT),[],nOvlp,dt);   spodFreq = spodFreq';  St = spodFreq*lenScales(1)/Uj;   
%  Plotting Modal Energy spectra vs Frequency(Hz)
fig_L = figure('name','Modal Energy Spectra');   set(fig_L,'Position',[100 300 998 359]);  ylabel('$SPOD \thinspace Mode \thinspace Energy$','Interpreter','latex');
subplot(121);   loglog(spodFreq,L,'LineWidth',1.2);    grid on;   xlim([300 spodFreq(end)-200]);     xlabel('$Hz$','Interpreter','latex');   box off; 
ylim([min(min(L)),max(L(:,1))]); ax = gca;   ax.FontSize =12;   title(['$TR \thinspace ',num2str(NTR),'.0 - NPR \thinspace ',num2str(NPR),'$'],'Interpreter','latex');  ax.TickLabelInterpreter = 'latex';
%   Plotting Modal Energy spectra vs Strouhal number(St)
subplot(122);   loglog(St,L,'LineWidth',1.2);grid on;   xlabel('$St$','Interpreter','latex');   ylabel('$SPOD \thinspace Mode \thinspace Energy$','Interpreter','latex');   ax = gca; 
xlim([0.01 St(end)-0.02]);   box off;   ylim([min(min(L)),max(L(:,1))]);   ax.FontSize = 12;   title(['$TR \thinspace ',num2str(NTR),'.0 - NPR \thinspace ',num2str(NPR),'$'],'Interpreter','latex');   ax.TickLabelInterpreter = 'latex';
%   Display mode & frequency count
disp(['No. of Frequencies         : ',num2str(length(spodFreq))]);
disp(['No. of modes per frequency : ',num2str(size(P,4))]);  figName = ['SPOD_Spectra' tag '-',condition{1}]; commandwindow;
%% =>    7.1 FIGURE SAVEl: SPOD MODE & FREQUENCY SPECTRA
GF_FigureSave(figName,driveOut,fig_L.Number)
%% 8.  * RE-DEFINING PARAMETERS FOR .gif GENERATION
SPOD_INP.L = L;   SPOD_INP.Freq = spodFreq;   SPOD_INP.figname = ['SPOD_vs_Acousctics' tag '-',condition{1}];   toc;
%% 9.    COMPARING WITH PROXIMITY TEST DATA
SPOD_Acoustics_Cmp(config,nozzle,condition,SPOD_INP);
%% 10. * TAKE INPUT TO PLOT RESULTS FOR A GIVEN MODE & FREQUENCY
Mode = 1;   prompt = [newline '->> Enter Frequency - '];   freqVal = input(prompt);
[a,~] = find(spodFreq<=freqVal);   freqPos = a(end);   figName = ['SPOD-Freq', tag, num2str(spodFreq(freqPos)),' Mode_', num2str(Mode)];
%   Display frequency specific energy for selected mode
figure(freqPos);   clim = [-0.005 0.005];   % set(gcf,'Position',[69 1 1852 1003]);
pcolor(Xn,Yn,real(squeeze(P(freqPos,:,:,Mode))));     shading interp,   axis equal,   colormap(custom_map3);   caxis(clim) 
xlabel(xName,'Interpreter','latex');   xlim(limX),    ax = gca;   ax.FontSize = 12;   ax.TickLabelInterpreter = 'latex';
ylabel(yName,'Interpreter','latex');   ylim(limY);    title(['$Freq. = ' num2str(spodFreq(freqPos)) '; Mode = ',num2str(Mode) '; \lambda = ' num2str(L(freqPos,Mode),'%.4g'),'$'],'Interpreter','latex');
%% =>    10.1 FIGURE SAVE: MODE ENERGY DISTRIBUTION
GF_FigureSave(figName,driveOut,freqVal);
%% =>  * 10.2 GENERATING & SAVING SPOD MODE SHIFT.gifs
% Animation parameters nt: No. of Frames, T: Time period, time: Total duration
nt = 80;    T = 1/spodFreq(freqPos);   time = linspace(0,T,nt); 
%  Adding above parameters to input Struct file 
SPOD_INP.time = time;  SPOD_INP.nt   = nt;    SPOD_INP.freqNo = freqPos;   SPOD_INP.map  = custom_map3; 
SPOD_INP.clim = clim;  SPOD_INP.Mode = Mode;  SPOD_INP.St  = St;          SPOD_INP.TYP = code;  
SPOD_INP.gifName = [figName '.gif'];          
GF_GIF_DisplayWrite(P,SPOD_INP); commandwindow;
%% 11. * SPATIAL FOURIER TRANSFORM
%  Re-arranging axis to start from the nozzle exit
[~,xStrt] = find(Xn==0);   xNew = Xn(xStrt:end); 
%  Taking the abolute value of the mode for input to Spatial FFT
spodMat = permute(P(:,:,:,Mode),[2,3,1]);  
%  Input Struct for Spatial FFT
sfftIn.Xn = xNew;       sfftIn.inpMat = spodMat(:,xStrt:end,freqPos);   sfftIn.scrFreq = spodFreq(freqPos); 
sfftIn.Yn = Yn;         sfftIn.lenScales = lenScales*1000;              sfftIn.NF = NF;
sfftIn.xName = xName;   sfftIn.yName = yName;                           sfftIn.config = config;
[sfftOut] = GF_SFFTCalc(sfftIn);
%% 12. * STREAMWISE FFT USING pspectrum
sFFTpS.Xn = xNew; sFFTpS.xName = xName;   sFFTpS.inpMat = spodMat(:,xStrt:end,freqPos);   sFFTpS.NF = NF;          
sFFTpS.Yn = Yn;   sFFTpS.yName = yName;   sFFTpS.scrFreq = spodFreq(freqPos);     sFFTpS.lenScales = lenScales;                                   
[sTpSOut] = GF_SFFTCalcPspec(sFFTpS);     figNo = sTpSOut.figNo;
figure(figNo);           ax = gca;     ax.ColorScale = 'log';    caxis([0.001 1]);   hold on;   set(gca,'Layer','top');
%% 13. * MODE PLOT: STREAM WISE SPATIAL FOURIER TRANSFORM AT SELECTED FREQUENCY
sfftPlot.X = sfftOut.xNew;    sfftPlot.dwStrmRes = sfftOut.dwStrmRes;   sfftPlot.upStrmRes = sfftOut.upStrmRes;
sfftPlot.Y = sfftOut.yNew;    sfftPlot.baseFFT = sfftOut.baseFFT;       sfftPlot.lenScales = lenScales;
sfftPlot.kD = sfftOut.kD;     sfftPlot.figNo = sfftOut.figNo;           sfftPlot.NF = NF;
sfftPlot.xName = xName;       sfftPlot.yName = yName;                   sfftPlot.fPeak = spodFreq(freqPos); 
sfftPlot.Uj    = Uj;          [kVal] = GF_SFFTModePlot(sfftPlot);       curntFig = gcf;
disp([newline,'=> Central Wave no.: ',num2str(kVal)]);
%% -------------------------------------------- Functions ---------------------------------------------------------------%%
%%      Function 1: COMPUTE BLOCKSIZE BASED ON ACQUISITION RATE
function [blkSize] = blkSizer(config)
    nozlConfig = {'C';'S';'S2';'SS2';'TR';'TR2';'TS';'TS1';'TS2';'TRV0'};
    frameRate  = [45E3;41E3;204E3;204E3;41E3;204E3;41E3;112E3;204E3; 45E3];
    fftBlks    = [900; 820; 2048; 2048; 820; 2048; 820; 2240; 2048;  900];
    schTestConfig = table(nozlConfig,frameRate,fftBlks);
    inDx = strcmp(schTestConfig.nozlConfig,config);
    blkSize = schTestConfig.fftBlks(inDx);
end
%%      Function 2: INSTANTANEOUS IMAGES .gif
function [figNo,imgNo] = plotInstantaneous(M2,Xn,Yn,lim_x,lim_y,xName,yName,bckLim)
    prompt = [newline,'->> Pick or Randomize(p/r) - '];   imgOrdr = input(prompt,'s');
    chc = 'y';    load('schInst.mat');   %#ok
    while strcmp(chc,'y')
        if strcmp(imgOrdr,'p')
%  Selecting a particular image
           prompt = [newline '->> Choose figure number(1:',num2str(size(M2,3)),') - ']; imgNo = input(prompt);
        else
%  Random image
           imgNo = randperm(size(M2,3),1);  disp(['Image number - ',num2str(imgNo)]);        
        end
%  Plotting figure
        GF_FigurePlot(M2(:,:,imgNo),Xn,Yn); figNo = get(gcf,'Number'); figure(figNo);  
%  Formatting
       xlim(lim_x); ylim(lim_y); caxis(bckLim); colormap(gray); colorbar off; ax = gca;
       ax.TickLabelInterpreter = 'latex'; set(gca,'FontSize',12); set(gcf,'Position',[490 440 560 420]); 
       xlabel(xName,'Interpreter','latex'); ylabel(yName,'Interpreter','latex'); %yticks([-2 -1 0 1 2]);
       prompt = [newline '->> Plot another image?(y/n) - '];  chc = input(prompt,'s');
    end
end

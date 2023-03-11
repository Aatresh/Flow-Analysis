% DESCRIPTION:
%  Code uses Fourier decomposition to seperate the upstream & downstream 
%  propagating parts of supersonic jets from Schlieren images by using
%  n-dimensional Fourier transform
% ----------------------------- Configuration Legend: Active Use -----------------------------------
% C      - CircularMedium --------------------------> Acquisition Rate: 45,000  Hz
% S      - SingleRectangular -----------------------> Acquisition Rate: 41,000  Hz
% S2/SS2 - SingleRectangular(Schlieren/Shadowgraph)-> Acquisition Rate: 204,800 Hz
% TR     - TwinRectagular --------------------------> Acquisition Rate: 41,000  Hz
% TR2    - TwinRectangular -------------------------> Acquisition Rate: 204,800 Hz
% TS     - TwinSquare ------------------------------> Acquisition Rate: 41,000  Hz
% TS1    - TwinSquare ------------------------------> Acquisition Rate: 112,000 Hz
% TS2    - TwinSquare ------------------------------> Acquisition Rate: 204,800 Hz
% Nozzle heights: C - 20.6502;  S - 12.945; TR - 12.19;  TS - 16.61
%---------------------------------------------------------------------------------------------------

tic;        fclose all;     clc;        clearvars;     set(0,'defaultfigurecolor',[1 1 1]);       
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';     code = 'SchFourier';
cd([root1 'Jet_Analysis\Global_Functions\']); addpath([root1 'Jet_Analysis\Schlieren & SPOD_Codes\']);       
load 'blckToRed.mat';    load schJet.mat;    load 'schFaMap.mat';    load blckToRed.mat;
tests = {'NPR_2p5_TR_1p0', 'NPR_2p6_TR_1p0', 'NPR_2p9_TR_1p0','NPR_3p0_TR_1p0',...
         'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_5p0_TR_1p0', 'NPR_6p0_TR_1p0'};

condition = tests(4);   config = 'C';   nozzle = 'Major';   NF = 'D';     

% Drive selection based on frame rate 
   [OutputStruct] = GF_DriveSelect(config,nozzle,code);   nozzle = OutputStruct.nozzle;   dt = OutputStruct.dt;
% Jet Parameters
   [Mj,Uj,NPR,NTR] = GF_Velocity(condition{1});           disp([newline '  ->> Acquisition Rate - ',num2str(1/dt),' Hz']);
   spdSound = 345.501931;
% Loading video matrix
   driveIn  = [OutputStruct.in_root,condition{1}(9:14),'\',condition{1}(1:7)  '\'];
   driveOut = [OutputStruct.out_root,condition{1}(9:14) '\',condition{1}(1:7) '\'];
   driveCheck(driveIn);   load([driveIn,condition{1}(1:7),'_DAT']);   load([driveIn,'X']);   load([driveIn,'Y']);
%% Background Subtraction
   bckgrndPath = OutputStruct.bckgrnd_root;    [M2,~,~] = bckgrndSub(bckgrndPath,Master_U,config);
% Axis definition & normalization
   [Xn,Yn,limX,limY,xName,yName,lenScales,figSize] = GF_AxisDefnSch(config,nozzle,NF,X,Y); 
   Xn = Xn';          Yn = Yn';        nozHt = lenScales(1);         textSize = 13;    
%% Normalize images with respective peak values after nozzle exit(see function for description)        
   [M2] = imageNorm(M2,Xn,Yn);         clear Master_U;                                     
% Mean Schlieren images & fluctuating components used for analysis
   schMean = mean(M2,3);               flucM2 = M2 - schMean;         clear M2;
%% Image Matrix Dimensions & Zeroing terms 
%  Dimensions of image matrix
   rowNo = size(flucM2,1);            colNo = size(flucM2,2);         imgNo = size(flucM2,3);                                
%  Half Column location for zeroing during spatial FFT
   halfCol = round(colNo/2);     
%% Temporal Fourier Transform Computation
   %  Block Size 
%    [blockSize] = blkSizer(config); 
%    tempFFT = fft(flucM2,blockSize,3);     freqRes = 1/dt/blockSize;    schFreq = 0:freqRes:1/dt-1;
%% FOURIER DECOMPOSITION: Computing n - dimensional Fourier transform by using overlapping blocks
% Block size, overlap amount & number of Blocks
   blockSize = size(flucM2,3)/2;              nOvlp = 0.5*blockSize;      halfImgNo = round(blockSize/2); 
   nBlks = (size(flucM2,3)- nOvlp)/(blockSize - nOvlp);    
% Separated component matrices
   dwnStrmComp = zeros(size(flucM2));         upStrmComp = zeros(size(flucM2)); 
% Spatial FFT Loop
   for iBlk = 1:nBlks
       offset   = min((iBlk-1)*(blockSize - nOvlp)+blockSize,size(flucM2,3))-blockSize;
       nImgs    = (1:blockSize) + offset;           
% n-Dimensional FFT
       blockFFT = fftn(flucM2(:,:,nImgs));    fwdFFT = blockFFT;    bckFFT = blockFFT;    clear blockFFT;       
% Extracting Downstream & Upstream Components
       fwdFFT(:,1:halfCol+1,1:halfImgNo)       = 0;      fwdFFT(:,halfCol+1:end,halfImgNo+1:end) = 0;
       bckFFT(:,1:halfCol,halfImgNo+1:end)     = 0;      bckFFT(:,halfCol+1:end,1:halfImgNo)     = 0;
% Inverse FFT to recompute real part & phase of flow 
       dwnStrmComp(:,:,nImgs(1):nImgs(end)) = ifftn(fwdFFT);     clear fwdFFT;
       upStrmComp(:,:,nImgs(1):nImgs(end))  = ifftn(bckFFT);     clear bckFFT;
   end;    disp([newline '  <---------- Spatial Decomposition Complete ---------->']); 
%% .gif DISPLAY
% ------------------------------ Full Axial Domain Input Struct Definition -------------------------
   prompt = [newline,'Specify starting image location(1-2000):- '];    strtLoc = input(prompt);
   frameNos = 50;          frameEnd = strtLoc+frameNos;         imgRange = strtLoc:1:frameEnd-1;         
   convctnDisply.nt      = frameNos;      convctnDisply.frameDt = 0.5;          convctnDisply.cLI   = [0.05 0.5];    
   convctnDisply.Xn      = Xn;            convctnDisply.Yn    = Yn;             convctnDisply.Sz    = [20 40 745 845];          
   convctnDisply.lim_x   = limX;          convctnDisply.lim_y = [0 2];          convctnDisply.driveOut = driveOut;
   convctnDisply.TYP     = 'CONVC-SEP';   convctnDisply.xName = xName;          convctnDisply.yName = yName; 
   convctnDisply.Ftag    = 'inst';        convctnDisply.Ctag  = 'comp';         convctnDisply.NF    = NF;
   convctnDisply.nozzle  = nozzle;        convctnDisply.cMap  = schFaMap;       convctnDisply.cLC   = [-0.05 0.05];          
   convctnDisply.Ftitle  = ['$\bf NPR \thinspace',num2str(NPR),': Instantaneous \thinspace Schlieren \thinspace Image$'];
   convctnDisply.Dtitle  = ['$\bf NPR \thinspace',num2str(NPR),':Downstream \thinspace Components$'];
   convctnDisply.Utitle  = ['$\bf NPR \thinspace',num2str(NPR),':Upstream \thinspace Components$'];              
   convctnDisply.gifName = ['FourSepAxLen_',nozzle,'_','NPR ',num2str(NPR),'-imgStrt-',num2str(strtLoc)];

   imgMat.Fluc = flucM2(:,:,imgRange)+schMean;    imgMat.DwnStrm = dwnStrmComp(:,:,imgRange);    
   imgMat.UpStrm = upStrmComp(:,:,imgRange);      GF_GIF_DisplayWrite(imgMat,convctnDisply);
%% ------------------------------ Near Nozzle Region input Struct Parameters ------------------------
   convctnDisply.lim_y   = [0 1.2];       convctnDisply.lim_x = [0 3];      convctnDisply.Sz = [20 40 405 845]; 
   convctnDisply.Ftitle  = ['$\bf NPR \thinspace',num2str(NPR),':NF \thinspace Instantaneous$'];
   convctnDisply.Dtitle  = ['$\bf NPR \thinspace',num2str(NPR),':NF \thinspace$'];
   convctnDisply.Utitle  = ['$\bf NPR \thinspace',num2str(NPR),':NF \thinspace Upstream$']; 
   convctnDisply.cMap    = schFaMap;       convctnDisply.frameDt = 1;
   convctnDisply.gifName = ['FourSepNF_',nozzle,'_NPR ',num2str(NPR),'-imgStrt-',num2str(strtLoc)];    GF_GIF_DisplayWrite(imgMat,convctnDisply);
%% Save Upstream & Downstream matrices
   prompt = [newline,'Save Component Matrices(y/n)- '];    chc = input(prompt,"s");
   if strcmp(chc,'y')
      save([driveIn,'Upstream_Compnent'],'upStrmComp');    save([driveIn,'Downstream_Compnent'],'dwnStrmComp');
   end
%% Extracting axial & temporal intensity data from a single Radial location
% Mean Schlieren image 
   GF_FigurePlot(schMean,Xn,Yn);   crntFig = gcf;   figNo1 = crntFig.Number;   aX = gca;   set(aX,'FontSize',12);
   aX.TickLabelInterpreter = 'latex';      colorbar off;       xlabel(xName,'Interpreter','latex');
   ylabel(yName,'Interpreter','latex');    hold on;            xlim(limX);   ylim(limY);   set(gcf,'Position', [25 430 560 420])          
   caxis([0 0.5]);                         colormap gray;      set(gca, 'layer', 'top');   title('$\bf Average \thinspace Schlieren \thinspace Image$','Interpreter','latex');
% Instantaneous Schlieren image 
   instFig = schMean+flucM2(:,:,100);
   GF_FigurePlot(instFig,Xn,Yn);   crntFig = gcf;   figNo2 = crntFig.Number;   aX = gca;   set(aX,'FontSize',12);
   aX.TickLabelInterpreter = 'latex';      colorbar off;       xlabel(xName,'Interpreter','latex');
   ylabel(yName,'Interpreter','latex');    hold on;            xlim(limX);   ylim(limY);   set(gcf,'Position',[590 430 560 420])         
   caxis([0 0.6]);                         colormap gray;      set(gca, 'layer', 'top');   title('$\bf Instantaneous \thinspace Schlieren \thinspace Image$','Interpreter','latex')
% Limits for separated components & contour limits
   limX2 = [0 0.008];                      limY2 = [0 5];      cLim = [-0.05 0.05];
% Taking radial location as input 
   prompt = [newline,'->> Enter radial location to extract data(-4 to 4)- '];  radLoc = input(prompt);
   if radLoc < 0
      [locTn,~] = find(Yn > radLoc);       locTn = locTn(1) - 1;     
   else
      [locTn,~] = find(Yn < radLoc);       locTn = locTn(end) + 1;  
   end;   disp([newline,'>> -- Radial Location: - ',num2str(Yn(locTn))])
% Plotting that line on average & instantaneous figures
   yLoctn = zeros(size(Xn));               yLoctn(:) = Yn(locTn);
   figure(figNo1);                         plot(Xn,yLoctn,'-.w','LineWidth',1.2) ; 
   figure(figNo2);                         plot(Xn,yLoctn,'-.w','LineWidth',1.2) ; 
% Extracting time series data for that radial location from fluctuating component
   timeSeries = flucM2(locTn,:,:);         timeSeries = reshape(timeSeries,size(timeSeries,2),size(timeSeries,3));
%% -----------------------------------  Plotting base time histories ------------------------------------
% Data recording time 
   recTime = (dt:dt:size(flucM2,3)*dt);   GF_FigurePlot(timeSeries',Xn,recTime);      baseFig = gcf;   axis normal;
   aX = gca;            set(aX,'FontSize',12);   aX.TickLabelInterpreter = 'latex';   caxis(cLim);
   colorbar off;        ylabel('$Time(s)$','Interpreter','latex');   xlabel(xName,'Interpreter','latex');
   colormap gray;       set(baseFig,'Position',[15 85 540 770]);     ylim(limX2);     xlim(limY2)
   title(['$Time \thinspace Series \thinspace at \thinspace Y = ',num2str(round(Yn(locTn),1)),'$'],'Interpreter','latex');
% --------------------------------- Downstream Time Series ----------------------------------------- 
   dwnStrmSeries = dwnStrmComp(locTn,:,:); 
   dwnStrmSeries = reshape(dwnStrmSeries,size(dwnStrmSeries,2),size(dwnStrmSeries,3));
   dwnStrmImage  = real(squeeze(dwnStrmSeries));
% Plot time series
   GF_FigurePlot(dwnStrmImage',Xn,recTime);   dwnStrm = gcf;   axis normal;     
   aX = gca;            set(aX,'FontSize',12);   aX.TickLabelInterpreter = 'latex';   caxis( cLim);
   colorbar off;        ylabel('$Time(s)$','Interpreter','latex');   xlabel(xName,'Interpreter','latex');
   colormap gray;       set(dwnStrm,'Position',[205 85 540 770]);    ylim(limX2);     xlim(limY2);         
   title(['$Downstream \thinspace Series \thinspace at \thinspace Y = ',num2str(round(Yn(locTn),1)),'$'],'Interpreter','latex');
% ---------------------------------- Upstream Time Series ------------------------------------------ 
   upStrmSeries = upStrmComp(locTn,:,:); 
   upStrmSeries = reshape(upStrmSeries,size(upStrmSeries,2),size(upStrmSeries,3));
   upStrmImage  = real(squeeze(upStrmSeries));
% Plot time series
   GF_FigurePlot(upStrmImage',Xn,recTime);   upStrm = gcf;   axis normal;     
   aX = gca;            set(aX,'FontSize',12);   aX.TickLabelInterpreter = 'latex';   caxis(cLim)
   colorbar off;        ylabel('$Time(s)$','Interpreter','latex');   xlabel(xName,'Interpreter','latex');
   colormap gray;       set(upStrm,'Position',[305 85 540 770]);     ylim(limX2);     xlim(limY2);
   title(['$Upstream \thinspace Series \thinspace at \thinspace Y = ',num2str(round(Yn(locTn),1)),'$'],'Interpreter','latex');
%--------------------------------- Extract edges from the images -----------------------------------
% Extracting edges from downstream components image
   edgesDwnStrm = double(edge(dwnStrmImage,'canny',[0 0.2]));   [X1,recTime1] = meshgrid(Xn',recTime);      
   figure;pcolor(X1,recTime1,edgesDwnStrm');     dwnEdg = gcf;       axis normal; shading flat;       
   aX = gca;            set(aX,'FontSize',12);   aX.TickLabelInterpreter = 'latex';   caxis([0 1]);
   colorbar off;        ylabel('$Time(s)$','Interpreter','latex');   xlabel(xName,'Interpreter','latex');
   colormap(blckToRed); set(dwnEdg,'Position',[205 85 540 770]);     ylim(limX2);     xlim(limY2);
   title(['$Downstream \thinspace Edges \thinspace at \thinspace Y = ',num2str(round(Yn(locTn),1)),'$'],'Interpreter','latex');
% Extracting edges from upstream components image
   edgesUpStrm  = double(edge(upStrmImage,'canny')); 
   figure;pcolor(X1,recTime1,edgesUpStrm');         upEdg = gcf;      axis normal;    shading flat;        
   aX = gca;            set(aX,'FontSize',12);   aX.TickLabelInterpreter = 'latex';   caxis([0 1]);
   colorbar off;        ylabel('$Time(s)$','Interpreter','latex');   xlabel(xName,'Interpreter','latex');
   colormap(blckToRed); set(upEdg,'Position',[305 85 540 770]);      ylim(limX2);     xlim(limY2);
   title(['$Upstream \thinspace Edges \thinspace at \thinspace Y = ',num2str(round(Yn(locTn),1)),'$'],'Interpreter','latex');
%% INPUT REGIONS FROM UPSTREAM & DOWNSTREAM COMPONENTS TO COMPUTE VELOCITIES
%% Downstream Ranges
   figure(dwnEdg.Number);       disp([newline,'->> Downstream Ranges']);
% Enter graphical input by clicking a point on the edge image
   [x,y,~] = ginput(1);         xStrtD = x(1);    timeStrtD = y(1); 
% Fixed interval for spatial & temporal end points - (Rectangular:Xs = 47;Ts = 20)
   dx = Xn(2) - Xn(1);          xEndD = xStrtD + (dx*47);    timeEndD = timeStrtD + (dt*20);
   disp([newline,'>>- Starting axial location - ',num2str(xStrtD)]);   disp([newline,'>>- Starting time - ',num2str(timeStrtD)]);
   disp([newline,'>>- Ending axial location   - ',num2str(xEndD)]);    disp([newline,'>>- Ending time   - ',num2str(timeEndD)]);
   figure(dwnEdg.Number);       ylim([timeStrtD timeEndD]);         xlim([xStrtD xEndD]); 
% Defining start & end values in time and axial direction - Downstream
   [~,xStLocD] = find(X1(1,:)<xStrtD);            [~,xEnLocD] = find(X1(1,:)>xEndD);
   [tStrtD,~]  = find(recTime1(:,1)<timeStrtD);   [tEndD,~]   = find(timeEndD>recTime1(:,1));
% Cropping axial & temporal scales to redefine grid size  - Downstream 
   X2D = X1(1,xStLocD(end)+1:xEnLocD(1)-1);       recTime2D   = recTime(1,tStrtD(end)+1:tEndD(end));
   [X3D,recTime3D] = meshgrid(X2D,recTime2D);
% Cropping downstream image to match scales
   cropdDwnStrm = edgesDwnStrm(xStLocD(end)+1:xEnLocD(1)-1,tStrtD(end)+1:tEndD(end));
%% Upstream Ranges
   figure(upEdg.Number);        disp([newline,'->> Upstream Ranges']);
 % Enter graphical input by clicking a point on the edge image
   [x,y,~] = ginput(1);         xStrtU = x(1);    timeStrtU = y(1); 
% Fixed interval for spatial & temporal end points - (Rectangular:Xs = 47;Ts = 20)
   dx = Xn(2) - Xn(1);          xEndU = xStrtU + (dx*47);    timeEndU = timeStrtU + (dt*20);
   disp([newline,'>>- Starting axial location - ',num2str(xStrtU)]);   disp([newline,'>>- Starting time - ',num2str(timeStrtU)]);
   disp([newline,'>>- Ending axial location   - ',num2str(xEndU)]);    disp([newline,'>>- Ending time   - ',num2str(timeEndU)]);
   figure(upEdg.Number);        ylim([timeStrtU timeEndU]);         xlim([xStrtU xEndU]);  
% Defining start & end values in time and axial direction - Upstream 
   [~,xStLocU] = find(X1(1,:)<xStrtU);            [~,xEnLocU] = find(X1(1,:)>xEndU);
   [tStrtU,~]  = find(recTime1(:,1)<timeStrtU);   [tEndU,~]   = find(timeEndU>recTime1(:,1));
% Cropping axial & temporal scales to redefine grid size  - Upstream
   X2U = X1(1,xStLocU(end)+1:xEnLocU(1)-1);        recTime2U   = recTime(1,tStrtU(end)+1:tEndU(end));
   [X3U,recTime3U] = meshgrid(X2U,recTime2U);
% Cropping upstream image to match scales
   cropdUpStrm = edgesUpStrm(xStLocU(end)+1:xEnLocU(1)-1,tStrtU(end)+1:tEndU(end));
%% ---------------------------------- Downstream Components ----------------------------------------
   disp([newline '>-----------------------------------------------------<']);
   patrnIdentify(recTime2D,recTime3D,X2D,X3D,cropdDwnStrm,'Downstream',blckToRed,schJet,nozHt,Uj,xName)
%% ----------------------------------- Upstream Components ----------------------------------------
   disp([newline '>-----------------------------------------------------<']);
   patrnIdentify(recTime2U,recTime3U,X2U,X3U,cropdUpStrm,'Upstream',blckToRed,schJet,nozHt,spdSound,xName)
   %% Velocity computation using cross correlation from axial locations
axRes = Xn(2) - Xn(1);    disp([newline,'->> Axial Resolution(axRes): ',num2str(axRes),xName(4)])
prompt = [newline,'->> Enter first  X - location(0-6) - '];                    inpX1   = input(prompt);  
prompt = [newline,'->> Enter axial steps to 2nd X - location(n*axRes) - '];    axSteps = input(prompt);
% Finding the exact column from the axial location
[xStrt,~,~,~] = precisLoc(inpX1,[],Xn,[]);         xEnd = xStrt + axSteps; 
disp([newline,'->> Second axial location: ',num2str(Xn(xEnd))]);              axDist = (Xn(xEnd)-Xn(xStrt))*nozHt;
disp([newline,'->> Axial distance(m) between points: ',num2str(axDist)]);      
% Function to compute velocity
crosCorelVelFind(xStrt,xEnd,dwnStrmImage,dt,axDist,Uj)

%% Functions:
%%      Function 1: DRIVE CHECK FOR VALID INPUT DRIVE
function driveCheck(driveIn)
    if isfolder(driveIn) == 0
       disp([newline,'->> Check nozzle & condition']); return;
    end
end
%%      Function 2: BACKGROUND SUBTRACTION
function [M2,schMeanOne,schMeanTwo] = bckgrndSub(bckgrndPath,Master_U,config)
    if isempty(bckgrndPath) ~= 1 
        disp([newline '->> Subtracting Background' newline]);
        Bckgrnd = load([bckgrndPath 'Bckgrnd_DAT']);    
        Bckgrnd = Bckgrnd.Master_U;                         schMeanOne = mean(Master_U,3);  
       if strcmp(config,'TS1') ~= 1
          M2 = Master_U./mean(Bckgrnd,3);                  schMeanTwo = mean(M2,3);
       end
    else
       schMeanTwo = mean(Master_U,3);                      schMeanOne = schMeanTwo; 
    end
end
%%      Function 3: COMPUTE BLOCKSIZE BASED ON ACQUISITION RATE
function [blkSize] = blkSizer(config)                                                              %#ok
    nozlConfig = {'C'; 'S'; 'S2';'SS2'; 'TR';'TR2';'TS';'TS1';'TS2';'TRV0'};
    frameRate  = [45E3;41E3;204E3;204E3;41E3;204E3;41E3; 112E3;204E3; 45E3];
    fftBlks    = [900; 820; 2048; 2048; 820; 2048; 820; 2240; 2048;  900];
    schTestConfig = table(nozlConfig,frameRate,fftBlks);
    inDx = strcmp(schTestConfig.nozlConfig,config);
    blkSize = schTestConfig.fftBlks(inDx);
end
%%      Function 4: IMAGE VALUE NORMALIZATION
function [Master_U] = imageNorm(Master_U,X,Y)
    % Bright peaks on the nozzle wall & the region outside the flow due to 
    % background subtraction disproportionately skew values making the 
    % rest of the image very dim thus lowering the dynamic range available
    % for flow recognition. The function therefore uses only the values 
    % after the nozzle exit in the range Y(-3.5 3,5) & X(0 5) to normalize 
    % the intensity for each image
    xZero = find(X==0);      xZero = xZero+2;      [~,E] = find(X<5);      xEnd = E(end)+1;       
    [S,~] = find(Y>-3.5);    yStrt = S(1);         [E,~] = find(Y<3.5);    yEnd = E(end);
    for ctr = 1:size(Master_U,3)
        Master_U(:,:,ctr) = Master_U(:,:,ctr)/max(max(Master_U(yStrt:yEnd,xZero:xEnd,ctr)));
    end;    disp([newline '   >---------- Normalization Complete ----------<']);
end
%%      Function 5: IDENTIFY PATTERNS IN THE FOURIER SEPERATED COMPONENTS
function patrnIdentify(recTime2,recTime3,X2,X3,schComponent,tag,blckToRed,schJet,nozHt,limitVel,xName)
% Plotting cropped image
    figure;   pcolor(X3,recTime3,schComponent');  a = gcf;            shading flat;    figNo1 = a.Number;   
    aX = gca;            set(aX,'FontSize',12);   aX.TickLabelInterpreter = 'latex';   caxis([0 1])
    colorbar off;        ylabel('$Time(s)$','Interpreter','latex');   xlabel(xName,'Interpreter','latex');
    colormap(blckToRed); set(a,'Position',[305 85 540 770]);                    
    title(['$Cropped \thinspace ',tag,' \thinspace Series$'],'Interpreter','latex');
% Hough transform computation to identify linear patterns 
    [H,theTa,disTnc] = hough(schComponent,'RhoResolution',0.5,'Theta',-90:0.3:89);  %default: -90:0.2:89
% Identifying peak values to extract angle(theTa) & distance(disTnc)
    peaks   = houghpeaks(H,20);  %default: (H,20)
% Plotting transform & distances
    GF_FigurePlot(H,theTa,disTnc);    a1 = gcf;    axis normal;    colormap(schJet);    caxis([0 30]); 
    set(gcf,'Position',[190 340 650 379]);         c = colorbar;   c.TickLabelInterpreter = 'latex';
    set(gca,'FontSize',12);          aX = gca;     aX.TickLabelInterpreter = 'latex';
    xlabel('$Theta(\theta)$','Interpreter','latex');    ylabel('$Distance(\rho)$','Interpreter','latex');
    title(['$Hough \thinspace Transform:',tag,'$'],'Interpreter','latex');
    title(c,'$H$','Interpreter','latex');     hold on;     set(a1,'Position',[65 270 700 510]);
% Plotting the peaks 
    plot(theTa(peaks(:,2)),disTnc(peaks(:,1)),'s','color','white');
% Extracting lines from the peak points
    lines = houghlines(schComponent,theTa,disTnc,peaks);
% Showing the match between pattern & Hough Lines
    figure(figNo1);   hold on;   locMat = zeros(numel(lines),4);   velVals = zeros(numel(lines),1); 
    for ctr = 1:numel(lines)
        x1 = lines(ctr).point1(1);     t1 = recTime2(x1);
        y1 = lines(ctr).point1(2);     p1 = X2(y1);            
        x2 = lines(ctr).point2(1);     t2 = recTime2(x2);
        y2 = lines(ctr).point2(2);     p2 = X2(y2); 
        locMat(ctr,:) = [p1,p2,t1,t2];
        p1 = p1*nozHt;                 p2 = p2*nozHt;
        lineSlope = (p2-p1)/(t2-t1);   velVals(ctr,1) = lineSlope;
    end
%  Removing infinite values
    [inDx1,~] = find(velVals==Inf);            [inDx2,~] = find(velVals==-Inf);
    velVals(inDx1) = [];                       velVals(inDx2) = [];
    locMat(inDx1,:) = [];                      locMat(inDx2,:) = []; 
    if strcmp(tag,'Downstream')
       [inDx3,~] = find(velVals<0);            velVals(inDx3) = [];       locMat(inDx3,:) = [];
       [inDx4,~] = find(velVals>limitVel);     velVals(inDx4) = [];       locMat(inDx4,:) = [];
       [inDx5,~] = find(velVals<0.3*limitVel); velVals(inDx5) = [];       locMat(inDx5,:) = [];  
    else
       [inDx3,~] = find(velVals>0);            velVals(inDx3) = [];       locMat(inDx3,:) = [];
       [inDx4,~] = find(velVals<limitVel*-1);  velVals(inDx4) = [];       locMat(inDx4,:) = [];
       velVals = abs(velVals);
    end
%  Plotting the remaining points on time history contour
    for ctr = 1:size(locMat,1)
        p1 = locMat(ctr,1);    p2 = locMat(ctr,2);    t1 = locMat(ctr,3);    t2 = locMat(ctr,4);
        plot([p1 p2],[t1 t2],'-.','Color','w','LineWidth', 2);
        disp([newline,'->> Velocity line(',num2str(ctr),'-',tag,'): ',num2str(velVals(ctr))]);
    end
    hold off;    disp([newline '>-----------------------------------------------------<']);
    disp([newline,'->> Mean Velocity(',tag,'): ',num2str(mean(velVals))]);
    disp([newline,'->> Max  Velocity(',tag,'): ',num2str(max(velVals))]);
    disp([newline,'->> Med  Velocity(',tag,'): ',num2str(median(velVals))]);
    disp([newline,'->> Mean Velocity Ratio(',tag,'): ',num2str(mean(velVals)/limitVel)]);
    disp([newline,'->> Peak Velocity Ratio(',tag,'): ',num2str(max(velVals)/limitVel)]);
    disp([newline,'->> Med  Velocity Ratio(',tag,'): ',num2str(median(velVals)/limitVel)]);
    disp([newline '>-----------------------------------------------------<']);
end
%%      Function 6: VELOCITY COMPUTING USING CROSS CORRELATION
function crosCorelVelFind(xCol1,xCol2,timeSeries,dt,axDist,Uj)
    sigOne = timeSeries(xCol1,:);    sigTwo = timeSeries(xCol2,:);
    [corel,Lag] = xcorr(sigOne,sigTwo,'normalized');    Lag = Lag * dt;    figure;    plot(Lag*dt,corel,'k','LineWidth',1.2);
    grid on;   box off;   xlabel('$\bf \Delta t(s)$','Interpreter','latex');   ylabel('$\bf \hat{C}$','Interpreter','latex');
    ax = gca;  ax.TickLabelInterpreter = 'latex';   ax.FontSize = 12;        
    [~,b] = max(corel);    travlTime = dt*Lag(b);       convecVelocity = axDist/ abs(travlTime);
    disp([newline,'->> Convective Velocity(m/s): ',num2str(convecVelocity)]);
    disp(['->> Convective Velocity(Uc/Uj): ',num2str(convecVelocity/Uj)]);
end
%%      Function 7: FIND PRECISE LOCATION
function [xCol,xVal,yRow,yVal] = precisLoc(inpX,inpY,Xn,Yn)
    xCol = 0;    xVal = 0;    yRow = 0;    yVal = 0;
%  Finding the precise X location match for input value
    if size(Xn) ~=0
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
    end
%  Finding the precise X location match for input value
    if size(Yn) ~= 0
       [y,~] = find(Yn == inpY);
       if isempty(y) && size(Yn) ~= 0
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
end
%% DEPRICATED CODE
% Background Subtraction
% if isempty(OutputStruct.bckgrnd_root) ~= 1 
%    disp([newline '->> Subtracting Background' newline]);
%    Bckgrnd = load([OutputStruct.bckgrnd_root 'Bckgrnd_DAT']);    
%    Bckgrnd = Bckgrnd.Master_U;                         schMeanOne = mean(Master_U,3);  
%    if strcmp(config,'TS1') ~= 1
%       M2 = Master_U./mean(Bckgrnd,3);                  schMeanTwo = mean(M2,3);
%    end
% else
%    Bckgrnd = [];  schMeanTwo = mean(Master_U,3);       schMeanOne = schMeanTwo; 
% end
%--------------------------------------------------------------------------------
% Spatial FFT Code v1
%    imgFourTran = fftn(flucMasterOne(:,:,1:1000));             % imgFourTran = fftshift(imgFourTran,3);
% %  Zeroing diagonal segments for downstream & upstream components 
%    fwdFourTrans = imgFourTran;                  bckFourTrans = imgFourTran;
% %  Downstream components
%    fwdFourTrans(:,1:halfCol,1:halfImg)         = 0;   % Top-Left     : Q2
%    fwdFourTrans(:,halfCol+1:end,halfImg+1:end) = 0;   % Bottom-Right : Q4
% %  Upstream components
%    bckFourTrans(:,1:halfCol,halfImg+1:end)     = 0;   % Bottom-Left  : Q1
%    bckFourTrans(:,halfCol+1:end,1:halfImg)     = 0;   % Top-Right    : Q3 
% %  Inverting fft to get image back respective components
%    imgDwnStrmComp = ifftn(fwdFourTrans);       imgUpStrmComp = ifftn(bckFourTrans); 
%--------------------------- Compute Slope from Seperated image sets--------------------------
%    velSet = zeros(4,1);    
%    for ctr = 1:4
%        prompt = [newline,'->> Y1 - '];       Y1 = (input(prompt)*nozHt);
%        prompt = [newline,'->> X1 - '];       X1 = input(prompt);
%        prompt = [newline,'->> Y2 - '];       Y2 = (input(prompt)*nozHt);
%        prompt = [newline,'->> X2 - '];       X2 = input(prompt);
% % Slope(velocity)
%        slopeVel = (Y2 - Y1)/(X2-X1);         disp([newline,'->> Velocity - ',num2str(slopeVel),'m/s']);
%        velSet(ctr) = slopeVel;
%    end
%    avgVel = mean(velSet,1);                disp([newline,'->> Velocity - ',num2str(avgVel),'m/s']);
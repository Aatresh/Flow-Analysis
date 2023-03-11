%% POD ANALYSIS FOR SCHLIEREN OR SHADOWGRAPH IMAGES
%  COMPUTES & SAVES POD MODE ENERGIES & TEMPORAL COMPONENTS FOR DIFFERENT JET
%  CONFIGURATIONS
fclose all;  clc;  clearvars;   set(0,'defaultfigurecolor',[1 1 1]);  code = 'Sch_PODCalc';
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';
cd([root1 'Jet_Analysis\Global_Functions\']); addpath([root1 'Jet_Analysis\POD_Codes\']);      
tests = {'NPR_2p5_TR_1p0', 'NPR_2p6_TR_1p0', 'NPR_2p9_TR_1p0','NPR_3p0_TR_1p0',...
         'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_4p5_TR_1p0','NPR_5p0_TR_1p0', 'NPR_6p0_TR_1p0'};

condition = [tests(1)];   config = 'C';   nozzle = 'Minor';   NF = 'D';    write = 'n';  modeBackup = 'n';

%  Drive slection based on configuration & frame rate
[OutputStruct] = GF_DriveSelect(config,nozzle,code);
%  Figure tags
if strcmp(NF,'D')
   tag = '(D)';   
else
   tag = '(h)';
end
%% LOADING DATA, BACKGROUND SUBTRACTION
for n = 1:length(condition)
    [Mj,Uj,NPR,NTR] = GF_Velocity(condition{n});
    drive_in = [OutputStruct.in_root condition{n}(9:14) '\' condition{n}(1:7)  '\'];
    driveOut = [OutputStruct.out_root '\' condition{n}(9:14) '\' condition{n}(1:7) '\'];
%  Loading video .mat file
    load([drive_in condition{n}(1:7) '_DAT']);   load([drive_in 'X']);  load([drive_in 'Y']);
    if isempty(OutputStruct.bckgrnd_root) ~= 1 
       disp([newline '->> Subtracting Background' newline]);
       Bckgrnd = load([OutputStruct.bckgrnd_root 'Bckgrnd_DAT']);           bckLim = [0.1 1.8];
       Bckgrnd = Bckgrnd.Master_U;                                          avgTwo = mean(Master_U,3);  
       M2 = Master_U./mean(Bckgrnd,3);      avgOne = mean(M2,3);
    else
       bckLim = [0 500];   M2 = Master_U;   avgOne = mean(Master_U,3);      avgTwo = avgOne; 
    end; disp([newline,'    <----   Done!   ---->']);
%% AXIS DEFINITIONS & POD COMPUTATION
    [Xn,Yn,limX,limY,xName,yName,lenScales,figSize] = GF_AxisDefnSch(config,nozzle,NF,X,Y);
    Xn = Xn';   Yn = Yn';
%  Average Schlieren
    schMean1 = mean(Master_U,3);    schMean2 = mean(M2,3);     imgNos = size(Master_U,3);
    rowNos   = size(Master_U,1);    colNos   = size(Master_U,2); 
%  Fluctuating component
    imgFluc = M2 - schMean2;        %imgFluc = Master_U - schMean1;                       
    imgTemp = permute(reshape(imgFluc,1,rowNos*colNos,imgNos),[2,3,1]);
%  Singular Value Decomposition
    disp('--> SVD Computation');  
    covarRho = imgTemp.'*imgTemp;    [eigenVec,eigenVal] = svd(covarRho); 
%  Mode Energies & POD Mode Extraction
    lambdaRho = diag(eigenVal);     disp('--> Computing POD Modes')
    phiBase = imgTemp*eigenVec;     clear covarRho Master_U imgFluc;
%  Normalizing(2-Norm) Spatial modes
    betaRho = zeros(size(phiBase));  wB = waitbar(0,'Computing 2-Norm for Spatial Modes');
    for ctr = 1:imgNos
        betaRho(:,ctr) = phiBase(:,ctr)/norm(phiBase(:,ctr),2);  waitbar(ctr/imgNos,wB);         
    end; close(wB);
%  Temporal co-efficient definition(can use either formulation
%     tempRho = betaRho'*imgTemp;         % 1 - Bhupa's method; 2 - Mitchell
    tempRho = eigenVec*norm(imgTemp*eigenVec);   clear imgTemp Bckgrnd;
% Image Normalization(don't use for SVD, causes low fluctiation modes to precipitate up)
    [M2] = imageNorm(M2,Xn,Yn);
%% IMAGE RECONSTRUCTION
%  Modes used for re-construction & No. of images to re-construct
    recModes = 5;     recImgs = imgNos; 
%  Re-constructed Image matrix
    imgRec = zeros(rowNos*colNos,recImgs);   wB = waitbar(0,['Re-Constructing Images using ',num2str(recModes),' Modes']);
    for ctr = 1:recImgs 
       imgRec(:,ctr) = betaRho(:,1:recModes)*tempRho(1:recModes,ctr); waitbar(ctr/imgNos,wB);
    end; close(wB);
% Re-arranging matrix order & adding mean
    newMaster = reshape(imgRec,[rowNos colNos recImgs]);
    newMaster = newMaster + schMean2;   clear imgRec;       
    gifImgs   = permute(newMaster,[3 1 2]);
% Creating input Struct file to tranfer input for .gif & image save functions
    SPOD_INP.nozzle = nozzle;   SPOD_INP.lim_x = limX;             SPOD_INP.lim_y = limY;
    SPOD_INP.TYP = 'POD';       SPOD_INP.Sz = [490 440 560 420];   SPOD_INP.driveOut = driveOut;
    SPOD_INP.Xn = Xn;           SPOD_INP.Yn = Yn;                  SPOD_INP.NF = NF;
    SPOD_INP.gifName = [condition{1} tag 'PODRec'];                SPOD_INP.config = config;
    SPOD_INP.bckLim = [0 2];
    % SPOD_INP.Sz = [676 88 1116 886];
%  Generating & saving  reconstruction .gif
    GF_GIF_DisplayWrite(gifImgs,SPOD_INP);
%% RE-NORMALIZING SPATIAL MODES TO RANGE BETWEEN -1 TO 1
    phiRho = zeros(size(phiBase));  
    wB = waitbar(0,'Re-Normalizing Spatial Modes');
    for ctr = 1:imgNos
  % Axial Modes - Finding positive & negative values
        [posVal,~] = find(phiBase(:,ctr)>0);   [negVal,~] = find(phiBase(:,ctr)<0);
  % Normalizing with peak values
        phiRho(posVal,ctr) = phiBase(posVal,ctr)/max(phiBase(posVal,ctr));  
        phiRho(negVal,ctr) = phiBase(negVal,ctr)/abs(min(phiBase(negVal,ctr)));
        clear posVal negVal;  waitbar(ctr/imgNos,wB);
    end; close(wB);
    phiRho = reshape(phiRho,[rowNos colNos imgNos]);         clear phiBase;
%% TEST PLOTTING MODES
    for ctr = 1:10
        figure; set(gca,'FontSize',13,'FontName','times');  axis equal;  
        pcolor(Xn,Yn,phiRho(:,:,ctr)); axis equal; shading interp; caxis([-0.5 0.5]); colormap bluewhitered;
        title(['$NPR:',num2str(NPR),'- Mode - ',num2str(ctr),'$'],'Interpreter','latex'); xlim(limX); ylim(limY);
        xlabel(xName,'Interpreter','latex'); ylabel(yName,'Interpreter','latex');  
    end
%% SAVING TO DRIVE
%  prompt = [newline '--> Save Data to Drive(y/n): ']; write = input(prompt,'s'); 
    if strcmp(write,'y')
       disp('<< - Saving Data to Drive - >>');  
       phiRho = phiRho(:,:,1:100);   
       save([driveOut 'modeEnergy-Rho'],'lambdaRho'); 
       save([driveOut 'PODModes-Rho'],'phiRho');      save([driveOut 'X'],'X');
       save([driveOut 'tempCoeff-Rho'],'tempRho');    save([driveOut 'Y'],'Y');
    end
% Saving all POD modes to Schlieren backup SSD (Data_Three)
    if strcmp(modeBackup,'y')
       disp('<< - Taking Mode Backup - >>'); 
       driveOut = ['Y:\Schlieren Data\POD_Reconstruction\' OutputStruct.out_root(74:end) condition{n}(9:14) '\' condition{n}(1:7) '\'];
       save([driveOut 'modeEnergy-Rho'],'lambdaRho');                             % Mode Energies
       save([driveOut 'PODModes-Rho'],'phiRho');      save([driveOut 'X'],'X');   % All POD modes  & X
       save([driveOut 'tempCoeff-Rho'],'tempRho');    save([driveOut 'Y'],'Y');   % Temporal modes & Y
       save([driveOut 'schMean'], 'schMean');         disp('< - Done ->');        % Average Schlieren 
   end
end
%% ------------------------------------ Functions ---------------------------------------------------------------%%
%%      Function 1: NORMALIZING IMAGE WITH HIGHEST INTENSITY VALUE
function [Master_U] = imageNorm(Master_U,X,Y)
    % Bright peaks on the nozzle wall & the region outside the flow due to 
    % background subtraction disproportionately skew values making the 
    % rest of the image very dim thus lowering the dynamic range available
    % for flow recognition. The function therefore uses only the values 
    % after the nozzle exit in the range Y(-3.5 3,5) & X(0 5) to normalize 
    % the intensity for each image
    xZero = find(X==0);      xZero = xZero+1;      [~,E] = find(X<5);      xEnd = E(end)+1;       
    [S,~] = find(Y>-3.5);    yStrt = S(1);         [E,~] = find(Y<3.5);    yEnd = E(end);
    for ctr = 1:size(Master_U,3)
        Master_U(:,:,ctr) = Master_U(:,:,ctr)/max(max(Master_U(yStrt:yEnd,xZero:xEnd,ctr)));
    end;    disp([newline '   >---------- Normalization Complete ----------<']);
end
%  SHUTTING DOWN PARALLEL POOL
% p = gcp;        delete(p);                          

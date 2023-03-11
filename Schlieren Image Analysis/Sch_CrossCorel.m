%% Code computes various types of autocorrelation & cross correlations
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
% Type of Analysis:
% 1. Auto Correlation at a single point
% 2. Cross correlation between Two single points
% 3. Cross correlation between Single point & Horizontal point array
% 4. Cross correlation between Single point & Vertical point array
% 5. Cross correlation between Two Linear point arrays
%% Main Code 
tic;   fclose all;   clc;   clearvars;   set(0,'defaultfigurecolor',[1 1 1]);   code = 'SchFourier';    
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';
driveCode = [root1 'Jet_Analysis\Global_Functions\']; cd(driveCode); addpath([root1 'Jet_Analysis\Schlieren & SPOD_Codes\']); 
colr = [0.9 0 0;0.9 0.6 0.1;0.8 0.5 0.8;0 0 0.9;0.6 0.6 0.6;0.9 0.2 0.5;0.2 0.2 0.2;0.5 0.9 0.6;0.2 0.6 0.9;];
load('custom_map3.mat');      load('schl_sc.mat');
tests = {'NPR_2p5_TR_1p0', 'NPR_2p6_TR_1p0', 'NPR_2p9_TR_1p0', 'NPR_3p0_TR_1p0',...
         'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_4p5_TR_1p0', 'NPR_5p0_TR_1p0', 'NPR_6p0_TR_1p0'};

condition = tests(1);    config = 'C';    nozzle = 'Major';    NF = 'D';   

%% SELECTING DRIVE BASED ON FRAME RATE 
[OutputStruct] = GF_DriveSelect(config,nozzle,code);   nozzle = OutputStruct.nozzle;   dt = OutputStruct.dt;
%   JET VELOCITY
[Mj,Uj,NPR,NTR] = GF_Velocity(condition{1});           disp([newline '--> Acquisition Rate - ',num2str(1/dt),' Hz']);
%  LOADING VIDEO MATRIX
driveIn  = [OutputStruct.in_root condition{1}(9:14) '\' condition{1}(1:7)  '\'];
driveOut = [OutputStruct.out_root condition{1}(9:14) '\' condition{1}(1:7) '\'];
driveCheck(driveIn);     load([driveIn condition{1}(1:7) '_DAT']);    load([driveIn 'X']);   load([driveIn 'Y']);

%% LOADING BACKGROUND FILE
if isempty(OutputStruct.bckgrnd_root) ~= 1 
   disp([newline '->> Subtracting Background' newline]);
   Bckgrnd = load([OutputStruct.bckgrnd_root 'Bckgrnd_DAT']);           cLim1 = [0 1.5];
   Bckgrnd = Bckgrnd.Master_U;                                          avgTwo = mean(Master_U,3);  
   M2 = Master_U./mean(Bckgrnd,3);  avgOne = mean(M2,3);
else
   Bckgrnd = [];  cLim1 = [0 500];  avgOne = mean(Master_U,3);          avgTwo = avgOne; 
end; disp([newline,'    <----   Done!   ---->']);  load schJet.mat;     
%  AXIS NORMALIZATION & LABELLING
[Xn,Yn,limX,limY,xName,yName,lenScales,figSize] = GF_AxisDefnSch(config,nozzle,NF,X,Y);   Xn = Xn';   Yn = Yn';

%% Main input strcut for function analysis
inputParams.imgMatOne = M2;        inputParams.imgMatTwo = Master_U;       inputParams.Xn   = Xn;       
inputParams.avgOne    = avgOne;    inputParams.avgTwo    = avgTwo;         inputParams.Yn   = Yn;       
inputParams.cLim1     = cLim1;     inputParams.limX      = limX;           inputParams.limY = limY; 
inputParams.xName     = xName;     inputParams.yName     = yName;          inputParams.dt   = dt;
inputParams.colr      = colr;      inputParams.flucMaster = M2 - avgOne;

%% Time axis for lags 
tTotal = size(M2,3)*dt;            tPos = 0:dt:tTotal;                     tNeg = fliplr(tPos)*-1;      
tLags  = cat(2,tNeg(2:end-1),tPos(1:end-1));        clear tPos tNeg;
%  Finding X = 0 & adjusting axis
[~,xStrt] = find(Xn==0);           xPos = Xn(xStrt:end);                   xNeg = fliplr(xPos)*-1; 
xLags = cat(2,xNeg(1:end-1),xPos(1:end));           clear xPos xNeg;
inputParams.tLags = tLags;         inputParams.xLags = xLags;              inputParams.xStrt = xStrt;
inputParams.xNew  = Xn(xStrt:end);

%% AVERAGE IMAGE
figure(1);   pcolor(Xn,Yn,avgOne);      shading interp;      set(gca,'FontSize',13);          ax = gca;  
ax.TickLabelInterpreter = 'latex';      axis equal;          xlim(limX);     ylim(limY);      colormap gray;  
xlabel(xName,'Interpreter','latex');    hold on;             set(gca, 'layer', 'top');        colorbar off;             
ylabel(yName,'Interpreter','latex');    caxis(cLim1);        set(gcf,'Position',figSize);
%% Display choice of analysis & trigger function selector
 compType = 3;                          compTypDisplay(compType);
%  Function Selector
if compType == 1
   singlPointAutcor(inputParams);  
elseif compType == 2
   twoPointCrosCorr(inputParams);
elseif compType == 3
   pointHLineCor(inputParams);
elseif compType == 4
   pointVLineCor(inputParams);
elseif compType == 5
   twinLineCor(inputParams)
end
%% IN-BUILT FUNCTIONS

%%      Function 1: DRIVE CHECK FOR VALID INPUT DRIVE
function driveCheck(driveIn)
    if isfolder(driveIn) == 0
       disp([newline,'->> Check nozzle & condition']); return;
    end
end
%%      Function 2: DISPLAY TYPE OF ANLYSIS CHOSEN
function compTypDisplay(compType)
    compList    = [1;2;3;4];
    compDescrip = {'Autocorrelation at Single point'; 'Cross correlation between Two single points';...
                   'Cross correlation between Single point & Horizontal array';...
                   'Cross correlation between Single point & Vertical array';...
                   'Cross correlation between Two Linear arrays'};
    inDx = compType == compList;
    disp([newline,'  ->> ',cell2mat(compDescrip(inDx)),' <<- ',newline]);
end
%%      Function 3: SINGLE POINT AUTO CORRELATION
function singlPointAutcor(inputParams)
    chc = 'y';   colrCtr = 1;   colr = inputParams.colr;   Xn = inputParams.Xn;   Yn = inputParams.Yn;   tLags = inputParams.tLags;
    while strcmp(chc,'y')
       prompt = [newline,'->> Enter Point X-Location - '];  inpX = input(prompt);
       prompt = [newline,'->> Enter Point Y-Location - '];  inpY = input(prompt);
%  Function to find exact location
       [xCol,~,yRow,~] = precisLoc(inpX,inpY,inputParams.Xn,inputParams.Yn);
%  Plotting marker lines & point locations
       pointMrk(1,xCol,yRow,Xn,Yn,'r',colrCtr);
%  Fluctuating component & Input Intensity Signal
       inpSignl = inputParams.flucMaster(yRow,xCol,:);       inpSignl = reshape(inpSignl,1,size(inputParams.flucMaster,3),[]);
%  Autocorrelation calculation, normalizing & plotting 
      [aCor,~] = xcorr(inpSignl,'coeff');      figure(2);      
      plot(tLags,aCor,'Color',colr(colrCtr,:),'LineWidth',1.2,'DisplayName',['$Location:(',num2str(round(Xn(xCol),2)),',',num2str(round(Yn(yRow),2)),')$']);
      xlabel('$\Delta \tau$','Interpreter','latex');        set(gca,'FontSize',13);    box off;       grid on;         hold on;
      ylabel('$R_{x_{1}x_{1}}$','Interpreter','latex');     ax = gca;    xlim([-35E-4 35E-4]);        c = legend;      c.EdgeColor = [1 1 1]; 
      ax.TickLabelInterpreter = 'latex';       set(gcf,'Position',[740 70 750 370]);   colrCtr = colrCtr+1;            c.Interpreter = 'latex';         
%  Conditon check to plot another point
      prompt = [newline,'->> Plot another point(y/n) - '];  chc = input(prompt,'s');
    end
end
%%      Function 4: POINT - POINT(Two Point)CROSS CORRELATION
function twoPointCrosCorr(inputParams)
    chc = 'y';        colrCtr = 1;      colr = inputParams.colr;
    Xn = inputParams.Xn;                Yn = inputParams.Yn;        tLags = inputParams.tLags;
    while strcmp(chc,'y')
%  Location of 1st point
       prompt = [newline,'->> Enter 1st Point X-Location - '];       inpX1 = input(prompt);
       prompt = [newline,'->> Enter 1st Point Y-Location - '];       inpY1 = input(prompt);
%  Loction of 2nd point
       prompt = [newline,'->> Enter 2nd Point X-Location - '];       inpX2 = input(prompt);
       prompt = [newline,'->> Enter 2nd Point Y-Location - '];       inpY2 = input(prompt);
%  Exact location of 1st point
       [xCol1,~,yRow1,~] = precisLoc(inpX1,inpY1,inputParams.Xn,inputParams.Yn);
%  Exact location of 2nd point
       [xCol2,~,yRow2,~] = precisLoc(inpX2,inpY2,inputParams.Xn,inputParams.Yn);
%  Plotting marker lines & point locations
       pointMrk(1,xCol1,yRow1,Xn,Yn,'r',1);                  pointMrk(1,xCol2,yRow2,Xn,Yn,'c',2);
%  Obtaining signals for both points
       inpSignl1 = inputParams.flucMaster(yRow1,xCol1,:);    inpSignl1 = reshape(inpSignl1,1,size(inputParams.flucMaster,3),[]);
       inpSignl2 = inputParams.flucMaster(yRow2,xCol2,:);    inpSignl2 = reshape(inpSignl2,1,size(inputParams.flucMaster,3),[]);
%  Computing cross correlation
       [crosCor,~] = xcorr(inpSignl1,inpSignl2,'coeff');     figure(2);      
       plot(tLags,crosCor,'Color',colr(colrCtr,:),'LineWidth',1.2,'DisplayName',['Set-',num2str(colrCtr)]);
       xlabel('$\Delta \tau$','Interpreter','latex');        set(gca,'FontSize',13);    box off;       grid on;         hold on;
       ylabel('$C_{x_{1}x_{2}}$','Interpreter','latex');     ax = gca;    xlim([-35E-4 35E-4]);        c = legend;      c.EdgeColor = [1 1 1]; 
       ax.TickLabelInterpreter = 'latex';       set(gcf,'Position',[740 70 750 370]);   colrCtr = colrCtr+1;            c.Interpreter = 'latex'; 
       title('$Two \thinspace Point \thinspace Cross-Correlation$','Interpreter','latex')
       prompt = [newline,'->> Plot another set(y/n) - '];  chc = input(prompt,'s');
    end
end
%%      Function 5: POINT - HORIZONTAL LINE CROSS CORRELATION
function pointHLineCor(inputParams)
    chc = 'y';                  xNew = inputParams.xNew;           ptCtr = 1;
    Xn = inputParams.Xn;        Yn = inputParams.Yn;               tLags = inputParams.tLags;
    while strcmp(chc,'y')
%  Location of 1st point
       prompt = [newline,'->> Enter Point X-Location(val>0) - '];    inpX = input(prompt);
       prompt = [newline,'->> Enter Point Y-Location - '];           inpY = input(prompt);
       prompt = [newline,'->> Enter Line  Y-Location - '];           inpLy = input(prompt);
%  Exact location of 1st point
       [xCol,~,yRow,~] = precisLoc(inpX,inpY,inputParams.Xn,inputParams.Yn);
%  Exact location of line
       [~,~,lRow,~]    = precisLoc(0,inpLy,inputParams.Xn,inputParams.Yn);
%  Plotting marker lines & point locations
       pointMrk(1,xCol,yRow,Xn,Yn,'r',ptCtr);
%  Plotting line
       yLoctn = zeros(size(Xn));          yLoctn(:) = Yn(lRow);
       figure(1);                         plot(Xn,yLoctn,'-.','Color',inputParams.colr(ptCtr,:),'LineWidth',1.2) ; 
%  Text indicating line
       text(3,Yn(lRow)+0.2,['$L:',num2str(ptCtr),'$'],'Color','k','FontSize',15,'Interpreter','latex');
%  Obtaining signals for point & line-plane(line though time)
       ptSignl  = inputParams.flucMaster(yRow,xCol,:);       ptSignl  = reshape(ptSignl,1,size(inputParams.flucMaster,3),[]);
       linSignl = inputParams.flucMaster(lRow,:,:);          linSignl = reshape(linSignl,size(inputParams.flucMaster,2),size(inputParams.flucMaster,3),[]);
       linSignl = linSignl(inputParams.xStrt:end,:);         % Omitting signals before X = 0
%  Creating a line-plane signal with values equal to the point signal
%      ptSignl  = repmat(ptSignl,size(linSignl,1),1);
%  Computing cross correlation
       ptLineCor = xcorr2(ptSignl,linSignl);   ptLineCor = ptLineCor/(max(max(ptLineCor)));  figure(1+ptCtr);
       pcolor(xNew,tLags,fliplr(ptLineCor'));         shading interp;    caxis([-1 1]);             colormap bluewhitered;       ylim([-0.001 0.001]);
       xlabel(inputParams.xName,'Interpreter','latex');          set(gca,'FontSize',13);    box off;       ax = gca;      
       ylabel('$\Delta \tau$','Interpreter','latex');            ax.TickLabelInterpreter = 'latex';        c = colorbar;       
       set(gcf,'Position',[870 70 620 470]);                     xlim([0 inputParams.limX(end)]);          c.TickLabelInterpreter = 'latex'; 
       title(c,'$C_{x -}$','Interpreter','latex');               ptCtr = ptCtr + 1;             
       title(['$\hat{R}_{(x,y)}:(',num2str(round(Xn(xCol),2)),',',num2str(round(Yn(yRow),2)),');\hat{L}_y=',num2str(round(Yn(lRow),2)),'$'],'Interpreter','latex');
       prompt = [newline,'->> Plot another set(y/n) - '];  chc = input(prompt,'s');
    end
end
%%      Function 6: POINT - VERTICAL   LINE CROSS CORRELATION
function pointVLineCor(inputParams)
    chc = 'y';                    tLags = inputParams.tLags;        ptCtr = 1;
    Xn = inputParams.Xn;          Yn = inputParams.Yn;        
    while strcmp(chc,'y')
%  Location of 1st point
       prompt = [newline,'->> Enter Point X-Location(val>0) - '];    inpX = input(prompt);
       prompt = [newline,'->> Enter Point Y-Location - '];           inpY = input(prompt);
       prompt = [newline,'->> Enter Line  X-Location - '];           inpLx = input(prompt);
%  Exact location of 1st point
       [xCol,~,yRow,~] = precisLoc(inpX,inpY,inputParams.Xn,inputParams.Yn);
%  Exact location of line
       [lCol,~,~,~]    = precisLoc(inpLx,0,inputParams.Xn,inputParams.Yn);
%  Plotting marker lines & point locations
       pointMrk(1,xCol,yRow,Xn,Yn,'r',ptCtr);
%  Plotting line
       xLoctn = zeros(size(Yn));          xLoctn(:) = Xn(lCol);
       figure(1);                         plot(xLoctn,Yn,'-.','Color',inputParams.colr(ptCtr,:),'LineWidth',1.2) ; 
%  Obtaining signals for point & line-plane(line though time)
       ptSignl  = inputParams.flucMaster(yRow,xCol,:);       ptSignl  = reshape(ptSignl,1,size(inputParams.flucMaster,3),[]);
       linSignl = inputParams.flucMaster(:,lCol,:);          linSignl = reshape(linSignl,size(inputParams.flucMaster,1),size(inputParams.flucMaster,3),[]);
%  Creating a line-plane signal with values equal to the point signal
%      ptSignl  = repmat(ptSignl,size(linSignl,1),1);
%  Computing cross correlation
       ptLineCor = xcorr2(ptSignl,linSignl);  ptLineCor = ptLineCor/(max(max(ptLineCor)));  figure(2+ptCtr);
       pcolor(Yn,tLags,fliplr(ptLineCor'));         shading interp;      caxis([-1 1]);             colormap bluewhitered;       ylim([-0.0008 0.0008]);
       xlabel(inputParams.yName,'Interpreter','latex');          set(gca,'FontSize',13);    box off;       ax = gca;      
       ylabel('$\Delta \tau$','Interpreter','latex');            ax.TickLabelInterpreter = 'latex';        c = colorbar;       
       set(gcf,'Position',[870 70 620 470]);                     xlim(inputParams.limY);    c.TickLabelInterpreter = 'latex'; 
       title(c,'$C_{x -}$','Interpreter','latex');               ptCtr = ptCtr + 1;               
       title(['$\hat{R}_{(x,y)}:(',num2str(round(Xn(xCol),2)),',',num2str(round(Yn(yRow),2)),';\hat{L}_x:',num2str(round(Xn(lCol),2)),')$'],'Interpreter','latex');
       prompt = [newline,'->> Plot another set(y/n) - '];  chc = input(prompt,'s');
    end
end
%%      Function 3.1: FIND PRECISE LOCATION
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
%%      Function 3.2: PLOTTING POINT & LOCATOIN MARKER LINES
function pointMrk(figNo,xCol,yRow,Xn,Yn,ptColr,ctr)   
%  Plotting marker lines
    yLoctn = zeros(size(Xn));            yLoctn(:) = Yn(yRow);
    figure(figNo);                       plot(Xn,yLoctn,'-.w','LineWidth',1.2) ; 
    xLoctn = zeros(size(Yn));            xLoctn(:,1) = Xn(xCol);           
    figure(figNo);                       plot(xLoctn,Yn,'-.w','LineWidth',1.2) ;          
%  Plotting point for reference
    figure(1);    plot(Xn(xCol),Yn(yRow),'o','MarkerFaceColor',ptColr,'MarkerEdgeColor','k','MarkerSize',6);
     if Yn(yRow) > 0
        adj = 0.2;
    else
        adj = -0.2;
    end
    text(Xn(xCol)+0.2,Yn(yRow)+adj,['$P:',num2str(ctr),'$'],'Color','k','FontSize',15,'Interpreter','latex');
end

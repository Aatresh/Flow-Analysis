%% COMPUTE UNCERTAINITY VALUES & STATISTICS CONVERGENCE
% Statistics are extracted from vector values & checked for uncertainity
% using Standard Deviation based confidence levels & for statistical convergence
% Reference: Lazar; A Practical Approach to PIV Uncertainity Analysis
% Validation Metrics: Confidence coefficient; Convergence coefficient
% confCoeff : Confidence Coefficient
% confCoeff = 1    -> 68% of all values lie within 1 Standard Deviation window
% confCoeff = 1.96 -> 95% of all values lie within 2 Standard Deviations
% covgCoeff : Level of convergence of statistics
%% Conditions
tic;    clearvars; clc;    fclose all;      set(0,'defaultfigurecolor',[1 1 1]); 
colr = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840] ;
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';  code = 'PIV_PODCalc';
cd([root1 'Jet_Analysis\Global_Functions\']); addpath([root1 'Jet_Analysis\PIV_Stiching_Codes\V2\']);
configTable;    disp([newline,'||--- Uncertainity & Convergence Computation ---||']);  
tests = {'NPR_2p5_TR_1p0', 'NPR_2p6_TR_1p0','NPR_2p9_TR_1p0','NPR_3p0_TR_1p0', 'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_5p0_TR_1p0',...
         'NPR_2p5_TR_1p1', 'NPR_3p0_TR_1p1', 'NPR_3p5_TR_1p1', 'NPR_3p6_TR_1p1', 'NPR_4p0_TR_1p1', 'NPR_4p5_TR_1p1',...
         'NPR_2p5_TR_2p0', 'NPR_3p0_TR_2p0', 'NPR_3p5_TR_2p0', 'NPR_3p6_TR_2p0', 'NPR_4p0_TR_2p0', 'NPR_4p5_TR_2p0',...
         'NPR_2p5_TR_2p6', 'NPR_3p0_TR_2p6','NPR_3p5_TR_2p6','NPR_3p6_TR_2p6', 'NPR_4p0_TR_2p6','NPR_4p5_TR_2p6'};

condition = tests(4);     nozzle = 'Major';      config = 'TR1';       
%% Input path & Image Filtering
[OutputStruct] = GF_DriveSelect(config,nozzle,code);     normFac = OutputStruct.dt;  
% Flow conditions
[~,Uj,NPR,NTR] = GF_Velocity(condition{1});
% Path for matrix files
driveMat  = [OutputStruct.in_root  condition{1}(9:14) '\' condition{1}(1:7) '\'];
% Path for velocity filters
driveFilt = [OutputStruct.bckgrnd_root condition{1}(9:14) '\' condition{1}(1:7) '\'];
% Loading data sets
load([driveMat 'U_vel']);       load([driveMat 'V_vel']);     
load([driveMat 'X_loc']);       load([driveMat 'Y_loc']);    
% Loading velocity filters & stitching/shift parameters 
velFil = load([driveFilt 'Velocity Filters_D.DAT']);
if contains(driveFilt,'Merged')
   load([driveFilt 'Img_Adjustments_Vx_' condition{1}(9:14) '_' condition{1}(1:7)]);
   imgStrt = Img_Adjustments(1);   imgEnd = Img_Adjustments(2);   imgShift = Img_Adjustments(3);
else
   load([driveFilt 'End_Values_AvgVx_' condition{1}(9:14) '_' condition{1}(1:7)]);
   imgStrt = Vals_TwoImg(1,1);     imgEnd = Vals_TwoImg(2,1);     imgShift = Vals_TwoImg(3,1);   zeroShift = Vals_TwoImg(1,3);
end
InputStruct = []; InputStruct.config = config;       %stchParm = Vals_TwoImg;   clear Vals_TwoImg;
% Assigning values for image filtering
InputStruct.Upos_C1 = velFil(1,1); InputStruct.Vpos_C1 = velFil(2,1); InputStruct.Vneg_C1 = velFil(3,1);   
InputStruct.Upos_C2 = velFil(1,2); InputStruct.Vpos_C2 = velFil(2,2); InputStruct.Vneg_C2 = velFil(3,2);
% Filtering images for both cameras
% camNos = size(U_vel,4);
for m = 1
  InputStruct.U = U_vel(:,:,:,m);   InputStruct.V = V_vel(:,:,:,m); InputStruct.n_cam = m;  
  [Master_U,Master_V,img_inf,bad_U,bad_V,U_Imgset,Vpos_imgset,Vneg_imgset] = GF_PIV_ImageSort(InputStruct);
  imgCounter(m) = size(Master_U,3);              %#ok
  disp([newline '->> Number of Images cut: ',num2str(size(U_vel,3)-imgCounter(m))]);
  disp([newline '>> Camera ',num2str(m),' Image Set after Sorting - ',num2str(imgCounter(m)),' Images' newline]);
end
% Adjusting image & axes based on stich/shift parameters
Master_U = Master_U(:,imgStrt:end-imgEnd,:);      Master_V = Master_V(:,imgStrt:end-imgEnd,:); 
% Reshaping & adjusting axes
X1 = reshape(X_loc(:,1),[size(U_vel,2),size(U_vel,1)]);    X1 = X1(:,1)';   
Y1 = reshape(Y_loc(:,1),[size(U_vel,2),size(U_vel,1)]);    Y1 = Y1(1,:)';
% Shifting & re-arranging axes according to parameters
X1 = X1(imgStrt:end-imgEnd);   X1 = X1 + abs(X1(1));       X1 = X1 + zeroShift;
if imgShift > 0
   Master_U = Master_U(1:end-imgShift,:,:);   Master_V = Master_V(1:end-imgShift,:,:);
   Y1 = Y1(1:end-imgShift);
else
   Master_U = Master_U(1-imgShift:end,:,:);   Master_V = Master_V(1-imgShift:end,:,:);
   Y1 = Y1(1-imgShift:end);
end
% Image size distribution
rowNos = size(Master_U,1);   colNos = size(Master_U,2);   imgNos = size(Master_U,3);
% Computing mean values
uBar = mean(Master_U,3);             vBar = mean(Master_V,3);
% Finding centerline
[jetEdges,Y] = GF_CenterFinder(uBar,Y1,config,nozzle,'Ux');
% Normalizing axes
Xn = X1/normFac;                     Yn = Y/normFac;
% Computing fluctuating components & additional statistics
uPrim = Master_U - uBar;             vPrim = Master_V - vBar;
uPrimSq = uPrim.^2;                  vPrimSq = power(vPrim,2);
TKE1 = ((mean(uPrimSq,3) + 2*mean(vPrimSq,3))/2)/Uj^2;     TKE2 = sqrt(TKE1*Uj^2)/Uj;
% Standard deviations of flucuting components
uPrimStd = std(uPrim,1,3);           vPrimStd = std(vPrim,1,3);
% Standard deviation of instantaneous velocity components
uInsStd = std(Master_U,1,3)/sqrt(imgNos);   vInsStd = std(Master_V,1,3)/sqrt(imgNos);
% Plot metrics - Axial
GF_FigurePlot(uInsStd,Xn,Yn);    title('$Axial \thinspace Instantaneous - \sigma_{u}$','Interpreter','latex');    set(gca,'FontSize',13);
xlabel('$X/D_e$','Interpreter','latex');    ylabel('$Y/D_e$','Interpreter','latex');    ax = gca;   ax.TickLabelInterpreter = 'latex';   c = colorbar;   c.TickLabelInterpreter = 'latex';
GF_FigurePlot(uPrimStd,Xn,Yn);   title('$Axial \thinspace Fluctuating - \sigma_{u''}$','Interpreter','latex'); set(gca,'FontSize',13);
xlabel('$X/D_e$','Interpreter','latex');    ylabel('$Y/D_e$','Interpreter','latex');    ax = gca;   ax.TickLabelInterpreter = 'latex';   c = colorbar;   c.TickLabelInterpreter = 'latex';
% Plot metrics - Radial
GF_FigurePlot(vInsStd,Xn,Yn);    title('$Radial \thinspace Instantaneous - \sigma_{v}$','Interpreter','latex'); set(gca,'FontSize',13);
xlabel('$X/D_e$','Interpreter','latex');    ylabel('$Y/D_e$','Interpreter','latex');    ax = gca;   ax.TickLabelInterpreter = 'latex';   c = colorbar;   c.TickLabelInterpreter = 'latex';
GF_FigurePlot(vPrimStd,Xn,Yn);    title('$Radial \thinspace Fluctuating - \sigma_{v''}$','Interpreter','latex'); set(gca,'FontSize',13);
xlabel('$X/D_e$','Interpreter','latex');    ylabel('$Y/D_e$','Interpreter','latex');    ax = gca;   ax.TickLabelInterpreter = 'latex';   c = colorbar;   c.TickLabelInterpreter = 'latex';
%% Temporary holding matrices to compute errors
uTemp = reshape(Master_U,rowNos*colNos,imgNos,[]);
vTemp = reshape(Master_V,rowNos*colNos,imgNos,[]);
%% Pick random image sets based on defined groups of image sets & confidence coefficient value
imgGrp = 5;    imgSets = imgGrp:imgGrp:1000;     zCoeff = 1.96; 
% Standard error of each set of random images exported to plot as .gif 
uErGif = zeros(rowNos,colNos,size(imgSets,2));   vErGif = zeros(rowNos,colNos,size(imgSets,2));
% Uncertainty of each set of random images exported to plot as .gif 
uUncerGif = zeros(rowNos,colNos,size(imgSets,2));   vUncerGif = zeros(rowNos,colNos,size(imgSets,2));
% Standard deviation & standard error matrices for each set of random images
uStdDv  = zeros(size(uTemp,1),size(imgSets,2));  vStdDv  = zeros(size(vTemp,1),size(imgSets,2));
uStdErr = zeros(size(uTemp,1),size(imgSets,2));  vStdErr = zeros(size(vTemp,1),size(imgSets,2));
% Peak & mean standard errors for each random image set
uStdErrPeak = zeros(1,size(imgSets,2));          vStdErrPeak = zeros(1,size(imgSets,2));
uStdErrMean = zeros(1,size(imgSets,2));          vStdErrMean = zeros(1,size(imgSets,2));
% Uncertainty matrices
uUncer = zeros(size(uTemp,1),size(imgSets,2));   vUncer = zeros(size(vTemp,1),size(imgSets,2));
% Uncertainty peaks & mean for each random image set 
uUncerPeak = zeros(1,size(imgSets,2));           vUncerPeak = zeros(1,size(imgSets,2));
uUncerMean = zeros(1,size(imgSets,2));           vUncerMean = zeros(1,size(imgSets,2));
% Peaks of Percentage deviation from mean values
uDevtnPeak = zeros(size(jetEdges,1),size(imgSets,2));
vDevtnPeak = zeros(size(jetEdges,1),size(imgSets,2));
tkeDevtnPeak = zeros(size(jetEdges,1),size(imgSets,2));
% Holding matrix for random set means
uRndMeanHold = zeros(size(uBar,1),size(uBar,2),size(imgSets,2));
vRndMeanHold = zeros(size(vBar,1),size(vBar,2),size(imgSets,2));
% Looping through & computing standard error for each image grouping
wB = waitbar(0,'Computing Error & Convergence Metrics');
for ctr = 1:size(imgSets,2)
% Picking a random set of images equal to the size of the image set N(5,10,15.. etc)
  N = imgSets(ctr);                  imgRange = randperm(imgNos,N);
  uRand = uTemp(:,imgRange);         vRand = vTemp(:,imgRange);
% Mean & fluctuating components
  uRndMean = mean(uRand,2);          uRndFluc = uRand - uRndMean;
  vRndMean = mean(vRand,2);          vRndFluc = vRand - vRndMean;
% TKE of random image set
  tkeRnd = 0.5*(mean(uRndFluc.^2,2) + 2*mean(vRndFluc.^2,2))/Uj^2;
% Standard deviation & standard error for the random image set
  uStdDv(:,ctr) = std(uRand,1,2);    uStdErr(:,ctr) = uStdDv(:,ctr)/sqrt(N);  
  vStdDv(:,ctr) = std(vRand,1,2);    vStdErr(:,ctr) = vStdDv(:,ctr)/sqrt(N);
% Exporting standard error to plor as .gif
  uErGif(:,:,ctr) = reshape(uStdErr(:,ctr),rowNos,colNos); 
  vErGif(:,:,ctr) = reshape(vStdErr(:,ctr),rowNos,colNos); 
% Peak & mean standard error
  uStdErrPeak(ctr) = max(uStdErr(:,ctr));           vStdErrPeak(ctr) = max(vStdErr(:,ctr)); 
  uStdErrMean(ctr) = mean(uStdErr(:,ctr),1);        vStdErrMean(ctr) = mean(uStdErr(:,ctr),1);
% Uncertainty/error of velocities based on 95% confidence interval 
  uUncer(:,ctr) = uStdErr(:,ctr)*zCoeff;            vUncer(:,ctr) = vStdErr(:,ctr)*zCoeff;
% Exporting standard error to plor as .gif
  uUncerGif(:,:,ctr) = reshape(uUncer(:,ctr),rowNos,colNos); 
  vUncerGif(:,:,ctr) = reshape(vUncer(:,ctr),rowNos,colNos); 
% Peak & mean uncertainty values across whole image
  uUncerPeak(ctr) = max(uUncer(:,ctr));             vUncerPeak(ctr) = max(vUncer(:,ctr));
  uUncerMean(ctr) = mean(uUncer(:,ctr));            vUncerMean(ctr) = mean(vUncer(:,ctr));
% Reshaping means for each random image set & storing in holding matrix
  uRndMean = reshape(uRndMean,rowNos,colNos);       vRndMean = reshape(vRndMean,rowNos,colNos);   tkeRnd = reshape(tkeRnd,rowNos,colNos);
  uRndMeanHold(:,:,ctr) = uRndMean;                 vRndMeanHold(:,:,ctr) = vRndMean;
% Digression Factor in Symmetry line & different shear layers
  uDevtn = zeros(size(jetEdges,1),size(uBar,2));    vDevtn = zeros(size(jetEdges,1),size(vBar,2));
  tkeDevtn = zeros(size(jetEdges,1),size(TKE1,2));
% Computing Digression Factor along centerline & symmetry lines from mean for each random image set
  for ctr1 = 1:size(uDevtn,1)
    uDevtn(ctr1,:)   = 100*((uBar(jetEdges(ctr1,1),:) - uRndMean(jetEdges(ctr1,1),:))./max(uBar(jetEdges(ctr1,1),:))); 
%     vDevtn(ctr1,:)   = 100*((vBar(jetEdges(ctr1,1),:) - vRndMean(jetEdges(ctr1,1),:))./vBar(jetEdges(ctr1,1),:));  
    tkeDevtn(ctr1,:) = 100*((TKE1(jetEdges(ctr1,1),:) - tkeRnd(jetEdges(ctr1,1),:))./TKE1(jetEdges(ctr1,1),:)); 
    % Peaks of % deviation
    uDevtnPeak(ctr1,ctr) = max(uDevtn(ctr1,:));     vDevtnPeak(ctr1,ctr) = max(vDevtn(ctr1,:));
    tkeDevtnPeak(ctr1,ctr) = max(tkeDevtn(ctr1,:));
  end
  clear uRndFluc uRndMean vRndFluc vRndMean;        waitbar(ctr/size(imgSets,2),wB);
end; close(wB);
% Note: Deviation/Digression factor for Vy is not straight forward to
% compute because of +/- values occuring in random image sets. Instead use
% tke to measure the combined error stemming from Vy as well
%% Plots & Results
% Peak Standard Error
figure(5);    subplot(121);    plot(imgSets,uStdErrPeak,'o','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k');   hold on;    grid on;
plot(imgSets,vStdErrPeak,'v','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','r');    box off;    c = legend('$v_x$','$v_y$');     c.EdgeColor = [1 1 1];    
xlabel(['$Image \thinspace Sets(\Delta n =',num2str(imgGrp),')$'],'Interpreter','latex');    set(gca,'FontSize',12);    ax = gca;         c.Interpreter = 'latex'; 
ylabel('${\dot \sigma}_{\varepsilon}(m/s)$','Interpreter','latex');    ax.TickLabelInterpreter = 'latex';  
title('$Peak \thinspace Velocity \thinspace Standard \thinspace Error$','Interpreter','latex'); 
% Mean Standard Error
subplot(122);    plot(imgSets,uStdErrMean,'o','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k');    hold on;        grid on;    box off;
plot(imgSets,vStdErrMean,'v','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','r');    c = legend('$v_x$','$v_y$');     c.EdgeColor = [1 1 1];   
xlabel(['$Image \thinspace Sets(\Delta n =',num2str(imgGrp),')$'],'Interpreter','latex');    set(gca,'FontSize',12);    ax = gca;         c.Interpreter = 'latex'; 
ylabel('${\bar \sigma}_{\varepsilon}(m/s)$','Interpreter','latex');    ax.TickLabelInterpreter = 'latex';  
title('$Mean \thinspace Velocity \thinspace Standard \thinspace Error$','Interpreter','latex'); 
set(gcf,'Position',[80 145 1235 410]); 
% Peak Uncertainty
figure(6);    subplot(121);    plot(imgSets,uUncerPeak,'o','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k');   hold on;    grid on;
plot(imgSets,vUncerPeak,'v','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','r');     box off;    c = legend('$v_x$','$v_y$');     c.EdgeColor = [1 1 1];    
xlabel(['$Image \thinspace Sets(\Delta n =',num2str(imgGrp),')$'],'Interpreter','latex');    set(gca,'FontSize',12);    ax = gca;         c.Interpreter = 'latex'; 
ylabel('$\dot \varepsilon(m/s)$','Interpreter','latex');    ax.TickLabelInterpreter = 'latex';  
title('$Peak \thinspace Velocity \thinspace Uncertainty$','Interpreter','latex'); 
% Mean Uncertainty
subplot(122);    plot(imgSets,uUncerMean,'o','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k');    hold on;        grid on;    box off;
plot(imgSets,vUncerMean,'v','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','r');     c = legend('$v_x$','$v_y$');     c.EdgeColor = [1 1 1];   
xlabel(['$Image \thinspace Sets(\Delta n =',num2str(imgGrp),')$'],'Interpreter','latex');    set(gca,'FontSize',12);    ax = gca;         c.Interpreter = 'latex'; 
ylabel('$\bar \varepsilon(m/s)$','Interpreter','latex');    ax.TickLabelInterpreter = 'latex';  
title('$Mean \thinspace Velocity \thinspace Uncertainty$','Interpreter','latex'); 
set(gcf,'Position',[80 145 1235 410]); 
% Peak Percentage Deviation - Digression Factor
plotPercentDevtn(uDevtnPeak,tkeDevtnPeak,imgSets,config,nozzle,imgGrp)
%% Generating GIFs for Error - Initializing new Struct
gIfGen = [];    gIfGen.Xn = Xn;    gIfGen.Yn = Yn;   
gIfGen.xLim = [0 5];    gIfGen.yLim = [-1.8 1.8];     gIfGen.imgSets = imgSets; 
gIfGen.driveOut = [root1 'PIV_Data\Validation\'];
gIfGen.driveCode = [root1 'Jet_Analysis\Global_Functions\'];
%% Axial Standard Error .gif
gIfGen.Input = uErGif;    gIfGen.gifname = ['Axial_StdError Evolution_',nozzle,'_',config];  
gIfGen.figTitle = '$\bar \sigma_{v_x} \thinspace Evolution[f_x(N_{images})]$';  
gIfGen.cBar = '$(m/s)$';  gIfGen.tag = 'StDeRR';      eRrGif(gIfGen); 
%% Radial Standard Error .gif
gIfGen.Input = vErGif;    gIfGen.gifname = ['Radial_StdError Evolution_',nozzle,'_',config]; 
gIfGen.figTitle = '$\bar \sigma_{v_y} \thinspace Evolution[f_x(N_{images})]$';  
gIfGen.cBar = '$(m/s)$';  gIfGen.tag = 'StDeRR';      eRrGif(gIfGen);  
%% Axial Uncertainty .gif
gIfGen.Input = uUncerGif; gIfGen.gifname = ['Axial_Uncertainty Evolution_',nozzle,'_',config];  
gIfGen.figTitle = '$\varepsilon_{v_x} \thinspace Evolution[f_x(N_{images})]$';  
gIfGen.cBar = '$(m/s)$';  gIfGen.tag = 'unCertn';     eRrGif(gIfGen); 
%% Radial Uncertainty .gif
gIfGen.Input = vUncerGif; gIfGen.gifname = ['Radial_Uncertainty Evolution_',nozzle,'_',config]; 
gIfGen.figTitle = '$\varepsilon_{v_y} \thinspace Evolution[f_x(N_{images})]$';  
gIfGen.cBar = '$(m/s)$';  gIfGen.tag = 'unCertn';     eRrGif(gIfGen);  
%% Convergence Rate - Use if needed
% chc = 'y';        ctr = 1;
% figure(10);    pcolor(Xn,Yn,uBar);    shading interp;    axis equal;    colormap jet;   ylim([-1.8 1.8]);    xlim([0 5]);    ax = gca;   
% ax.TickLabelInterpreter = 'latex';    c = colorbar;      set(ax,'FontSize',12);         c.TickLabelInterpreter = 'latex';    
% title('$Average \thinspace Axial \thinspace Velocity$','Interpreter','latex');          set(gcf,'Position',[55 405 620 410]);
% xlabel('$X/D_e$','Interpreter','latex');    ylabel('$Y/D_e$','Interpreter','latex');    title(c,'$\bar U_x$','Interpreter','latex');
% ax.TickLength = [0 0.025];    hold on;      set(gca,'Layer','top');
% while strcmp(chc,'y')
%   prompt = [newline,'->> Enter Point X-Location(val>0): '];    inpX = input(prompt);
%   prompt = '->> Enter Point Y-Location: ';                     inpY = input(prompt);
% % Exact location of 1st point
%   [xCol,xVal,yRow,yVal] = precisLoc(inpX,inpY,Xn,Yn);
% % Plotting marker lines & point locations
%   pointMrk(10,xCol,yRow,Xn,Yn,'g',ctr);                        meanConvgOrd = 0;                                               
% % Compute convergence rate for that location
%   for ctr2 = 1:size(imgSets,2)
%     uConvgRt = 100*abs(uRndMeanHold(yRow,xCol,ctr2) - uBar(yRow,xCol))/uBar(yRow,xCol);
%     figure(11);    plot(imgSets(ctr2),uConvgRt,'o','MarkerSize',6,'MarkerFaceColor',colr(ctr,:),'MarkerEdgeColor','k');    grid on;    hold on;
%     box off;   xlabel(['$Image \thinspace Sets(\Delta n =',num2str(imgGrp),')$'],'Interpreter','latex');    ax = gca;      ax.TickLabelInterpreter = 'latex';
%   end
%   prompt = [newline,'->> plot another point(y/n): '];          chc = input(prompt,'s');   ctr = ctr+1;
% end

%% Functions:
%%      Function 1: GIF Showing Error Development
function eRrGif(gIfGen)
fig_no = 50;     figure(fig_no);    figMat = gIfGen.Input;        imgNo = size(figMat,3); 
Xn = gIfGen.Xn;        Yn = gIfGen.Yn;      xname = '$X/D_e$';    yname = '$Y/h$'; 
xLim = gIfGen.xLim;    yLim = gIfGen.yLim;  dispTime = 0.08;      gifTime = 0.08;    imgSets = gIfGen.imgSets;
if strcmp(gIfGen.tag,'StDeRR')
   clim = [0 3.5];
elseif strcmp(gIfGen.tag,'unCertn')
   clim = [0 5];
end
  for ctr = 1:imgNo
    pcolor(Xn,Yn,figMat(:,:,ctr));       xlabel(xname,'Interpreter','latex');    ylabel(yname,'Interpreter','latex');   ax = gca;
    shading interp;   colormap jet;      axis equal;    caxis(clim);    c = colorbar;        title(c,gIfGen.cBar,'Interpreter','latex');
    c.TickLabelInterpreter = 'latex';    xlim(xLim);    ylim(yLim);     title(gIfGen.figTitle,'Interpreter','latex');             
    ax.TickLabelInterpreter = 'latex';   set(gca,'FontSize',12);          
    txT = text(3.5,1.5,['$N_{img}=',num2str(imgSets(ctr)),'$'],'Color','k','FontSize',12,'Interpreter','latex');
    txT.BackgroundColor = [1 1 1];       drawnow;       pause(dispTime);    
    frame1 = getframe(figure(fig_no));   im{ctr} = frame2im(frame1);     hold off;
  end
% Saving to Disk
prompt2 = [newline '>> Save .gif?(y/n) - '];  chc = input(prompt2,'s');   gifname = gIfGen.gifname;
if strcmp(chc,'y')
   cd(gIfGen.driveOut); 
   for idx = 1:size(im,2)
     [A,map] = rgb2ind(im{idx},256);
     if idx == 1
        imwrite(A,map,[gifname,'.gif'],'gif','LoopCount',Inf,'DelayTime',gifTime);
     else
        imwrite(A,map,[gifname,'.gif'],'WriteMode','append','DelayTime',gifTime);
     end
   end; cd(gIfGen.driveCode);
else
    return;
end
end
%%      Function 2: Plot the Digression Factors
function plotPercentDevtn(uDevtnPeak,tkeDevtnPeak,imgSets,config,nozzle,imgGrp)
if strcmp(config(1),'S') || strcmp(config(1),'C') || strcmp(nozzle,'Minor')
   lineNames = {'CL';'SL_{l}';'SL_{u}'};
else
   lineNames = {'SyL';'CL1';'CL2';'SL_{l1}';'SL_{u1}';'SL_{l2}';'SL_{u2}'};
end
for ctr = 1:length(lineNames)
  figure(ctr+6);    plot(imgSets,uDevtnPeak(ctr,:),'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k');    hold on;    grid on;    box off; 
  plot(imgSets,tkeDevtnPeak(ctr,:),'o','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','c');    box off;      ax = gca;   ax.TickLabelInterpreter = 'latex';           
  xlabel(['$Image \thinspace Sets(\Delta n =',num2str(imgGrp),')$'],'Interpreter','latex');            c = legend('$v_x$','$tke$');     c.EdgeColor = [1 1 1]; 
  ylabel('$\Omega_u, \Omega_{tke}(\%)$','Interpreter','latex');                                           c.Interpreter = 'latex';  ax.YScale = 'log';             
  title(['$Digression \thinspace Factor(\Omega):',lineNames{ctr},'$'],'Interpreter','latex');          set(gca,'FontSize',12);   set(gcf,'Position',[130 220 665 370]); 
end 
end
%%      Function 3: Find the precise location of given points
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
   [preY,~] = find(Yn<inpY);       [postY,~] = find(Yn>inpY);
   preDif = inpY - Yn(preY(1));    postDif = Yn(postY(end)) - inpY;
   if preDif < postDif
      yRow = preY(1);       yVal = round(Yn(preY(1)),2);
   else
      yRow = postY(end);    yVal = round(Yn(postY(end)),2);
   end
else
   yRow = y;   yVal = Yn(y);
end
end
%%      Function 4: Plotting point & Location Marker Lines
% function pointMrk(figNo,xCol,yRow,Xn,Yn,ptColr,ctr)   
% % Plotting marker lines
%   yLoctn = zeros(size(Xn));            yLoctn(:) = Yn(yRow);
%   figure(figNo);                       plot(Xn,yLoctn,'-.y','LineWidth',1.2) ; 
%   xLoctn = zeros(size(Yn));            xLoctn(:,1) = Xn(xCol);           
%   figure(figNo);                       plot(xLoctn,Yn,'-.y','LineWidth',1.2) ;          
% % Plotting point for reference
%   figure(figNo);    plot(Xn(xCol),Yn(yRow),'o','MarkerFaceColor',ptColr,'MarkerEdgeColor','k','MarkerSize',6);
%   if Yn(yRow) > 0
%      adj = 0.2;
%   else
%      adj = -0.2;
%   end
%   text(Xn(xCol)+0.2,Yn(yRow)+adj,['$P:',num2str(ctr),'$'],'Color','c','FontSize',15,'Interpreter','latex');
% end



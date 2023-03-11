%% POD ANALYSIS FOR SINGLE/TWIN JET 
% Computes the POD modes from vector files
clearvars; clc; fclose all; set(0,'defaultfigurecolor',[1 1 1]); code = 'PIV_PODCalc';
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\'; root2 = [root1 'PIV_Data\'];
cd([root1 'Jet_Analysis\Global_Functions\']); addpath([root1 'Jet_Analysis\POD_Codes\']); configTable;
tests = {'NPR_2p5_TR_1p0','NPR_2p6_TR_1p0','NPR_2p9_TR_1p0','NPR_3p0_TR_1p0', 'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0','NPR_5p0_TR_1p0', ...
         'NPR_2p5_TR_3p0','NPR_3p0_TR_3p0','NPR_3p6_TR_3p0','NPR_4p0_TR_3p0','NPR_4p5_TR_3p0'};

condition = [tests(3)];    nozzle = 'Minor';    config = 'S2';    verify = 'of';    testTKE = 'of';

[OutputStruct] = GF_DriveSelect(config,nozzle,code);   Deq = OutputStruct.dt;   InputStruct = [];
[Mj,Uj,NPR,TR] = GF_Velocity(condition{1});            InputStruct.config = config; 
% LOADING RAW DATA SETS
drive_in = [OutputStruct.in_root condition{1}(9:14) '\' condition{1}(1:7) '\'];
load([drive_in 'U_vel']); load([drive_in 'V_Vel']); load([drive_in 'X_loc']); load([drive_in 'Y_loc']);
X = reshape(X_loc(:,1),size(U_vel,2),size(U_vel,1));     X = X(:,1)';
%% DEFINING/LOADING IMAGE FILTERS
% Limits on U & V velocity components for image fitering : Load velocity filters if present 
clc; driveFilter = [OutputStruct.bckgrnd_root condition{1}(9:14) '\' condition{1}(1:7) '\'];
filt_name = 'Velocity Filters_D.DAT';             fid = fopen([driveFilter filt_name]);
if fid > -1
   disp('-> Loading Velocity Filters');   filtrs = textscan(fid,'%f %f','Delimiter',',');   camNo = 1; 
   InputStruct.Upos_C1  = filtrs{1}(1);   InputStruct.Vpos_C1 = filtrs{1}(2);  InputStruct.Vneg_C1 = filtrs{1}(3); 
   InputStruct.Upos_C2  = filtrs{2}(1);   InputStruct.Vpos_C2 = filtrs{2}(2);  InputStruct.Vneg_C2 = filtrs{2}(3);
else
   InputStruct.Upos_C1 = 2.5;   InputStruct.Vpos_C1 = 25;    InputStruct.Vneg_C1 = 21;    clc;
   InputStruct.Upos_C2 = [];    InputStruct.Vpos_C2 = [];    InputStruct.Vneg_C2 = [];   camNo = 1;
end
Y = reshape(Y_loc(:,1),size(U_vel,2),size(U_vel,1));     Y = Y(1,:)';   imgNos = size(U_vel,3);
%  IMAGE FILTERING
for m = camNo
  InputStruct.U = U_vel(:,:,:,m);   InputStruct.V = V_vel(:,:,:,m); InputStruct.n_cam = m;
  [Master_U,Master_V,img_inf,bad_U,bad_V,U_Imgset,Vpos_imgset,Vneg_imgset] = GF_PIV_ImageSort(InputStruct);
   img_counter(m) = size(Master_U,3);              %#ok
   disp([newline '>> Number of Images Cut-off - ',num2str(size(U_vel,3)-img_counter(m))]);
   disp([newline '>> Camera ',num2str(m),' Image Set After Sorting - ',num2str(img_counter(m)),' IMAGES' newline]);
%  PLOTTING THE IAMGE CUTOFF FOR U - VELOCITY
   figure; set(gcf,'Position',[70 420 1400 370])
   subplot(131);  plot(U_Imgset);   xlim([0 imgNos]);    title(['Ux peaks. Camera - ',num2str(m)]);
   box off;  grid on;  hold on; 
%  PLOTTING THE IMAGE CUTOFF FOR +ve V - VELOCITY
   subplot(132); plot(Vpos_imgset); xlim([0 imgNos]);  title(['Vy +ve peaks. Camera - ',num2str(m)]);
   box off;  grid on;
%  PLOTTING THE IMAGE CUTOFF FOR -ve V - VELOCITY
   subplot(133); plot(Vneg_imgset); xlim([0 imgNos]);  title(['Vy -ve peaks. Camera - ',num2str(m)]);
   box off;  grid on;
end
%% POD ANALYSIS
disp('--> Starting POD Analysis');          
%  ROW, COLUMN & IMAGE COUNT
rowNos = size(Master_U,1);      colNos = size(Master_U,2);   imgNos = size(Master_U,3);
%  RECALCULATING MEAN VELOCITY AFTER SORTING
U_bar = mean(Master_U,3);       V_bar = mean(Master_V,3);       
%  FLUCTUATING COMPONENTS
uPrim = Master_U - U_bar;       vPrim = Master_V - V_bar;
%  REARRANGING EACH IMAGE INTO ROW MATRIX
uTemp = permute(reshape(uPrim,1,rowNos*colNos,imgNos), [2,3,1]);
vTemp = permute(reshape(vPrim,1,rowNos*colNos,imgNos), [2,3,1]);
%  CREATING COVARIANCE MATRICES
covarU = (uTemp.'*uTemp);       covarV = (vTemp.'*vTemp);
%  SINGLE VALUE DECOMPOSITION
disp('-->> SVD Computation');          
[leftSinglrU,eigU] = svd(covarU);  [leftSinglrV,eigV] = svd(covarV);
%  EXTRACTING EIGEN VALUES - ENERGY CONTENT IN MODES
lambdaX = diag(eigU);           lambdaY = diag(eigV);
%  COMPUTING POD ENERGY MODES(NON - NORMALIZED)
phiXBase = uTemp*leftSinglrU;   phiYBase = vTemp*leftSinglrV;
%  NORMALIZING(2-NORM) SPATIAL MODES  
phiRecU = zeros(size(uTemp)); phiRecV  = zeros(size(vTemp));
for ctr = 1:imgNos
  phiRecU(:,ctr) = phiXBase(:,ctr)/norm(phiXBase(:,ctr),2);
  phiRecV(:,ctr) = phiYBase(:,ctr)/norm(phiYBase(:,ctr),2);
end; clear phiXBase phiYBase;
% TEMPORAL COEFFICIENT
tempUNorm = phiRecU'*uTemp;  tempVNorm = phiRecV'*vTemp; disp('   <- - Done! - ->');
%% RECONSTRUCTING IMAGES FROM POD MODES & COMPARING WITH SOURCE IMAGE TKEs FOR VERIFICATION
if strcmp(verify,'on')
% PERCENTAGE OF IMAGES TO PLOT
   percntImgs = round(0.5*imgNos);   wB = PoolWaitbar(imgNos, 'Reconstructing Approximation');
   nModes = 1:percntImgs; recUNorm = zeros(size(uTemp));    recVNorm = zeros(size(vTemp));
   for m = 1:imgNos
        recUNorm(:,m) = phiRecU(:,nModes)*tempUNorm(nModes,m);          
        recVNorm(:,m) = phiRecV(:,nModes)*tempVNorm(nModes,m);      
        increment(wB);
   end
   recUNorm = reshape(recUNorm,rowNos,colNos,[]);  uRecNorm = recUNorm + U_bar;
   recVNorm = reshape(recVNorm,rowNos,colNos,[]);  vRecNorm = recVNorm + V_bar;
%  COMPARING THE TKE OF RECONSTRUCTED IMAGES WITH BASE TKE
   baseTKE  = 0.5*(mean(uPrim.^2,3) + 2*mean(vPrim.^2,3));
   reconTKE = 0.5*(mean(recUNorm.^2,3) + 2*mean(recVNorm.^2,3));  diff = 100*((baseTKE - reconTKE)./baseTKE); filtTKE = 0.5*(baseTKE+reconTKE);
   GF_FigurePlot(baseTKE); fig1 = gcf; title('$Turbulence - Source \thinspace Images$','Interpreter','latex'); caxis([0 5000]);
   GF_FigurePlot(reconTKE);fig2 = gcf; title('$Turbulence - Reconstructed \thinspace Images$','Interpreter','latex'); caxis([0 5000]);
   GF_FigurePlot(diff);    fig3 = gcf; title('$\%Variation \thinspace in \thinspace TKE $','Interpreter','latex'); caxis([0 100]);
   GF_FigurePlot(filtTKE); fig4 = gcf; title('$Average \thinspace of \thinspace Source \thinspace \& \thinspace Reconstructed \thinspace TKE $','Interpreter','latex'); caxis([0 5000]);
   fig1.Position = [74 475 735 330];  fig2.Position = [74 100 735 330]; fig3.Position = [800 475 735 330]; fig4.Position = [800 100 735 330];
end
%% RE-NORMALIZING SPATIAL MODEA TO RANGE BETWEEN -1 TO 1
phiX = zeros(size(phiRecU));   phiY = zeros(size(phiRecV));
wB = waitbar(0,'Re-Normalizing Modes');
for ctr = 1:imgNos
% Axial Modes - Finding positive & negative values
  [posVal,~] = find(phiRecU(:,ctr)>0);   [negVal,~] = find(phiRecU(:,ctr)<0);
% Normalizing with peak values
  phiX(posVal,ctr) = phiRecU(posVal,ctr)/max(phiRecU(posVal,ctr));  
  phiX(negVal,ctr) = phiRecU(negVal,ctr)/abs(min(phiRecU(negVal,ctr)));
  clear posVal negVal;
% Radial Modes - Finding positive & negative values
  [posVal,~] = find(phiRecV(:,ctr)>0);   [negVal,~] = find(phiRecV(:,ctr)<0);
% Normalizing with peak values
  phiY(posVal,ctr) = phiRecV(posVal,ctr)/max(phiRecV(posVal,ctr));  
  phiY(negVal,ctr) = phiRecV(negVal,ctr)/abs(min(phiRecV(negVal,ctr)));
  clear posVal negVal; waitbar(ctr/imgNos,wB);
end; close(wB);
%% TEST PLOTTING MODES
prompt = ['--> Enter number of modes(1-',num2str(imgNos),')- ']; modeNo = input(prompt);
for ctr = 1:modeNo
  GF_FigurePlot(reshape(phiX(:,ctr),rowNos,colNos));fig1 = gcf; box off; colormap bluewhitered; 
  title(['$\bf\phi_x \thinspace No:',num2str(ctr),'$'],'Interpreter','latex'); fig1.Position = [100 430 560 290];
  GF_FigurePlot(reshape(phiY(:,ctr),rowNos,colNos));fig2 = gcf; box off; colormap bluewhitered; 
  title(['$\bf\phi_y \thinspace No:',num2str(ctr),'$'],'Interpreter','latex'); fig2.Position = [900 430 560 290];
end
%% PHASE PORTRAITS - Axial
mU1 = 1;     mU2 = 2;     tempMeanU = mean(sqrt(tempUNorm(mU1,:).^2+tempUNorm(mU2,:).^2));    
figure;plot(tempUNorm(mU1,1:2:end)/tempMeanU,tempUNorm(mU2,1:2:end)/tempMeanU,'s','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k'); box off;grid on;
mV1 = 2;     mV2 = 3;     tempMeanV = mean(sqrt(tempVNorm(mV1,:).^2+tempVNorm(mV2,:).^2));    
figure;plot(tempUNorm(mV1,1:2:end)/tempMeanV,tempVNorm(mV2,1:2:end)/tempMeanV,'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k'); box off;grid on;

%% PHASE PORTRAITS - Radial
mV1 = 1;     mV2 = 2;     tempMeanV = mean(sqrt(tempVNorm(mV1,:).^2+tempVNorm(mV2,:).^2));    
figure;plot(tempVNorm(mV1,1:2:end)/tempMeanV,tempVNorm(mV2,1:2:end)/tempMeanV,'o','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','k'); box off;grid on;
%% WRITING MODE DATA
prompt = [newline '>> WRITE MODE DATA(y/n)? - ']; write = input(prompt,'s');
drive_out = [OutputStruct.out_root condition{1}(9:14) '\' condition{1}(1:7) '\'];
if strcmp(write,'y')
   % AXIAL COMPONENTS
   save([drive_out 'modeEnergy-X'],'lambdaX');    save([drive_out 'PODModes-X'], 'phiX');
   save([drive_out 'tempCoeff-X'], 'tempUNorm');  save([drive_out 'ReconVar-X'], 'phiRecU');
   save([drive_out 'AvgVx'], 'U_bar');
   % RADIAL COMPONENTS
   save([drive_out 'modeEnergy-Y'],'lambdaY');    save([drive_out 'PODModes-Y'], 'phiY');
   save([drive_out 'tempCoeff-Y'], 'tempVNorm');  save([drive_out 'ReconVar-Y'], 'phiRecV');
   save([drive_out 'AvgVy'], 'V_bar');
end
%% TESTING TKE WEIGHTED AVERAGING
if strcmp(testTKE,'on')
   percntImgs = [0.2 0.75 0.85];   wB = waitbar(0,'Reconstructing Approximation');
   recUNorm = zeros(rowNos*colNos,imgNos,3);    recVNorm = zeros(rowNos*colNos,imgNos,3); 
   for ctr = 1:size(percntImgs,2)
     nModes = 1:round(percntImgs(ctr)*imgNos); 
     parfor m = 1:imgNos
        recUNorm(:,m,ctr) = phiRecU(:,nModes)*tempUNorm(nModes,m);      %#ok    
        recVNorm(:,m,ctr) = phiRecV(:,nModes)*tempVNorm(nModes,m);      %#ok
     end; waitbar(ctr/size(percntImgs,2),wB);
   end; close(wB);
   recUNorm = recUNorm(:,:,1)*1.1+recUNorm(:,:,2)+recUNorm(:,:,3); recUNorm = recUNorm/3;
   recVNorm = recVNorm(:,:,1)*1.1+recVNorm(:,:,2)+recVNorm(:,:,3); recVNorm = recVNorm/3;
   recUNorm = reshape(recUNorm,rowNos,colNos,[]);  uRecNorm = recUNorm + U_bar;
   recVNorm = reshape(recVNorm,rowNos,colNos,[]);  vRecNorm = recVNorm + V_bar;
%  COMPARING THE TKE OF RECONSTRUCTED IMAGES WITH BASE TKE
   baseTKE  = 0.5*(mean(uPrim.^2,3) + 2*mean(vPrim.^2,3));
   reconTKE = 0.5*(mean(recUNorm.^2,3) + 2*mean(recVNorm.^2,3));  diff = 100*((baseTKE - reconTKE)./baseTKE); 
   GF_FigPlot(baseTKE/Uj^2); fig1 = gcf; title('$Turbulence - Source \thinspace Images$','Interpreter','latex'); caxis([0 0.03]);
   GF_FigPlot(reconTKE/Uj^2);fig2 = gcf; title('$Turbulence - Reconstructed \thinspace Images$','Interpreter','latex'); caxis([0 0.03]);
   figure; contourf(diff,50); fig3 = gcf;colormap jet; title('$\%Variation \thinspace in \thinspace TKE $','Interpreter','latex'); caxis([0 100]); colorbar;
   fig1.Position = [74 475 735 330];  fig2.Position = [74 100 735 330]; fig3.Position = [800 475 735 330]; 
end
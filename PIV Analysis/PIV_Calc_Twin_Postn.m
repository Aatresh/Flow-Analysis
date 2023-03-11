%% CALCULATES THE FINAL PIV IMAGES BY FILTERING BAD VECTORS ACROSS IMAGES FOR TWO CAMERA POSITIONS
% Configuration Legend: Twin Square
% TS1  - TS-L50T16 : Computed{TwinCam} 
% TS2  - TS-L50M16 : Computed{Merged}  
% TS3  - TS-L50M16D: Davis{Merged}
clearvars; clc;    fclose all;      set(0,'defaultfigurecolor',[1 1 1]);
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';  code = 'PIV_VectCalc';
cd([root1 'Jet_Analysis\Global_Functions\']); addpath([root1 'Jet_Analysis\PIV_Stiching_Codes\V2\']);
disp([newline,'  Vector Calculation based on position of the Cameras']); configTable;
tests = {'NPR_2p4_TR_1p0', 'NPR_2p5_TR_1p0', 'NPR_2p6_TR_1p0', 'NPR_2p9_TR_1p0', 'NPR_3p0_TR_1p0', 'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_5p0_TR_1p0',...
         'NPR_2p5_TR_1p1', 'NPR_3p0_TR_1p1', 'NPR_3p5_TR_1p1', 'NPR_3p6_TR_1p1', 'NPR_4p0_TR_1p1', 'NPR_4p5_TR_1p1',...
         'NPR_2p5_TR_2p0', 'NPR_3p0_TR_2p0', 'NPR_3p5_TR_2p0', 'NPR_3p6_TR_2p0', 'NPR_4p0_TR_2p0', 'NPR_4p5_TR_2p0',...
         'NPR_2p5_TR_2p6', 'NPR_3p0_TR_2p6', 'NPR_3p5_TR_2p6', 'NPR_3p6_TR_2p6', 'NPR_4p0_TR_2p6', 'NPR_4p5_TR_2p6'};

condition = tests(6);     nozzle = 'Minor';      config = 'C4';   

[OutputStruct] = GF_DriveSelect(config,nozzle,code); 
count = 1;      imgCounter = [];       InputStruct = [];    
[~,Uj,NPR,NTR] = GF_Velocity(condition{1});
if strcmp(OutputStruct.bckgrnd_root(1),'D')
   type = {'AvgVx_1.DAT';'AvgVy_1.DAT';'TurbKineticE_1.DAT';'ReyStresXX_1.DAT';'ReyStresYY_1.DAT';'ReyStresXY_1.DAT'};
else
   type = {'AvgVx_2.DAT';'AvgVy_2.DAT';'TurbKineticE_2.DAT';'ReyStresXX_2.DAT';'ReyStresYY_2.DAT';'ReyStresXY_2.DAT'};
end
driveIn  = [OutputStruct.in_root  condition{1}(9:14) '\' condition{1}(1:7) '\'];
driveOut = [OutputStruct.out_root condition{1}(9:14) '\' condition{1}(1:7) '\'];  nozzle = OutputStruct.nozzle;
%  LOADING DATA
load([driveIn 'U_vel']);     load([driveIn 'V_vel']);
load([driveIn 'X_loc']);     load([driveIn 'Y_loc']);  imgNos = size(U_vel,3);
Flow = zeros(size(U_vel,1),size(U_vel,2),6*max(OutputStruct.dt));  InputStruct.config = config;
%% LIMITS ON U AND V VELOCITY FOR IMAGE FILTERING FOR CAMERA 1(C1) & CAMERA 2(C2) FOR EACH POSITION(D/U)
InputStruct.Upos_C1 = 2.0;  InputStruct.Vpos_C1 = 25;     InputStruct.Vneg_C1 = 30;   StdDev_Cutoff = 2.5;
InputStruct.Upos_C2 = 2.5;  InputStruct.Vpos_C2 = 40;     InputStruct.Vneg_C2 = 45;   clc;   tic;
%  FILTERING & COMPUTING FLOW METRICS
for m = OutputStruct.dt
%  ADDING OTHER INPUTS TO STRUCT FOR SORTING
   InputStruct.U = U_vel(:,:,:,m);   InputStruct.V = V_vel(:,:,:,m); InputStruct.n_cam = m;  
%  IMAGE SORTING FUNCTION
   [Master_U,Master_V,img_inf,bad_U,bad_V,U_Imgset,Vpos_imgset,Vneg_imgset] = GF_PIV_ImageSort(InputStruct);
   imgCounter(m) = size(Master_U,3);              %#ok
   disp([newline '->> Number of Images cut: ',num2str(size(U_vel,3)-imgCounter(m))]);
   disp([newline '>> Camera ',num2str(m),' Image Set after Sorting - ',num2str(imgCounter(m)),' Images' newline]);
%  PLOTTING THE IMAGE CUTOFF FOR U - VELOCITY
   figure; set(gcf,'Position',[70 420 1400 370])
   subplot(131);  plot(U_Imgset);   xlim([0 imgNos]);    title(['Ux peaks. Camera - ',num2str(m)]);
   box off;  grid on;  hold on; 
%  PLOTTING THE IMAGE CUTOFF FOR +ve V - VELOCITY
   subplot(132); plot(Vpos_imgset);   xlim([0 imgNos]);  title(['Vy +ve peaks. Camera - ',num2str(m)]);
   box off;  grid on;
%  PLOTTING THE IMAGE CUTOFF FOR -ve V - VELOCITY
   subplot(133); plot(Vneg_imgset);   xlim([0 imgNos]);  title(['Vy -ve peaks. Camera - ',num2str(m)]);
   box off;  grid on;
%  FLUCTUATING COMPONENT TO COMPUTE ReyStrXY
   U_prim = Master_U - mean(Master_U,3);         V_prim = Master_V - mean(Master_V,3);
%  STANDARD DEVIATION FILTERING
   [FilterOutput] = StdDev_Filter(Master_U,Master_V,StdDev_Cutoff);
%  TUBULENCE KINETIC ENERGY COMPUTATION
   TKE1 = sqrt(FilterOutput.TKE)/Uj;        
   FilterOutput.TKE      = FilterOutput.TKE/Uj^2; 
   FilterOutput.ReyStrXX = FilterOutput.ReyStrXX/Uj^2; 
   FilterOutput.ReyStrYY = FilterOutput.ReyStrYY/Uj^2; 
   FilterOutput.ReyStrXY = FilterOutput.ReyStrXY/Uj^2; 
%  REYNOLD'S STRESSES
   Rxy = (1/(size(U_prim,3)-1)) * sum(U_prim.*V_prim,3);
   if m == 1
%  CAMERA 1
      Flow(:,:,1) = FilterOutput.AvgUnew;   Flow(:,:,4) = FilterOutput.ReyStrXX; 
      Flow(:,:,2) = FilterOutput.AvgVnew;   Flow(:,:,5) = FilterOutput.ReyStrYY; 
      Flow(:,:,3) = FilterOutput.TKE;       Flow(:,:,6) = FilterOutput.ReyStrXY; 
   else
%  CAMERA 2
      Flow(:,:,7) = FilterOutput.AvgUnew;   Flow(:,:,10) = FilterOutput.ReyStrXX; 
      Flow(:,:,8) = FilterOutput.AvgVnew;   Flow(:,:,11) = FilterOutput.ReyStrYY; 
      Flow(:,:,9) = FilterOutput.TKE;       Flow(:,:,12) = FilterOutput.ReyStrXY; 
   end
   Axes = cat(2,X_loc(:,1),Y_loc(:,1),Y_loc(:,1));       [X,Y,~] = GF_TrixSort(Axes);
%  PLOTTING FLOW METRICS - TURBULENCE
   figure; set(gcf,'Position',[960 90 490 750]);
   subplot(311);  pcolor(X,Y,FilterOutput.TKE); axis equal;     shading interp;     colormap jet;   title('TKE - Type 1');
   colorbar;      caxis([0 max(max(FilterOutput.TKE))]);   xlim([min(X) max(X)]);  ylim([min(Y) max(Y)]);
%  PLOTTING FLOW METRICS - AXIAL VELOCITY
   subplot(312);  pcolor(X,Y,FilterOutput.AvgUnew); axis equal;   shading interp;     colormap jet;
   colorbar;     caxis([0 max(max(FilterOutput.AvgUnew))]);   xlim([min(X) max(X)]); ylim([min(Y) max(Y)]);
%  PLOTTING FLOW METRICS - RADIAL VELOCITY
   subplot(313);  pcolor(X,Y,FilterOutput.AvgVnew); axis equal;   shading interp;     colormap jet;
   colorbar;    caxis([min(min(FilterOutput.AvgVnew)) max(max(FilterOutput.AvgVnew))]);   xlim([min(X) max(X)]); ylim([min(Y) max(Y)]);
%  PLOTTING FLOW METRICS - TURBULENCE TYPE 2
   figure; pcolor(X,Y,TKE1);        axis equal;    shading interp;     colormap jet;  title('TKE - Type 2');
   colorbar;      caxis([0 0.5]);   xlim([min(X) max(X)]);  ylim([min(Y) max(Y)]);
%  PLOTTING REYNOLD'S STRESSES
   figure; set(gcf,'Position',[960 90 490 750]);
   subplot(311);  pcolor(X,Y,FilterOutput.ReyStrXX);   axis equal;   shading interp;     colormap jet;  title('Reynold''s Stress R_{xx}')
   colorbar;      caxis([0 max(max(FilterOutput.ReyStrXX))]);   xlim([min(X) max(X)]);  ylim([min(Y) max(Y)]);
   subplot(312);  pcolor(X,Y,FilterOutput.ReyStrYY);   axis equal;   shading interp;     colormap jet;  title('Reynold''s Stress R_{yy}')
   colorbar;      caxis([0 max(max(FilterOutput.ReyStrXX))]);   xlim([min(X) max(X)]);  ylim([min(Y) max(Y)]);
   subplot(313);  pcolor(X,Y,FilterOutput.ReyStrXY);   axis equal;   shading interp;     colormap jet;  title('Reynold''s Stress R_{xy}')
   colorbar;      xlim([min(X) max(X)]);  ylim([min(Y) max(Y)]);     toc;
end
%% WRITING TO FILE
commandwindow; prompt = ([newline '->> Save Data to Disk(y/n)?: ']); write = input(prompt,'s');
if strcmp(write,'y')
   cd(driveOut);  data = cell(6,size(Flow,3)/6);  matrixSize = size(Flow,1)*size(Flow,2);
%  DEFINING CELLS WITH ZERO MATRICES    
   for ctr = 1:size(data,1)*size(data,2)
     data{ctr} = zeros(size(Flow,1)*size(Flow,2),3);
   end
%  ADDING X & Y CO-ORDINATES TO CELL MATRICES
   for m = OutputStruct.dt
     for n = 1:6
       data{n,m}(:,1) = X_loc(:,m);   data{n,m}(:,2) = Y_loc(:,m);
     end
   end
%  ADDING U, V & TKE TO CELL MATRICES
   for count = 1:size(Flow,3) 
     data{count}(:,3) = reshape(Flow(:,:,count).',[],1);
   end
   m = max(OutputStruct.dt);
%  REARRANGING & WRITING TO .DAT FILE
   for x = 1:6
     if m>1
        tmat = zeros(matrixSize*m+4,3);
        tmat(4:matrixSize+3,:)    = data{x,1};
        tmat(matrixSize+5:end,:)  = data{x,2};
        dlmwrite(type{x},tmat,'precision','%.5f','delimiter',' ');
     else
        tmat = zeros(matrixSize*m+3,3);
        tmat(4:matrixSize+3,:)    = data{x,1};
        dlmwrite(type{x},tmat,'precision','%.5f','delimiter',' ');
     end
   end 
%  LOGGING THE FILTERING PARAMETERS & IMAGE VECTOR COUNT
   removd = imgNos - imgCounter;           limits = zeros(4,2);
   limits(1,1) = InputStruct.Upos_C1;      limits(1,2) = InputStruct.Upos_C2;
   limits(2,1) = InputStruct.Vpos_C1;      limits(2,2) = InputStruct.Vpos_C2;
   limits(3,1) = InputStruct.Vneg_C1;      limits(3,2) = InputStruct.Vneg_C2;
   limits(4,1) = size(FilterOutput.AvgUnew,1);            
   limits(4,2) = size(FilterOutput.AvgUnew,2);
%  WRITING ALL THE DATA TO DISK LOCATION
   if strcmp(OutputStruct.bckgrnd_root(1),'D')
      dlmwrite('Removd_Imgs_D.DAT',removd);   dlmwrite('Velocity Filters_D.DAT',limits);
   else
      dlmwrite('Removd_Imgs_U.DAT',removd);   dlmwrite('Velocity Filters_U.DAT',limits);
   end; disp('->> Done! <<-'); 
end
%% STANDARD DEVIATOIN FILTERING 
% Algorithm: 
%1. Check image set at every location for outliers based on standard deviation criteria
%2. Record position of all non outlier images in to matrices for average statistics computation
%3. Compute statistics based on position history obtained in above step
%4. Repeat process with a range of weighted standard deviations for Rxy computation
function [FilterOutput] = StdDev_Filter(U_temp,V_temp,StdDev_Cutoff)
disp([newline '<< -- Standard Deviation Filtering  -- >>']);
% Compute average,standard deviation
Avg_U = mean(U_temp,3);      StdDev_U = std(U_temp,0,3);
Avg_V = mean(V_temp,3);      StdDev_V = std(V_temp,0,3);
% Change format to 2D matrices
StdDev_U1 = reshape(StdDev_U,size(StdDev_U,1)*size(StdDev_U,2),1);
StdDev_V1 = reshape(StdDev_V,size(StdDev_V,1)*size(StdDev_V,2),1);
U_temp1 = reshape(U_temp,size(U_temp,1)*size(U_temp,2),size(U_temp,3));
V_temp1 = reshape(V_temp,size(V_temp,1)*size(V_temp,2),size(V_temp,3));
Avg_U1 = reshape(Avg_U,size(Avg_U,1)*size(Avg_U,2),1);
Avg_V1 = reshape(Avg_V,size(Avg_V,1)*size(Avg_V,2),1);
% Limits for standard deviation
  % Axial velocity limits
    diffPos_U = Avg_U1 + StdDev_Cutoff*StdDev_U1;
    diffNeg_U = Avg_U1 - StdDev_Cutoff*StdDev_U1;
  % Radial velocity limits
    diffPos_V = Avg_V1 + StdDev_Cutoff*StdDev_V1;
    diffNeg_V = Avg_V1 - StdDev_Cutoff*StdDev_V1;
% Loop for statistics computation
% New statistics after image filtering
AvgUnew = zeros(size(Avg_U1));          UprimSqMean = zeros(size(Avg_U1));
AvgVnew = zeros(size(Avg_V1));          VprimSqMean = zeros(size(Avg_V1));
ReyStrXX = zeros(size(Avg_U1));         %ReyStrXY = zeros(size(Avg_U1));
ReyStrYY = zeros(size(Avg_V1));         
posTracker_U = zeros(size(U_temp1));    posTracker_V = zeros(size(V_temp1));                    
% Axial velocity filtering
wB = waitbar(0,'Primary Statistics Filtering');
for posCounter  = 1:size(U_temp1,1)
  [~,imgNos] = find(U_temp1(posCounter,:) <= diffPos_U(posCounter) & U_temp1(posCounter,:) >= diffNeg_U(posCounter));
  % Average axial velocity at the point
  AvgUnew(posCounter) = mean(U_temp1(posCounter,imgNos));
  % Mean square of fluctuating component at the point
  UprimSqMean(posCounter) = mean((U_temp1(posCounter,imgNos) - AvgUnew(posCounter)).^2);
  % Compliant image tracking for each location
  posTracker_U(posCounter,1:size(imgNos,2)) = imgNos; 
  Uprim1 = U_temp1(posCounter,imgNos) - AvgUnew(posCounter); 
  % Reynolds stress of axial fluctuating component
  ReyStrXX(posCounter) = (1/(size(imgNos,2)-1))*sum(Uprim1.^2);     clear imgNos;
  % Radial velocity filtering
  [~,imgNos] = find(V_temp1(posCounter,:) <= diffPos_V(posCounter) & V_temp1(posCounter,:) >= diffNeg_V(posCounter));
  % Average axial velocity at the point
  AvgVnew(posCounter) = mean(V_temp1(posCounter,imgNos));
  % Mean square of fluctuating component at the point
  VprimSqMean(posCounter) = mean((V_temp1(posCounter,imgNos) - AvgVnew(posCounter)).^2);
  % Compliant image tracking for each location
  posTracker_V(posCounter,1:size(imgNos,2)) = imgNos;
  Vprim1 = V_temp1(posCounter,imgNos) - AvgVnew(posCounter);
  % Reynolds stress of radial fluctuating component
  ReyStrYY(posCounter) = (1/(size(imgNos,2)-1))*sum(Vprim1.^2); waitbar(posCounter/size(U_temp1,1),wB);
end; close(wB);
% Turbulence kinetic energy computation
TKE = (UprimSqMean + 2*VprimSqMean)/2;
% Generate a range of standard deviations and average Rxy over that range
StdDevRange = linspace(StdDev_Cutoff-0.5,StdDev_Cutoff*10,20);
% Weighting window to apply over Standard Deviation Range. Rxy calculated
% at lower levels get higher weighting
window = 1.35 - 0.5*cos((pi/4)*(0:size(StdDevRange,2)-1)/(size(StdDevRange,2)-1));
window = fliplr(window);
ReyStrXY    = zeros(size(Avg_U1,1),size(StdDevRange,2)); wB = waitbar(0,'R_{xy} Computation');
for ctr = 1:size(StdDevRange,2)
  % Limits for standard deviation
  % Axial velocity limits
    diffPos_U = Avg_U1 + StdDevRange(ctr)*StdDev_U1;
    diffNeg_U = Avg_U1 - StdDevRange(ctr)*StdDev_U1;
  % Radial velocity limits
    diffPos_V = Avg_V1 + StdDevRange(ctr)*StdDev_V1;
    diffNeg_V = Avg_V1 - StdDevRange(ctr)*StdDev_V1;
  for posCounter  = 1:size(U_temp1,1)
      [~,imgNos] = find(U_temp1(posCounter,:) <= diffPos_U(posCounter) & U_temp1(posCounter,:) >= diffNeg_U(posCounter));
     % Average axial velocity at the point
     AvgUnew(posCounter) = mean(U_temp1(posCounter,imgNos));
     % Compliant image tracking for each location
     Uprim1 = U_temp1(posCounter,imgNos) - AvgUnew(posCounter);   clear imgNos;
     % Radial velocity filtering
     [~,imgNos] = find(V_temp1(posCounter,:) <= diffPos_V(posCounter) & V_temp1(posCounter,:) >= diffNeg_V(posCounter));
      % Average axial velocity at the point
     AvgVnew(posCounter) = mean(V_temp1(posCounter,imgNos));
     Vprim1 = V_temp1(posCounter,imgNos) - AvgVnew(posCounter);
     % Compare sizes to use same number of images across U' & V'
     if size(Uprim1,2) < size(Vprim1,2)
        Vprim1 = Vprim1(1,1:size(Uprim1,2));   NoImgs = size(Uprim1,2);
     elseif size(Vprim1,2) < size(Uprim1,2)
        Uprim1 = Uprim1(1,1:size(Vprim1,2));   NoImgs = size(Vprim1,2);
     else
        NoImgs = size(imgNos,2);
     end
     % Reynold's stress for cross - component fluctuations
     ReyStrXY(posCounter,ctr) = window(ctr)*(1/(NoImgs-1))*sum(Uprim1.*Vprim1);
  end; waitbar(ctr/size(StdDevRange,2),wB);
end; close(wB);
ReyStrXY = mean(ReyStrXY,2);
% Reshaping all matrices
FilterOutput.AvgUnew     = reshape(AvgUnew,size(Avg_U,1),size(Avg_U,2));
FilterOutput.AvgVnew     = reshape(AvgVnew,size(Avg_U,1),size(Avg_U,2));
FilterOutput.UprimSqMean = reshape(UprimSqMean,size(Avg_U,1),size(Avg_U,2));
FilterOutput.VprimSqMean = reshape(VprimSqMean,size(Avg_U,1),size(Avg_U,2));
FilterOutput.ReyStrXX    = reshape(ReyStrXX,size(Avg_U,1),size(Avg_U,2));
FilterOutput.ReyStrYY    = reshape(ReyStrYY,size(Avg_U,1),size(Avg_U,2));
FilterOutput.ReyStrXY    = reshape(ReyStrXY,size(Avg_U,1),size(Avg_U,2));
FilterOutput.TKE         = reshape(TKE,size(Avg_U,1),size(Avg_U,2));  
end

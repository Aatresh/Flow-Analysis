%% STICHTING CODE TO COMBINE IMAGES FROM POSITION - 1
clc;clearvars;fclose all; set(0,'defaultfigurecolor',[1 1 1]); code = 'PIV_Stitch';
root = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';
cd([root 'Jet_Analysis\Global_Functions\']); addpath([root 'Jet_Analysis\PIV_Stiching_Codes\V2\']);
tests = {'NPR_2p4_TR_1p0','NPR_2p5_TR_1p0','NPR_2p6_TR_1p0','NPR_2p9_TR_1p0','NPR_3p0_TR_1p0', 'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0','NPR_5p0_TR_1p0', ...
         'NPR_2p5_TR_3p0','NPR_3p0_TR_3p0','NPR_3p6_TR_3p0','NPR_4p0_TR_3p0','NPR_4p5_TR_3p0'};
       
condition = tests(1);   nozzle = 'Major';   config = 'C9';   resTyp = 1;   datSave = 'n';  configTable;

%% IMAGE SMOOTHING CONDITION & VALUE CORRECTION 
smoothing = 0;      cVal1 = 0.000;       cVal2 = 0.000;
colRange  = 2; %No. of Colums to be smoothed    
smothCol  = 2; %No. of Colums used for smoothing 
%  FOLDER SELECTION 
[OutputStruct] = GF_DriveSelect(config,nozzle,code);   Deq = OutputStruct.dt;  
driveIn = [OutputStruct.in_root condition{1}(9:14) '\' condition{1}(1:7) '\'];
% DISTANCE FROM NOZZLE EXIT(in mm)
dist_exit = 0.7;
%%  IMAGE CUT & SHIFT PARAMETERS
%  for y.. (-) shifts image up.. (+) shifts image down
x1_del1 = 6;   x1_del2 = 7;    y1_del1 = 0;
  
x2_del1 = 5;   x2_del2 = 5;   y2_del1 = 0;
%% MATRIX SIZE & RESULT TYPE SELECTION
[type,col,matrixSize,clim,txtFormat,lineSkip,saveName] = Matrix_Size(driveIn,resTyp);
%% IMAGE LOAD & PREP - POSITION 1
camPostn = 1;               
% Load main data file & split into respective cameras
textFilename = [type,'_',num2str(camPostn), '.dat'];   fid = fopen([driveIn textFilename]);   
dataMain =  textscan(fid,txtFormat,'headerlines',lineSkip,'delimiter',' ');   
dataMain = cat(2,dataMain{1},dataMain{2},dataMain{3});
dataCam1 = dataMain(1:matrixSize,:);     dataCam2 = dataMain(matrixSize+2:end,:);
%  PREPARING FIRST CAMERA IMAGE
x1 = dataCam1(:,1);           y1 = dataCam1(:,2);                   
x_co1 = x1(1:col,1);          y_co1 = y1(1:col:matrixSize,1);     clear x1 y1;       
z1 = dataCam1(:,3)-cVal1;     zCam1 = vec2mat(z1,col);               clear z1;
%  REARRANGING AXES, CROPPING & ROTATING IMAGE
xAxisC1 = x_co1';                           x1_range = x1_del1:size(xAxisC1,2)-x1_del2;   
xAxisC1 = xAxisC1(x1_range);                
if xAxisC1(1) < 0
   xAxisC1 = xAxisC1 + abs(xAxisC1(1));
else
   xAxisC1 = xAxisC1 - xAxisC1(1);
end
zCam1 = circshift(zCam1,[y1_del1,0]);       xResCam1 = xAxisC1(3) - xAxisC1(2); %Axial Resolution - Cam:1
%  PREPARING SECOND CAMERA IMAGE
x2 = dataCam2(:,1);           y2 = dataCam2(:,2);       
x_co2 = x2(1:col,1);          y_co2 = y2(1:col:matrixSize,1);     clear x2 y2;                
z2 = dataCam2(:,3)-cVal2;     zCam2 = vec2mat(z2,col);               clear z2;
%  REARRANGING AXES, CROPPING & ROTATING IMAGE
xAxisC2 = x_co2';                           x2_range = x2_del1:size(xAxisC2,2)-x2_del2;           
xAxisC2 = xAxisC2(x2_range); 
%  RESOLUTION CHECK & AXIAL RE-SCALING
[xAxisC2] = resolCheck(xAxisC2,xResCam1,xAxisC1(end));
zCam2 = circshift(zCam2,[y2_del1,0]);           
%  TAKING MAX & MIN OF Y - SHIFT FOR VERTICAL CROP FROM BOTH CAMERA IMAGES         
if min([y1_del1 y2_del1])>0
   low = 0;
else    
   low = min([y1_del1 y2_del1]);
end
high = max([y1_del1 y2_del1]);              y1_range = 1-low:size(zCam1,1)-high;        
Z1 = zCam1(1-low:end-high,x1_range);        Z2 = zCam2(1-low:end-high,x2_range);
%  MERGING BOTH CAMERA AXES & ADDING THE NOZZLE DISTANCE CORRECTION
xFull = [xAxisC1,xAxisC2];  yFull = y_co1(y1_range);   xFull = xFull + dist_exit;
%  MERGING BOTH CAMERA IMAGES
ZFull = cat(2,Z1,Z2);
disp([newline '>> Vector Field Resolution(RxC): ',num2str(size(ZFull,1)),' X ',num2str(size(ZFull,2))]);
%% IMAGE SMOOTHING BY AVERAGING THE BOUNDARY BETWEEN IMAGES
if smoothing == 1
   temp1 = zeros(size(ZFull,1),colRange*2);
for i = 1:size(ZFull,1)
  avg1  = mean(diff([Z1(i,end-colRange:end) Z2(i,1:colRange)]));
  for ctr = 1:colRange*2
    temp1(i,ctr) =  Z1(i,end-smothCol)+avg1*ctr;
  end
end
   ZFull(:,size(Z1,2)-(colRange-1):size(Z1,2)+colRange)=temp1;
else
   colRange = NaN;    smothCol = NaN;
end
%% PLOT STICHED IMAGE
figure(1000);  pcolor(xFull/Deq,yFull/Deq,ZFull); set(gcf,'Position',[33 244 1280 443]);
shading interp;  colormap(jet);  colorbar; caxis(clim);  xlim([0 max(xFull/Deq)]);
% figure(2000);  pcolor(xFull,yFull,ZFull); set(gcf,'Position',[33 244 1280 443]);
% shading interp;  colormap(jet);  colorbar; caxis(clim);  xlim([0 max(xFull)]);
%  PLOT CENTERLINE
Plot_Center(xFull,yFull,ZFull,resTyp,Deq,driveIn);
%% VALUES OF CUT LENGTHS FOR EACH IMAGE
Vals_TwoImg = [];
Vals_TwoImg(:,1) = [x1_del1,x1_del2,y1_del1];      %First image correction
Vals_TwoImg(:,2) = [x2_del1,x2_del2,y2_del1];      %Second image correction
Vals_TwoImg(:,3) = dist_exit;                      %Distance from nozzle exit
Vals_TwoImg(:,4) = [cVal1,cVal2,NaN];              %Magnitude Correction
Vals_TwoImg(:,5) = [smoothing,colRange,smothCol];  %Smoothing toggle & colummn ranges 
%% SAVE DATA TO DISK
if strcmp(datSave,'y')
   saveX = repmat(xFull',size(yFull,1),1);           N = size(xFull,2);
   saveY = yFull(repmat(1:size(yFull,1),N,1),:);     saveZ = reshape(ZFull.',[],1);
   saveas(figure(1000),[driveIn condition{1}(1:7) '_' condition{1}(9:14) '_' nozzle '_' type '.jpg']);
   saveAll = [saveX saveY saveZ];
   save([driveIn condition{1}(1:7) '_' condition{1}(9:14) '_' nozzle '_' type '.dat'],'saveAll','-ascii');
%  SAVING DATA FOR EACH CONDITION
   file_out_data = [saveName condition{1}(9:14) '_' condition{1}(1:7)];
   save([driveIn file_out_data],'Vals_TwoImg');
end 
%% MATRIX SIZE FUNCTION
function [type,col,matrixSize,clim,txtFormat,lineSkip,saveName] = Matrix_Size(driveIn,resTyp)
resultNumb = [1;2;3;4;5;6];     
resultName = {'AvgVx';'AvgVy';'TurbKineticE';'ReyStresXX';'ReyStresYY';'ReyStresXY'};
inDx = find(resultNumb == resTyp);    type = resultName{inDx}; saveName = ['End_Values_',resultName{inDx},'_'];
if contains(driveIn,'Davis')  
   fid_mat = fopen([driveIn type '.DAT']);   mat_siz = textscan(fid_mat,'%s%s%s%s%s',5); 
   col = mat_siz{4}(4); row = mat_siz{5}(4);  clear fid_vel;  txtFormat = '%n%n%n%n%n'; 
   col = str2double(col{1}(3:end-1)); row = str2double(row{1}(3:end)); lineSkip = 3;
   climLwVal  = [0;-30;0;0;0;-1500];    climUpVal  = [550;30;5500;5500;2500;1500];
else
   fid_mat = fopen([driveIn 'Velocity Filters_D.DAT']);  mat_siz = textscan(fid_mat,'%d%d','headerlines',3,'Delimiter',',');
   row = mat_siz{1};  col = mat_siz{2};   txtFormat = '%n%n%n';  lineSkip = 3;
   climLwVal  = [0;-30;0;0;0;-0.01];    climUpVal  = [550;30;0.03;0.03;0.015;0.01];
end
clim = [climLwVal(inDx) climUpVal(inDx)];           matrixSize = row*col;
end
%% CENTERLINE VALUE PLOT
function Plot_Center(x_full,y_full,Z_full,result_typ,Deq,driveIn)
x = 1;          
strtPnt    = [5;10;5;5;5;100];
resultNumb = [1;2;3;4;5;6]; 
if contains(driveIn,'Davis') 
   limitVals = [250;10;800;800;400;900];
else
   limitVals = [250;5;0.005;0.01;0.001;0.003];
end
inDx      = resultNumb == result_typ;
chckLimit = limitVals(inDx);
strtLimit = strtPnt(inDx);
while isnan(Z_full(x,strtLimit))
    x = x+1;
end
while Z_full(x,strtLimit) < chckLimit 
    x = x+1;
end;       jetTop = x;    
while Z_full(x,strtLimit) > chckLimit
    x = x+1;
end;       jetBot = x-1;
jetCentr = round((jetTop+jetBot)/2);   locTn1 = num2str(round(y_full(jetCentr)/Deq,2)); locTn2 = num2str(round(y_full(jetCentr),2));
figure(2); subplot(211);plot(x_full/Deq,Z_full(jetCentr,:),'LineWidth',1.2); xlim([0 max(x_full/Deq)]); grid on; hold on;
title(['$Location-',locTn1,'$'],'Interpreter','latex');
subplot(212); plot(x_full,Z_full(jetCentr,:),'LineWidth',1.2); xlim([0 max(x_full)]); title(['$Location-',locTn2,'$'],'Interpreter','latex'); 
set(gcf,'Position',[100 100 1095 650]); grid on; hold on;
end
%% CHECK IF AXIAL RESOLUTION MATCHES FIRST CAMERA
function [newAxis] = resolCheck(xAxis,xResCam1,endVal)
if xAxis(1) > 0
   xAxis = xAxis - xAxis(1);
else
   xAxis = xAxis + abs(xAxis(1));
end
xRes = xAxis(3) - xAxis(2);
if xRes > xResCam1
   % Re-creating the axial co-ordinate with same axial resolution as Camera - 1
   newAxis = 0:xResCam1:xAxis(end);      newAxis = newAxis + (endVal+xResCam1);
else
   % Keep the axial resolution same & add term for length continuation
   newAxis = xAxis + (endVal+xResCam1);
end
end
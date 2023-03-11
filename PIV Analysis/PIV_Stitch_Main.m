%% CODE FOR STITCHING PIV IMAGES FOR ALL FOUR CAMERAS
clc;    clearvars;  fclose all; set(0,'defaultfigurecolor',[1 1 1]); code = 'PIV_Stitch';
root = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';
cd([root 'Jet_Analysis\Global_Functions\']); addpath([root 'Jet_Analysis\PIV_Stiching_Codes\V2\']);
tests = {'NPR_2p5_TR_1p0','NPR_2p9_TR_1p0','NPR_3p0_TR_1p0', 'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0','NPR_5p0_TR_1p0', ...
         'NPR_2p5_TR_3p0','NPR_3p0_TR_3p0','NPR_3p6_TR_3p0','NPR_4p0_TR_3p0','NPR_4p5_TR_3p0'};
       
condition = tests(3);   config = 'TS-L50T16';   nozzle = 'Major';   result_typ = 3;   data_save = 'y';

%  FOLDER SELECTION 
[OutputStruct] = GF_Drive_Select(config,nozzle,code);       Deq = OutputStruct.dt;  
drive_in = [OutputStruct.in_root condition{1}(9:14) '\' condition{1}(1:7) '\'];
% DISTANCE FROM NOZZLE EXIT(in mm)
dist_exit = -1;
%%  IMAGE CUT & SHIFT PARAMETERS
%  for y.. (-) shifts image up.. (+) shifts image down
x1_del1 = 1;   y1_del1 = 0;
x1_del2 = 10;  

x2_del1 = 15;   y2_del1 = 0;
x2_del2 = 20;  

x3_del1 = 10;   y3_del1 = -5;
x3_del2 = 1;  

x4_del1 = 10;   y4_del1 = -5;
x4_del2 = 15;  
%%  SMOOTHING ON(1)/OFF(0)
smoothing = 1;
%  AXIS INITIAL POINT FROM NOZZLE EXIT
start_mm = dist_exit+x1_del1*0.63;
%  TKE CORRECTION PER IMAGE
val00 = 0;
val0  = 0;
val1  = 0.003;
val2  = 0.002;
%  TYPE OF RESULT BEING STITCHED
if result_typ == 1
   type = 'AvgVx';
elseif result_typ == 2
   type = 'AvgVy';
elseif result_typ == 3
   type = 'TurbKineticE';
elseif result_typ == 4
   type = 'ReyStresXX';
elseif result_typ == 5
   type = 'ReyStresYY';     
else
   type = 'ReyStresXY';
end
%% MATRIX SIZE - BOTTOM TWO CAMERAS
position = 1;
[row1,col1,matrix_size1,~,~,~] = Matrix_Size(drive_in,config,result_typ,position);

% fid_mat = fopen([drive_in 'Velocity Filters_D.DAT']);     mat_siz = textscan(fid_mat,'%d%d','headerlines',3,'Delimiter',',');
% row1 = mat_siz{1};  col1 = mat_siz{2};    matrix_size1 = row1*col1;    clear fid_mat;
% IMAGE LOAD & PREP - POSITION 1 - BOTTOM TWO CAMERAS

textFilename = [type,'_',num2str(position), '.dat'];
fid=fopen([drive_in textFilename]);
data1=textscan(fid,'%n%n%n%n',matrix_size1,'headerlines',3,'delimiter',' ');
fid=fopen([drive_in textFilename]);
data2=textscan(fid,'%n%n%n%n',matrix_size1,'headerlines',matrix_size1+4,'delimiter',' ');
data1 = cat(2,data1{1},data1{2},data1{3});
data2 = cat(2,data2{1},data2{2},data2{3});
%  PREPARE FIRST IMAGE
x1 = data1(:,1);                                        x_co1 = x1(1:col1,1);
x_length1 = abs(x1(1)-x1(col1));                        y1 = data1(:,2);
y_co1 = y1(1:col1:matrix_size1,1);                      z1 = data1(:,3)-val00;
Z1a = vec2mat(z1,col1);                                 x1_range = x1_del1+1:size(Z1a,2)-x1_del2;
Z1a = circshift(Z1a,[y1_del1,0]);                     
%  PREPARE SECOND IMAGE
x2 = data2(:,1);                                        x_co2 = x2(x2_del1:col1,1);
x_length2 = abs(x2(1)-x2(col1-x2_del1));                y2 = data2(:,2);
y_co2 = y2(1:col1:matrix_size1,1);                      z2 = data2(:,3)-val0;
Z2a = vec2mat(z2,col1);                                 x2_range = x2_del1+1:size(Z2a,2)-x2_del2;
Z2a = circshift(Z2a,[y2_del1,0]);
%% MATRIX SIZE - TOP TWO CAMERAS
position = 2;
[row2,col2,matrix_size2,clim,~,~] = Matrix_Size(drive_in,config,result_typ,position);

% fid_mat = fopen([drive_in 'Velocity Filters_U.DAT']);     mat_siz = textscan(fid_mat,'%d%d','headerlines',3,'Delimiter',',');
% row2 = mat_siz{1};  col2 = mat_siz{2};    matrix_size2 = row2*col2;
%  IMAGE LOAD & PREP - POSITION - 2 (TOP IMAGE SETS)

textFilename = [type,'_',num2str(position), '.dat'];
fid=fopen([drive_in textFilename]);
data3=textscan(fid,'%n%n%n%n',matrix_size2,'headerlines',3,'delimiter',' ');
fid=fopen([drive_in textFilename]);
data4=textscan(fid,'%n%n%n%n','headerlines',matrix_size2+4,'delimiter',' ');
data3 = cat(2,data3{1},data3{2},data3{3});
data4 = cat(2,data4{1},data4{2},data4{3});
%  PREPARE THIRD IMAGE
x3 = data3(:,1);                                        x_co3 = x3(x3_del1:col2,1);
y3 = data3(:,2);                                        y_co3 = y3(1:col2:matrix_size2,1);
z3 = data3(:,3)-val1;                                   Z3a = vec2mat(z3,col2);
x3_range = x3_del1+1:size(Z3a,2)-x3_del2;               Z3a = circshift(Z3a,[y3_del1,0]);
%  PREPARE FOURTH IMAGE
x4 = data4(:,1);                                        x_co4 = x4(x4_del1:col2,1);
y4 = data4(:,2);                                        y_co4 = y4(1:col2:matrix_size2,1);
z4 = data4(:,3)-val2;                                   Z4a = vec2mat(z4,col2);
x4_range = x4_del1+1:size(Z4a,2)-x4_del2;               Z4a = circshift(Z4a,[y4_del1,0]);
if min([y1_del1 y2_del1 y3_del1 y4_del1])>0
   low = 0;
else    
   low = min([y1_del1 y2_del1 y3_del1 y4_del1]);
end
high = max([y1_del1 y2_del1 y3_del1 y4_del1]);          y1_range = 1-low:size(Z1a,1)-high;
Z1 = Z1a(1-low:end-high,x1_range);                      Z2 = Z2a(1-low:end-high,x2_range);
Z3 = Z3a(1-low:end-high,x3_range);                      Z4 = Z4a(1-low:end-high,x4_range);
x_mm = x_length1 + 3*x_length2; %x_length2*3 when using non-shielded
x_size = length(x1_range)+length(x2_range)+length(x3_range)+length(x4_range);
x_full = start_mm:(x_mm-start_mm)/(x_size-1):x_mm;      y_full = y_co1(y1_range);
%% ROW SIZE CHECK BETWEEN CAMERA IMAGES
[Z1,Z2,Z3,Z4] = RowSizeCheck(Z1,Z2,Z3,Z4);
%% COMBINING MARTICES & SMOOTHING IMAGE BORDERS
Z_full = cat(2,Z1,Z2,Z3,Z4);        x_full = x_full(1:size(Z_full,2));   y_full = y_full(1:size(Z_full),1);
%  IMAGE SMOOTHING VIA BOUNDARY BLUR(AVERAGING GROUPS OF COLUMNS AT IMAGE EDGES)
if smoothing == 1
   for i=1:size(Z_full,1)
     avg1 = mean(diff([Z1(i,end-3) Z1(i,end-2) Z1(i,end-1) Z1(i,end) Z2(i,1) Z2(i,2) Z2(i,3) Z2(i,4)]));
        temp1(i,1:8)=[Z1(i,end-3) Z1(i,end-3)+avg1 Z1(i,end-3)+avg1*2 Z1(i,end-3)+avg1*3 Z1(i,end-3)+avg1*4 Z1(i,end-3)+avg1*5 Z1(i,end-3)+avg1*6 Z1(i,end-3)+avg1*7];
        
     avg2 = mean(diff([Z2(i,end-3) Z2(i,end-2) Z2(i,end-1) Z2(i,end) Z3(i,1) Z3(i,2) Z3(i,3) Z3(i,4)]));
        temp2(i,1:8)=[Z2(i,end-3) Z2(i,end-3)+avg2 Z2(i,end-3)+avg2*2 Z2(i,end-3)+avg2*3 Z2(i,end-3)+avg2*4 Z2(i,end-3)+avg2*5 Z2(i,end-3)+avg2*6 Z2(i,end-3)+avg2*7];
        
     avg3 = mean(diff([Z3(i,end-3) Z3(i,end-2) Z3(i,end-1) Z3(i,end) Z4(i,1) Z4(i,2) Z4(i,3) Z4(i,4)]));
        temp3(i,1:8)=[Z3(i,end-3) Z3(i,end-3)+avg3 Z3(i,end-3)+avg3*2 Z3(i,end-3)+avg3*3 Z3(i,end-3)+avg3*4 Z3(i,end-3)+avg3*5 Z3(i,end-3)+avg3*6 Z3(i,end-3)+avg3*7];
   end
  Z_full(:,size(Z1,2)-3:size(Z1,2)+4)=temp1;
  Z_full(:,size(Z1,2)+size(Z2,2)-3:size(Z1,2)+size(Z2,2))=temp2(:,1:4);
  Z_full(:,size(Z1,2)+size(Z2,2)+size(Z3,2)-3:size(Z1,2)+size(Z2,2)+size(Z3,2)+4)=temp3;
else
end
size(x_full);  size(y_full); size(Z_full)
%%  PLOT RESULT IMAGE
figure(1000); pcolor(x_full/Deq,y_full/Deq,Z_full); axis equal; shading interp; colormap(jet); colorbar;
ylim([min(y_full/Deq) max(y_full/Deq)]);   set(gcf,'Position',[33 244 1280 443]);
switch result_typ
  case 1 
     caxis([0 max(max(Z_full))]*1.2);
  case 2
     caxis([-max(max(Z_full)) max(max(Z_full))]);
  case 3
     caxis([0 0.04]);   % caxis([0 0.01]);   
     %caxis([0 6000]);
  case 4
     caxis([0 5]);
  case 5
     caxis([0 5]);
  case 6
     caxis([-2 2]);
end
%  PLOT CENTERLINE
Plot_Center(x_full,Z_full,result_typ,Deq);
%%  REARRANGE & SAVE AXES AND DATA TO .TXT FILE
if strcmp(data_save,'y')
   saveX = repmat(x_full',size(y_full,1),1);          N = size(x_full,2);
   saveY = y_full(repmat(1:size(y_full,1),N,1),:);    saveZ = reshape(Z_full.',[],1);
   size(saveX);  size(saveY);  size(saveZ);
   saveAll = [saveX saveY saveZ];  size(saveAll);     
   save([drive_in condition{1}(1:7) '_' condition{1}(9:14) '_' nozzle '_' type '.dat'],'saveAll','-ascii');
   screensize = get( groot, 'Screensize' );
%  SAVE RESULT IMAGE
   saveas(figure(1000),[drive_in condition{1}(1:7) '_' condition{1}(9:14) '_' nozzle '_' type '.jpg'])
%  SAVING STICH VALUES
   Vals = [];
   Vals(:,1) = [x1_del1,x1_del2,y1_del1];      %FIRST IMAGE CUTS
   Vals(:,2) = [x2_del1,x2_del2,y2_del1];      %SECOND IMAGE CUTS
   Vals(:,3) = [x3_del1,x3_del2,y3_del1];      %THIRD IMAGE CUTS
   Vals(:,4) = [x4_del1,x4_del2,y4_del1];      %FOURTH IMAGE CUTS
   Vals(:,5) = start_mm - (x1_del1*0.63);      %DIST. FROM EXIT TO FIRST IMAGE
   Vals(:,6) = [val0,val1,val2];               %MAGNITUDE CORRECTION
   if result_typ == 1
      file_out_data = ['End_Values_AvgVx_' condition{1}(9:14) '_' condition{1}(1:7)];
   elseif result_typ == 2
      file_out_data = ['End_Values_AvgVy_' condition{1}(9:14) '_' condition{1}(1:7)];
   else
      file_out_data = ['End_Values_TurbKE_' condition{1}(9:14) '_' condition{1}(1:7)];
   end
   save([drive_in file_out_data],'Vals');
end
%% CENTERLINE VELOCITY
function Plot_Center(x_full,Z_full,result_typ,De)
x = 1;
if result_typ == 1
   while Z_full(x,20) < 250
         x = x+1; 
   end
   jet_top = x;    
   while Z_full(x,20) > 250
         x = x+1;
   end
   jet_bot = x-1;
   centr = round((jet_top+jet_bot)/2);
   figure; subplot(211);plot(x_full/De,Z_full(centr,:),'LineWidth',1.2); xlim([0 max(x_full/De)]); grid on; 
   subplot(212); plot(x_full,Z_full(centr,:),'LineWidth',1.2); xlim([0 max(x_full)]); 
   set(gcf,'Position',[350 120 1200 780]); grid on;
end
end
%% MATRIX SIZE FUNCTION
function [row,col,matrix_size,clim,txtFormat,lineSkip] = Matrix_Size(drive_in,config,result_typ,position)
if strcmp(config(end),'D')
   if result_typ == 1 
      type = ['AvgVx_' num2str(position)];            clim = [0 550];
   elseif result_typ == 2
      type = ['AvgVy_' num2str(position)];            clim = [-30 30];
   else
      type = ['TurbKineticE_' num2str(position)];     clim = [0 9000];
   end
   fid_mat = fopen([drive_in type '.DAT']);   mat_siz = textscan(fid_mat,'%s%s%s%s%s',5); 
   col = mat_siz{4}(4); row = mat_siz{5}(4);  clear fid_vel;  txtFormat = '%n%n%n%n%n'; 
   col = str2double(col{1}(3:5)); row = str2double(row{1}(3:5)); lineSkip = 3;
else
   if result_typ == 1
      type = ['AvgVx_' num2str(position)];         clim = [0 550];
   elseif result_typ == 2
      type = ['AvgVx_' num2str(position)];         clim = [-30 30];
   else
      type = ['AvgVx_' num2str(position)];         clim = [0 0.04];
   end
   fid_mat = fopen([drive_in 'Velocity Filters_D.DAT']);  mat_siz = textscan(fid_mat,'%d%d','headerlines',3,'Delimiter',',');
   row = mat_siz{1};  col = mat_siz{2};   txtFormat = '%n%n%n';  lineSkip = 0;
end
matrix_size = row*col; disp([newline '>> VECTOR FIELD RESOLUTION(RxC) : ',num2str(row),' X ',num2str(col)]);
end
%% CHECK ROW SIZE BETWEEN CAMERA IMAGES
function [Z1,Z2,Z3,Z4] = RowSizeCheck(Z1,Z2,Z3,Z4)
%  CHECK ROW SIZE BETWEEN CAMERA 1 & 2
if size(Z1,1) < size(Z2,1)
   Z2 = Z2(1:size(Z1,1),:);
elseif size(Z1,1) > size(Z2,1)
   Z1 = Z1(1:size(Z2,1),:);
end
%  CHECK ROW SIZE BETWEEN CAMERA 2 & 3
if size(Z2,1) < size(Z3,1)
   Z3 = Z3(1:size(Z2,1),:);
elseif size(Z2,1) > size(Z3,1)
   Z2 = Z2(1:size(Z3,1),:);
end
%  CHECK ROW SIZE BETWEEN CAMERA 3 & 4
if size(Z3,1) < size(Z4,1)
   Z4 = Z4(1:size(Z3,1),:);
elseif size(Z3,1) > size(Z4,1)
   Z3 = Z3(1:size(Z4,1),:);
end
end
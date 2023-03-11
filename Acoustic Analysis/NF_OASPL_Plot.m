%% PLOTTING NEARFIELD OASPL
% CONDITIONS - CODE PLOTS MULTIPLE PLOTS DEPENDING ON THE SPECIFIED
% CONDITIONS
tic; set(0,'defaultfigurecolor',[1 1 1]);
clc; clearvars; close all;  root = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\Near_Field_Data\';
tests = {'NPR_2p5_TR_1p0', 'NPR_2p9_TR_1p0','NPR_3p0_TR_1p0', 'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_4p5_TR_1p0','NPR_5p0_TR_1p0', ...
         'NPR_2p5_TR_2p8','NPR_3p0_TR_2p8','NPR_3p6_TR_2p8','NPR_4p0_TR_2p8','NPR_4p5_TR_2p8','NPR_5p0_TR_2p8'};

condition = [tests(1)];        nozzle  =   {'Major'};   jet_typ = 'S';

switch jet_typ
  case 'S' 
    drive = [root 'Single_Rectan\Results\'];   De = 0.813; lim_x = [-0.51 44]; lim_y = [0 22]; Hold = zeros(36,64);
  case  'TR'
    drive = [root 'Twin_Jet\Results\'];     De = 0.727; lim_x = [-0.7 45];  lim_y = [0 25]; Hold = zeros(36,64);
  case  'TR2'
    drive = [root 'Twin_Jet_V2\Results\'];  De = 0.727; lim_x = [-0.7 30];  lim_y = [-1.95 1.95]; Hold = zeros(14,64);
end
% VARIABLES 
c = 0;      smooth = 'y';   textsize = 13;                
% PLOT PARAMATERS
mic_rem = [3,4];                                         
%% JASA SPECIFIC
% cl = colorbar;  title(cl,'\bfOASPL(dB)'); xlabel('\bfX/D_e'); ylabel('\bfY/D_e'), ylabel('\bfZ/D_e');
% caxis([120 145]); set(gcf,'Position',[325 315 675 295]); xlim([-5 44]);
%%
for m = 1:length(nozzle)
 for n = 1:size(condition,1)*size(condition,2)
   drive_in  = [drive nozzle{m} '\' condition{n}(9:14) '\' condition{n}(1:7) '\'];     c = 0;
   drive_out = [drive 'Peak_plots\' nozzle{m} '\' condition{n}(9:14) '\' condition{n}(1:7) '\']; 
%  OASPL MATRIX FORMULATION
   if strcmp(jet_typ,'T2')
      [Hold,~] = NF_MxMatch_V2(condition,drive_in,'OASPL',[],[],[],[]);
   else
      [Hold,~] = NF_MxMatch(condition,drive_in,'OASPL',[],[],[],[]);
   end
%  DEFINING AXES
   [X,Y,Y_off] = NF_Axis_Def(jet_typ,nozzle{m});   
   if strcmp(jet_typ,'S') 
      r = 49:64;    Y(:,r) = Y(:,r) - 0.3875;
   end
     data = flipud(Hold);  clear Hold;  
%  REMOVING THE BAD MICS
   for ctr = mic_rem 
     if ctr ~= 0
        rez_data = data(:,[1:ctr-1 ctr+1:ctr+15 ctr+17:ctr+31 ctr+33:ctr+47 ctr+49:64]);
        rez_X = X(:,[1:ctr-1 ctr+1:ctr+15 ctr+17:ctr+31 ctr+33:ctr+47 ctr+49:64]);
        rez_Y = Y(:,[1:ctr-1 ctr+1:ctr+15 ctr+17:ctr+31 ctr+33:ctr+47 ctr+49:64]);
     else
       rez_X = X;   rez_Y = Y;  rez_data = data; break;
     end
   end
% data(1,1:64) = 0;%making the bottom row equal to zero to reduce saturation
   if strcmp(condition{n}(9:14),'TR_1p0')                  %dB RANGE FOR CONTOUR PLOT
      lim_u   = 145;        lim_d = 110;
   else
      lim_u   = 155;        lim_d = 130;
   end
   if strcmp(smooth,'y')
      sc_1 = 0.65;                                         %SCALING FACTOR FOR CONTOUR SMOOTHING
      rez_data1 = imresize(rez_data,sc_1,'bilinear');
      rez_X1 = (imresize(rez_X,sc_1,'bilinear'))/De;
      rez_Y1 = (imresize(rez_Y,sc_1,'bilinear'))/De;
      sc_2 = size(data);
      rez_data2 = imresize(rez_data1,sc_2,'bilinear');
      rez_X2 = imresize(rez_X1,sc_2,'bilinear')+(rez_X(1)/De-rez_X1(1));
      rez_Y2 = imresize(rez_Y1,sc_2,'bilinear'); 
      if strcmp(jet_typ,'T')
         rez_Y2 = rez_Y2 + (Y_off - rez_Y2(1));
      elseif strcmp(jet_typ,'S')
         rez_Y2 = rez_Y2 + (Y_off - rez_Y2(1,2));
      end
   else
      rez_data2 = rez_data; rez_X2 = rez_X/De; rez_Y2 = rez_Y/De + Y_off;
   end
   if (length(condition)>2)    
       subplot(2,2,n)
   else
       subplot(1,size(condition,2),n)
   end
   if strcmp(nozzle{m},'Minor')
      yname = 'Y/D_e';
   else
      yname = 'Z/D_e';
   end
   figure(1); %imshow(img_maj); hold on;
   contourf(rez_X2,rez_Y2,rez_data2,10),axis equal, %axis([-0.38 max(max(X))/D 0 max(max(Y))/D+1.598]);
   oaspl = colorbar;   title(oaspl, 'OASPL(dB)');   colormap jet ;   caxis([lim_d lim_u]);
   xlim(lim_x);  ylim(lim_y); box off;     xlabel('X/D_e'); ylabel(yname);
   grid on;    set(gca, 'layer', 'top');    hold on;       set(gca,'FontSize',textsize,'FontName','times')
%    Hold = zeros(64,36); %RESETTING THE HOLDING MATRIX TO ZERO FOR NEXT LOOP%%% 
    %% OBTAIN DATA FROM CONTOUR PLOT
%   [c,h] = contourf(rez_X2,rez_Y2,rez_data2,10);
%   h.LevelList = round(h.LevelList,0); 
%   v = [120:3:130,130:2:140];
%   clabel(c,h,'FontSize',11,'FontWeight','bold');
 end
end 
%% FAR FIELD ANGLE MARKERS
commandwindow;  prompt4 = [newline '>> ADD FARFIELD ANGLES(y/n) - '];angl_on = input(prompt4,'s');
if strcmp(angl_on,'y')
   NF_Angle_Plot_V2(1,jet_typ);  xlim([-5 46]);     ylim(lim_y);
end
if n<2
   set(gcf,'Position',[750 450 650 380]);
else
   set(gcf,'Position',[420 493 1244 375]);
end
commandwindow; %prompt = [newline '>> SAVE IMAGE?(y/n) - '];  data_save = input(prompt,'s');
%% DATA SAVE
% prompt = [newline '>> SAVE PLOT DATA(y/n) - '];        data_save = input(prompt,'s');   
data_save = 'n';
if strcmp(data_save,'y')
   cd(drive_out);    
   saveas(figure(1),['OASPL_AC_' nozzle{n} '_' condition{n}(9:14) '.png']); 
end
drive_code = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\Jet_Analysis\Near_Field_Codes\';
cd(drive_code);
%% EXTRACTING DATA FROM THE RESIZED CONTOUR PLOT
%      Directvty = data;      %SAVING DATA FOR EACH CONDITION
%      file_out_data = ['Directvty_' freqs{fr} '_' condition{n}(9:14) '_' condition{n}(1:7)];
%      save([drive_out file_out_data],'Directvty');

H = findobj(figure(1),'type','contour');
x_data = get(H,'Xdata');
y_data = get(H,'Ydata');
oaspl_data = get(H,'Zdata');
figure(2);contourf(x_data,y_data,oaspl_data,10);axis equal; grid on;
oaspl1 = colorbar;   title(oaspl1, 'OASPL(dB)');   colormap jet ;   caxis([lim_d lim_u]);
xlim(lim_x);  ylim(lim_y); box off;     xlabel('X/D_e'); ylabel(yname); set(gca,'FontSize',textsize,'FontName','times')
toc

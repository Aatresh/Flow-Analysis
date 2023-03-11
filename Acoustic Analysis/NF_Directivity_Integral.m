%% CODE PLOTS THE OASPL FOR THE SPECIFIED FREQUENCY RANGE THROUGH NUMERICAL
%  INTEGRATION OVER SPECIFIED FREQUENCY RANGE FOR A GIVEN FLOW CONDITION
%  AVOIDING MICS 3,4 AND OPTIONAL WAVELET SMOOTHING

clc; clearvars; close all;  set(0,'defaultfigurecolor',[1 1 1]);
root = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';
cd([root 'Jet_Analysis\Global_Functions\']);       addpath([root 'Jet_Analysis\Near_Field_Codes\']);
tests = {'NPR_2p5_TR_1p0','NPR_2p9_TR_1p0','NPR_3p0_TR_1p0','NPR_3p6_TR_1p0','NPR_4p0_TR_1p0','NPR_5p0_TR_1p0', ...
         'NPR_2p5_TR_2p0','NPR_3p0_TR_2p8','NPR_3p6_TR_2p8','NPR_4p0_TR_2p8','NPR_5p0_TR_2p8'};

condition = [tests(3)];     nozzle = {'Minor'};     jet_typ = 'S';

switch jet_typ
  case 'S' 
    drive = [root 'Near_Field_Data\Single_Rectan\Results\'];    Deq = 0.813; lim_x = [-0.51 44]; lim_y = [0 22]; rowNo = 36; colNo = 64;
    prompt = [newline '>> PLOT OASPL?(y/n) - '];  oaspl = input(prompt,'s');   scl = 'holdr';
  case  'T'
    drive = [root 'Near_Field_Data\Twin_Rectan\Results\'];      Deq = 0.727; lim_x = [-0.7 45];  lim_y = [0 25]; rowNo = 36; colNo = 64;
    prompt = [newline '>> PLOT TRUE OASPL OR SCALE TO SINGLE JET(t/s)'];  scl = input(prompt,'s');
  case  'T2'
    drive = [root 'Near_Field_Data\Twin_Rectan_V2\Results\'];   Deq = 0.727; lim_x = [-0.7 40];  lim_y = [-1.95 1.95]; rowNo = 14; colNo = 64;
end
%% JASA SPECIFIC
% cl = colorbar;  title(cl,'\bfSPL(dB)'); xlabel('\bfX/D_e'); ylabel('\bfY/D_e'), ylabel('\bfZ/D_e');
% caxis([90 140]); set(gcf,'Position',[325 315 675 297]); xlim([-5 44]);
%%  CONTROL PARAMETERS
offset    = 1;                                     % DOMAIN OFFSET AT 30De 
smooth    = 'y';    sc_1 = 0.6;                    % SCALING FACTOR FOR CONTOUR SMOOTHING  
wv_smooth = 'n';                                   % WAVELET SMOOTHING
mic_rem   = [3,4];                                 % REMOVING DATA FROM MIC 3 or 4
lim_u = 140;        lim_d = 90;                    % LIMITS 
% VARIABLES
a = 0; b = 0; y = 0; u = 0; ct = 0;   chc = 'n';  textsize = 13;   
Hold = zeros(36,64);        %  HOLDING MATRIX FOR SPLs
if strcmp(oaspl,'y')
   fRange = 300:50:100000; fRange = fRange/50;  bar_title = 'OASPL(dB)';
else
   prompt  = [newline '>> PLOT RANGE OR SINGLE FREQUENCY(r/s) - ']; chc = input(prompt,'s');
   if strcmp(chc,'s')
      prompt = [newline '>> ENTER FREQUENCY - '];      fRange = input(prompt);
   else
      prompt = [newline '>> ENTER STARTING FREQUENCY(Hz) - '];    lim1 = input(prompt);
      prompt = [newline '>> ENTER ENDING FREQUENCY(Hz) - '];      lim2 = input(prompt);
      fRange = lim1:50:lim2;                       %MATRIX CONTAINING THE FREQUENCY RANGE TO COMPUTE OASPL

   end;     bar_title = 'SPL(dB)';      fRange = fRange/50;     oaspl = 'n';
end
Matrix = zeros(rowNo,colNo,length(fRange));        %MATRIX HOLDING FREQUENCY BANDS
tic
for n = 1:length(condition)
  drive_in = [drive nozzle{1} '\' condition{n}(9:14) '\' condition{n}(1:7) '\']; 
  drive_out = [drive 'Peak_plots\'  nozzle{1} '\' condition{n}(9:14) '\' condition{n}(1:7) '\'];
  [Mj,Uj,NPR,NTR] = GF_Velocity(condition{n});
%  MATRIX FORMULATION
  if strcmp(jet_typ,'T2')
     [Matrix,freq] = NF_MxMatch_V2(condition,drive_in,'SPL',fRange,jet_typ,oaspl,wv_smooth);
  else
     [Matrix,freq] = NF_MxMatch(condition,drive_in,'SPL',fRange,jet_typ,oaspl,scl,wv_smooth);
  end
%  CALCULATING THE PRESURE FROM SPL VALUES
  Hold1 = 10.^(Matrix/10);        Hold2=(Hold1.*(20e-6)^2./50)/2;
  df = freq(fRange)';             data = zeros(size(Hold1,1),size(Hold1,2));
  if length(fRange) == 1                           %NUMERICAL INTEGRATION TO COMPUTE OASPL
     data = Matrix;
  else
   wB = waitbar(0,'Computing Frequency Integral SPL');
     for a = 1:size(Hold1,1)
       for b = 1:size(Hold1,2)
         data(a,b)=trapz(df,Hold2(a,b,:));  
         data(a,b)=10*log10(data(a,b)/(20e-6)^2);         waitbar(a/size(Hold1,1),wB);
       end
     end
  end
%  DEFINING AXES
  [X,Y,Y_off] = NF_Axis_Def(jet_typ,nozzle{n}); 
%  ADDING OFSET TO ACCOUNT FOR TRAVERSE MOVEMENT
  if strcmp(jet_typ,'S')
     rowNo = 49:64;    Y(:,rowNo) = Y(:,rowNo) - 0.3875;
  end
%  REMOVING SELECTED MICS
  for ctr = mic_rem 
     if ctr ~= 0
        rez_data = data(:,[1:ctr-1 ctr+1:ctr+15 ctr+17:ctr+31 ctr+33:ctr+47 ctr+49:64]);
        rez_X = X(:,[1:ctr-1 ctr+1:ctr+15 ctr+17:ctr+31 ctr+33:ctr+47 ctr+49:64]);
        rez_Y = Y(:,[1:ctr-1 ctr+1:ctr+15 ctr+17:ctr+31 ctr+33:ctr+47 ctr+49:64]);
     else
       rez_X = X;   rez_Y = Y ;  rez_data = data; break;
     end
  end
  if strcmp(nozzle{n},'Minor') 
     yname = 'Y/D_e';
  else
     yname = 'Z/D_e';
  end
%  SETTING SUBPLOT CONDITIONS   
  if length(condition)>2 
     subplot(2,2,n)
  else
     subplot(1,size(condition,2),n)
  end
  hold on
  if strcmp(smooth,'y')
     rez_data1 = imresize(rez_data,sc_1,'bilinear');
     rez_X1 = (imresize(rez_X,sc_1,'bilinear'))/Deq;
     rez_Y1 = (imresize(rez_Y,sc_1,'bilinear'))/Deq;
     sc_2 = size(data);
     rez_data2 = imresize(rez_data1,sc_2,'bilinear');
     rez_X2 = imresize(rez_X1,sc_2,'bilinear')+(rez_X(1)/Deq-rez_X1(1));
     rez_Y2 = imresize(rez_Y1,sc_2,'bilinear');
     if strcmp(jet_typ,'T')
         rez_Y2 = rez_Y2 + (Y_off - rez_Y2(1));
     elseif strcmp(jet_typ,'S')
         rez_Y2 = rez_Y2 + (Y_off - rez_Y2(1,2));
      end
  else
     rez_data2 = rez_data; rez_X2 = rez_X/Deq; rez_Y2 = rez_Y/Deq + Y_off;
  end; rez_Y2 = flipud(rez_Y2);
  if strcmp(jet_typ,'T') && strcmp(oaspl,'y')
     disp([newline '>> OSASPL of Twin Jet scaled to Single Jet']);
  else
     disp([newline '>> OASPL - Single Jet']);
  end
  figure(1)
  contourf(rez_X2,rez_Y2,rez_data2,10);   hold on ;    
  colormap jet;     caxis([lim_d lim_u]);  oaspl = colorbar;      title(oaspl,bar_title,'FontSize',13);
  axis equal; xlim(lim_x); ylim(lim_y);   %xticks(-0.38:9:42);
  xlabel('X/D_e');  ylabel(yname); set(gca,'FontSize',textsize,'FontName','times')
  grid on;          set(gca, 'layer', 'top');
end
toc
%% FARFEILD ANGLE MARKERS
commandwindow;  prompt4 = [newline '>> ADD FARFIELD ANGLES(y/n) - '];angl_on = input(prompt4,'s');
if strcmp(angl_on,'y')
   NF_Angle_Plot_V2(1,jet_typ);figure(1);   xlim([-5 46]);     ylim(lim_y);
end
if n<2 && strcmp(jet_typ,'T2')~=1
   set(gcf,'Position',[750 450 650 380]);
elseif n >2 && (strcmp(jet_typ,'T') || strcmp(jet_typ,'T')) == 1
   set(gcf,'Position',[420 493 1244 375]);
elseif n<2 && strcmp(jet_typ,'T2')
   axis normal; set(gcf,'Position',[75 520 925 285]);
end
commandwindow; %prompt = [newline '>> SAVE IMAGE?(y/n) - '];  data_save = input(prompt,'s');
%% DISPLAY St RANGE
if strcmp(chc,'r') 
   St1 = lim1*(Deq*0.0254)/Uj;    St2 = lim2*(Deq*0.0254)/Uj;
   disp([newline,'>> Strouhal Range - ',num2str(St1),' - ',num2str(St2)]);
elseif strcmp(chc,'s')
   St = freq(fRange)*(Deq*0.0254)/Uj; disp([newline,'>> Strouhal number - ',num2str(St)]);
end
%% IMAGE SAVE
data_save = 'n';
if strcmp(data_save,'y')
% %    drive_out = 'X:\OneDrive - University of Cincinnati\Working_Directory\Publicaitons & Conferences\Journal Drafts\Screech Characterization\';
%    drive_code = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\Jet_Analysis\Near_Field_Codes\';
%    cd(drive_out);    
%    if strcmp(chc,'s')
%       saveas(figure(1),['NF_' nozzle{n} '_FR_' num2str(freq(frange)) '.png']);
%    elseif strcmp(chc,'n') && strcmp(scl,'y')
%       saveas(figure(1),['OASPL_' nozzle{n} '_' condition{n}(9:14) '.png']); 
%    else
%       saveas(figure(1),['NF_' nozzle{n} '_' condition{n}(9:14)  ' FR_' num2str(lim1) '_' num2str(lim2) '.png']);
%    end
%    cd(drive_code);
%% DATA SAVE FOR SCHLIEREN CALIBRATION
%  OASPL IS PLOTTED WITHIN THE CHAMBER RANGE - 300 TO 100,000 Hz
drive_out = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\Schlieren\Single_Jet\OASPL\';
OASPL = ['OAPSL_CR_' nozzle{n} '_' condition{n}(9:14) '_' condition{n}(1:7)];
save([drive_out OASPL],'rez_data2');  save([drive_out 'X'],'X');  save([drive_out 'Y'],'Y');
%      Directvty = data;      %SAVING DATA FOR EACH CONDITION
%      file_out_data = ['Directvty_' freqs{fr} '_' condition{n}(9:14) '_' condition{n}(1:7)];
%      save([drive_out file_out_data],'Directvty');
end
beep


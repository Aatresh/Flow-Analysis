%% PLOTS THE FAR FIELD SPECTRUM AS ANGLE VS STROUHAL NUMBER
% ALWAYS PLOT IN PAIRS FOR TEXT TO APPEAR
clc;    clearvars;  set(0,'defaultfigurecolor',[1 1 1]);    code = 'FF_Plot';
root = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';  cd([root 'Jet_Analysis\Global_Functions\']);
addpath([root 'Jet_Analysis\Far_Field_Acoustic_Codes\']);   load('schFaMap.mat');  
tests = {'NPR_2p5_TR_1p0','NPR_2p6_TR_1p0','NPR_2p9_TR_1p0','NPR_3p0_TR_1p0','NPR_3p6_TR_1p0','NPR_4p0_TR_1p0','NPR_4p5_TR_1p0','NPR_5p0_TR_1p0', ...
         'NPR_2p5_TR_2p0','NPR_3p0_TR_2p0','NPR_3p6_TR_2p0','NPR_4p0_TR_2p0','NPR_4p5_TR_2p0','NPR_5p0_TR_2p0'};

condition = [tests(4)];   config = 'TS';   Nozzle = {'Major';'Minor'};


prompt = [newline '>> Show Titles?(y/n) - '];       shwTitle   = input(prompt,'s');
prompt = [newline '>> Show Noise Types?(y/n) - '];  shwNoisTyp = input(prompt,'s');
k = 1;   textsize = 13;

for m = 1:length(Nozzle)
  for n = 1:length(condition)
    [Mj,Uj,NPR,TR] = GF_Velocity(condition{n});
%   Input & Output Drive Selection 
    [OutputStruct] = GF_DriveSelect(config,Nozzle{m},code);   nozzle = OutputStruct.nozzle;
%   Input Drive;                                              Far field angles                                           
    driveIn = OutputStruct.in_root;                           ffAngles = fliplr(OutputStruct.bckgrnd_root);
%   Far field mic distance;                                   Equivalent Diameter for St scaling
    ffDistnc = OutputStruct.dt.ffDistnc;                      Deq = OutputStruct.dt.Deq; 
    fileName = ['FFT_Nick_nofilter_N' condition{n}(9:14) '_' condition{n}(1:7)];
%   Loading SPL File;                                                                       
    val = load([driveIn 'FFT_Nick_nofilter_N' condition{n}(9:14) '_' condition{n}(1:7)]); 
%   SPL Spectra;                         Frequency Range          
    splSpectra = val.PHI_yy_amp_avg;     splFreq = val.freq;     
%   Nyquist Cut off;                     Subplot counter          
    nyqCut = size(splSpectra,2)/2;       subLoc = 1;              
%   Converting frequency to Strouhal 
    splFreq = splFreq*(Deq/Uj);   
%   Plotting pcolor 
    figure(1)
    if length(condition) > 1 || length(Nozzle)>1
       subplot(length(Nozzle),length(condition),k); 
    end
    pcolor(splFreq(1:nyqCut),ffAngles,splSpectra(:,1:nyqCut));    shading interp;    colormap(schFaMap);
    c = colorbar;   hold on;    yticks([45 70 90 110 132 152]);   xticks([0.5 1 1.5 2 2.5 3 3.5 4 4.5]); 
    title(c,'SPL(dB)','Interpreter','latex');                     set(gca,'FontSize',textsize);
    ylabel('$Angles$','Interpreter','latex');                     xlabel('$St$','Interpreter','latex');
    aX = gca;       aX.TickLabelInterpreter = 'latex';            c.TickLabelInterpreter = 'latex';
    if TR < 2
       xlim([0 3]);     caxis([70 110]);
    else
       xlim([0 3]);     caxis([65 120]);
    end
%   Display condition on top of contour
    if strcmp(shwTitle,'y')
      condName = ['NPR - ',num2str(NPR)];   title(condName,'Interpreter','latex');      
      cL = legend(Nozzle{m});               cL.Interpreter = 'latex';        cL.EdgeColor = [1 1 1];
    end;                                    k = k+1;
    if strcmp(shwNoisTyp,'y')
       if NPR == 2.9  
          text(0.2,105, 'H1','Color','k','FontSize',12,'FontWeight','bold','Interpreter','latex');
          text(0.58,105,'H2','Color','k','FontSize',12,'FontWeight','bold','Interpreter','latex');
          text(1,60,    'H3','Color','k','FontSize',12,'FontWeight','bold','Interpreter','latex');
       elseif NPR == 3 
          text(0.13,105,'H1','Color','k','FontSize',12,'FontWeight','bold','Interpreter','latex');
          text(0.55,105,'H2','Color','k','FontSize',12,'FontWeight','bold','Interpreter','latex');
          text(0.92,60, 'H3','Color','k','FontSize',12,'FontWeight','bold','Interpreter','latex');
       end
       text(0.57,144,'\bf LST','Color','k','FontSize',12,'Interpreter','latex');
       text(2,120,  '\bf BBSN','Color','k','FontSize',12,'Interpreter','latex');
    end
  end
end
%  Change plot size & Position based on number of conditions being plot
if length(nozzle) > 1 && length(condition) == 1
   set(gcf,'Position',[140 90 570 740]);
elseif length(condition)>1 && length(nozzle) > 1
   set(gcf,'Position',[123 142 1043 777]);   
elseif length(condition) > 1 && length(nozzle) == 1
   set(gcf,'Position',[130 450 1260 305]);
end
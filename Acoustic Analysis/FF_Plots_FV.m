%% PLOTS THE FAR-FIELD ACOUSTIC PLOTS FOR SINGLE JET or TWIN JET
% All results are scaled to a distance of 44.715 inches from nozzle centerline
% Plot combinations: Multi NPR    - Single Nozzle 
%                  : Multi Nozzle - Single NPR 
% Plot markers based on config: C - Circle; S - Diamond; TR - Pyramid; TS - Square             
clc;   clearvars;   set(0,'defaultfigurecolor',[1 1 1]);   
root = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';   code = 'FF_Plot';
cd([root 'Jet_Analysis\Global_Functions\']);       addpath([root 'Jet_Analysis\Far_Field_Acoustic_codes\']);
% Test conditions
tests = {'NPR_2p5_TR_1p0', 'NPR_2p6_TR_1p0', 'NPR_2p9_TR_1p0', 'NPR_3p0_TR_1p0', 'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_4p5_TR_1p0','NPR_5p0_TR_1p0',...
         'NPR_2p5_TR_2p0', 'NPR_2p9_TR_2p0', 'NPR_3p0_TR_2p0', 'NPR_3p6_TR_2p0', 'NPR_4p0_TR_2p0', 'NPR_4p5_TR_2p0', 'NPR_5p0_TR_2p0'};
%          'NPR_2p6_TR_1p0', 'NPR_2p7_TR_1p0','NPR_3p9_TR_1p0',};

condition = [tests(1)];    config = 'S';    Nozzle = {'Minor'};

% Sequence of angles to plot                  Pressure Statistics
angSeq = [1,13,16];                              statNo = 1;         
% Choice of pressure statistics to plot;      Font Size for plots
[statName,yVals] = statChoice(statNo);        textsize = 13;
% Write data to text file                     Frequency prompt      
write = 'of';                                 prompt = [newline '->> St or Hz(1/2) - ']; 
St = input(prompt);                           colorCtr = 1;  % Counter to change color
% Marker Set                                  Config Set
mrkrSet = {'o';'d';'^';'s'};                  configSet = {'C';'S';'TR';'TS'};
inDx = strcmp(config,configSet);              plotMrkr = cell2mat(mrkrSet(inDx));
% Plot Colormap (colr13 in GF_ColorTester)
colr = [0.9 0 0;0.9 0.6 0.1;0.8 0.5 0.8;0 0 0.9;0.6 0.6 0.6;0.9 0.2 0.5;0.2 0.2 0.2;0.5 0.9 0.6;0.2 0.6 0.9;];
% Default Color Set
% colr = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840;0.3 0.3 0.3];
% Inputs for color selector
selectrInp.colr  = colr;                      selectrInp.Nozzle    = Nozzle;
selectrInp.tests = tests;                     selectrInp.condtnCtr = length(condition);
%---------------------------------- Computing & Plotting Loop --------------------------------------
for m = 1:length(Nozzle)
 for n = 1:length(condition)
   [Mj,Uj,NPR,TR] = GF_Velocity(condition{n});
%  Input & Output Drive Selection 
   [OutputStruct] = GF_DriveSelect(config,Nozzle{m},code);   nozzle = OutputStruct.nozzle;
%  Legend tag, Plot title & Marker Color based on test conditions
   selectrInp.nozzle = nozzle;                selectrInp.condition = condition{n};
   selectrInp.NPR    = NPR;                   selectrInp.colorCtr  = colorCtr;
   selectrInp.yVals  = yVals;   
   [condName,plotName,colorCtr,mrkrFaceColr,lineColr] = colorSelector(selectrInp);
%  Input Drive;                                              Far field angles                                           
   driveIn = OutputStruct.in_root;                           ffAngles = fliplr(OutputStruct.bckgrnd_root);
%  Far field mic distance;                                   Equivalent Diameter for St scaling
   ffDistnc = OutputStruct.dt.ffDistnc;                      Deq = OutputStruct.dt.Deq; 
%---------------------------------- Plotting Pressure Statictics -----------------------------------
%  Loading Statistics file
   val = load([driveIn statName 'N' condition{n}(9:14) '_' condition{n}(1:7)]);   presStat = struct2array(val);
%  Setting limits for stat value;                            YAxis Label                                             
   limY  = yVals.limY;                                       yName = yVals.yName;                                    
%  Plotting Pressure Statistics
   figure(1);   plot(ffAngles,presStat(:,1),['-.',plotMrkr],'color',lineColr,'LineWidth',1.2,'MarkerSize',5,'MarkerFaceColor',mrkrFaceColr,'DisplayName',condName);
   xlim([40 160]);   xlabel('$Azimuthal \thinspace Angle(\psi)$','Interpreter','latex');   xticks([40 60 80 100 120 140 160]);  
   ylim(limY);       ylabel(['$',yName,'$'],'Interpreter','latex');                        aX = gca;
   title(['$',plotName,'$'],'Interpreter','latex');   box off;   grid on;       set(gca,'FontSize',textsize);                               
   aX.TickLabelInterpreter = 'latex';   hold on;   cL = legend;  cL.EdgeColor = [1 1 1];   cL.Interpreter = 'latex';                     
%---------------------------------- Plotting Frequency Spectra -------------------------------------
%  Loading SPL File;                                                                       
   val = load([driveIn 'FFT_Nick_nofilter_N' condition{n}(9:14) '_' condition{n}(1:7)]); 
%  SPL Spectra;                         Frequency Range          
   splSpectra = val.PHI_yy_amp_avg;     splFreq = val.freq;     
%  Nyquist Cut off;                     Subplot counter          
   nyqCut = size(splSpectra,2)/2;       subLoc = 1;              
%  Converting frequency to Strouhal 
   if St == 1 
      splFreq = splFreq*(Deq/Uj);      
   end
%  Spectra plotting Loop
   figure(2);
   for ctr = 1:length(angSeq)
%  FF Distance Correction
     if (ffDistnc ~= 44.715)
         splSpectra = splSpectra - abs(20*log10(ffDistnc/44.715));
     end                                        
%    Plot Spectra
     subplot(1,length(angSeq),ctr);     semilogx(splFreq(1:nyqCut),splSpectra(angSeq(ctr),1:nyqCut),'color',lineColr,'LineWidth',1.2,'DisplayName',condName);
     hold on;    box off;     grid on;  title(['$\psi = ',num2str(ffAngles(angSeq(ctr))),'^o$'],'Interpreter','latex');   aX = gca;
     cL = legend;             cL.EdgeColor = [1 1 1];     cL.Interpreter = 'latex';      cL.FontSize = 11;
     if St == 1 
        xlim([0.02 5.5]);               xlabel('$St$','Interpreter','latex');
     else
        xlim([500 10^5]);               xlabel('$Hz$','Interpreter','latex');
     end;    ylim([60 130]);            ylabel('$SPL(dB)$','Interpreter','latex');
     set(gca,'FontSize',textsize);      aX.TickLabelInterpreter = 'latex';          subLoc = subLoc+1; 
%    Peak Noise levels for each angle     
     [maxdB,freqLoc] = max(splSpectra(angSeq(ctr),1:nyqCut));
     disp([newline '->> Noise direction(angle)   - ', num2str(ffAngles(angSeq(ctr)))]);
     disp([        '->> Peak Noise Level(dB/SPL) - ', num2str(maxdB)]);
     disp([        '->> Peak Frequency(Hz/St)    - ', num2str(splFreq(freqLoc))]);   
     if subLoc == length(angSeq)
        subLoc = 1;
     end
   end;   figure(2);       set(gcf,'Position',[15 290 1520 375]); 
%  Integrated OASPL Calculation
   [oasplIntgrl] = oasplCalc(splSpectra);
%---------------------- Plotting Pressure Statistics in Polar Co-ordinates -------------------------
%  Convert angles to radians from degrees
   angleRad = deg2rad(ffAngles);
%  Plot Vertically 
   figure(3);               set(gcf,'Position',[810 130 585 420]);
   polarplot(angleRad,presStat(:,1),['-.',plotMrkr],'color',lineColr,'LineWidth',1.2,'MarkerSize',5,'MarkerFaceColor',mrkrFaceColr,'DisplayName',condName);
   thetalim([45 152]);                  rlim([100 135]);                hold on;    grid on;     pAx1 = gca; 
   pAx1.ThetaZeroLocation = 'bottom';   pAx1.ThetaDir = 'clockwise';
   rticks([100 110 120 130 140]);       thetaticks([45 60 75 90 105 116 132 152]); 
   pAx1.ThetaTickLabel = {'$45^o$', '$60^o$', '$75^o$', '$90^o$', '$105^o$', '$116^o$', '$132^o$', '$152^o$'};
   set(gca,'FontSize',textsize);        pAx1.TickLabelInterpreter = 'latex';        cL = legend; cL.EdgeColor = [1 1 1];   cL.Interpreter = 'latex';
   txtBox1 = text(1.33,145,'$Azimuthal \thinspace Angle(\psi)$','FontSize',textsize,'Interpreter','latex');  txtBox1.Rotation = 90;   
   txtBox2 = text(-3.24,118,'$OASPL(dB)$','FontSize',textsize,'Interpreter','latex');            cL.Location = 'southeastoutside';   
%  Plot Horizontally
   figure(4);               set(gcf,'Position',[95 345 790 420]);
   polarplot(angleRad,presStat(:,1),['-.',plotMrkr],'color',lineColr,'LineWidth',1.2,'DisplayName',condName,'MarkerSize',5,'MarkerFaceColor',mrkrFaceColr);
   thetalim([45 152]);                  rlim([100 135]);                hold on;    grid on;     pAx2 = gca; 
   pAx2.ThetaZeroLocation = 'left';     pAx2.ThetaDir = 'clockwise';
   rticks([100 105 110 115 120 125 130]);   thetaticks([45 60 75 90 105 120 136 152]);
   pAx2.ThetaTickLabel = {'$45^o$', '$60^o$', '$75^o$', '$90^o$', '$105^o$', '$120^o$', '$136^o$', '$152^o$'};
   set(gca,'FontSize',textsize);        pAx2.TickLabelInterpreter = 'latex';        cL = legend; cL.EdgeColor = [1 1 1];
   cL.Location = 'northeastoutside';    cL.Interpreter = 'latex';       cl.FontSize = 11;
   text(-3.3,114,'\it OASPL(dB)','FontSize',textsize,'Interpreter','latex');                   
   pAx2.OuterPosition = [0.0076   -0.1067    0.7282    1.0000];
   pAx2.Position      = [0.1023    0.0033    0.5644    0.8150];
   text(1.45,140,'$Azimuthal \thinspace Angle(\psi)$','FontSize',textsize,'Interpreter','latex');
 end
end
%% WRITING DATA TO A TXT FILE
if strcmp(write,'on')
%  OASPL & SPECTRUM FILES
   OASPL = data(:,1);
   drive_out = ['Z:\Updates\Stanford\Sent_Data\Acoustics_Data\Twin_Jet\OASPL\' nozzle{1} '\'];
   name_1 = ['OASPL_' condition{1} '_' nozzle{1}];
   name_2 = ['Spectrum_' condition{1} '_' nozzle{1}];
   cd(drive_out);
   dlmwrite([name_1 '.txt'],OASPL,'precision','%.4f');
   dlmwrite([name_2 '.txt'],PHI_yy_amp_avg,'precision','%.4f');
%  ANGLES FILE
   Angles = angles2';
   dlmwrite('Mic Angles.txt',Angles,'precision','%.4f');
end
%% FUNCTIONS:
% PRESSURE STATICTICS SELECTOR
function [statName,yVals] = statChoice(statNo)
statOrder  = [1;2;3;4;5]; % OASPL, P skewness, P kurtosis, dpdt Skewness, dpdt Kurtosis
statDisply = {'OASPL';'Pressure Skewness';'Pressure Kurtisos';'dpdt Skewness';'dpdt Kurtosis'};
figTitles  = {'OASPL';'Pressure \thinspace Skewness';'Pressure \thinspace Kurtosis';'dP/dt \thinspace Skewness';'dP/dt \thinspace Kurtosis'};
yAxTitles  = {'dB';'{\tilde \mu}_4(P)';'{\beta_4}(P)';'{\tilde \mu}_4(dP/dt)';'{\beta_4}(dP/dt)'};
varNames   = {'spl_pres_';'skew_pres_';'kurt_pres_';'skew_dpdt_';'kurt_dpdt_'};
yLimsUp    = [135;0.5;4;1;4;14];          yLimsDown = [100;-0.1;3;0;3];
inDx       = statNo == statOrder;         statName  = varNames{inDx};
yVals.limY = [yLimsDown(inDx) yLimsUp(inDx)]; 
yVals.plotName  = figTitles{inDx};        yVals.yName = yAxTitles{inDx}; 
disp(['->>',statDisply{inDx}]);
end
%% OASPL FROM SPL INTEGRATION
function [oasplIntgr] = oasplCalc(splSpectra)
% Power spectral density calculation
PSD = 10.^(splSpectra/10);                PSD =(PSD.*(20e-6)^2./50)./2;
% Frequency Range                         Integrated OASPL                        
df = 50:50:204800;                        oasplIntgr = zeros(size(PSD,1),size(PSD,2));
% Numerical Integration to compute OASPL
for ctr = 1:size(PSD,1)                         
  oasplIntgr(ctr,:) = trapz(df,PSD(ctr,:));   
  oasplIntgr(ctr,:) = 10*log10(oasplIntgr(ctr,:)/(20e-6)^2);
end
end
%% COLOR SELECTOR
function [condName,plotName,colorCtr,mrkrFaceColr,lineColr] = colorSelector(selectrInp)
% Possible nozzle configs                                                Index of the input nozzle
nozleConfigs = {'Major';'Minor';'Circular'};                             inDx1 = strcmp(selectrInp.nozzle,nozleConfigs);
% Index of condition in tests(relevant range: tests(1:7))
inDx2 = strcmp(selectrInp.tests,selectrInp.condition);
% FaceColor-Major:White,Circ:Plot color,Minor:Black                      MarkerFace Color
colorSet = [1,1,1;selectrInp.colr(inDx2,:);1,0,0];                       mrkrFaceColr = colorSet(inDx1,:); 
% Check for multi nozzle condition
if length(selectrInp.Nozzle) > 1 || selectrInp.condtnCtr == 1
   condName = ['NPR - ',num2str(selectrInp.NPR),':',selectrInp.nozzle];  plotName = selectrInp.yVals.plotName;
   colorCtr = selectrInp.colorCtr;
   if inDx1(2) == 1
      lineColr = selectrInp.colr(inDx2,:)*0.1;
   else
      lineColr = selectrInp.colr(inDx2,:);
   end
else
   condName = ['NPR - ',num2str(selectrInp.NPR)];                        lineColr = selectrInp.colr(selectrInp.colorCtr,:);               
   plotName = [selectrInp.yVals.plotName,' - ',selectrInp.nozzle];       colorCtr = selectrInp.colorCtr + 1;
end 
end
%% NOZZLE PLOT COMBINATIONS
% function [nozzle] = nozzleCombinations(plotFormat)
% switch plotFormat
%     case 'A'
%      nozzle = {'Minor','Major_Offset','Major'};           % FORMAT - 1  
%      colr = [0.1660 0.6740 0.1880;0.4940 0.1840 0.5560;0.8500 0.3250 0.0980];
%   case 'B'
%      nozzle = {'Major_Offset','Minor','Major'};           % FORMAT - 2 
%      colr = [0.4940 0.1840 0.5560;0.1660 0.6740 0.1880;0.8500 0.3250 0.0980];
%   case 'C'
%      nozzle = {'Minor','Major'};                          % FORMAT - 3 
% %      colr = [0 0.4470 0.7410;0 0.8 0.9];              %NPR 2.5
% %      colr = [0.4940 0.1840 0.5560;0.8 0.5 0.7];       %NPR 3.67
% %      colr = [0.4660 0.6740 0.1880;0.6 0.9 0.6];       %NPR 4.0
% %      colr = [0.6350 0.0780 0.1840;0.9 0.4 0.55];      %NPR 5.0
%      colr = [1 0 0;0 0 1];                      %COMPARISON TO AUTOCORRELATION
% %      colr = [0 0.5 0;0.8500 0.3250 0.0980];
% %      colr = [0.1660 0.6740 0.1880;0.8500 0.3250 0.0980];
%   case 'D'
%      nozzle = {'Minor'};
% %      colr = [0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];
% %      colr = [0 0.6 1; 0.9 0.3 0];
% %      colr = [1 0 0;0 0 1;0.6 0.6 0.6;];                      %COMPARISON TO AUTOCORRELATION
%   case 'E'
%      nozzle = {'Major'};
%      colr = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840];
% %      colr = [0 0.4470 0.7410;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840];
% %      colr = [0 0.4470 0.7410;0.4940 0.1840 0.5560;0.6350 0.0780 0.1840];
% %      colr = [0 0.6 1; 0.9 0.3 0];
%   otherwise
%      disp('Check Format');
% end
% end

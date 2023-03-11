tic
clearvars; clc; fclose all;
NPRs  = {'NPR_2p5';'NPR_2p9';'NPR_3p0';'NPR_3p5';'NPR_3p6'};
Temps = {'NTR_1p0';};
stats_name = {'spl_pres'};  stats_name2 = {'OASPL'};
title_npr = {'NPR = 2.5';'NPR = 2.9';'NPR = 3.0';'NPR = 3.5';'NPR = 3.67'};
temps = [1];        nprs = [3];                    %%%FLOW CONDITIONS
ang_seq = [16];     orint = 2;                  %%%ANGLE TO PLOT AND ORIENTATION 1 - MAJOR; 2 - MINOR
mode = 'St';                                  %%%PLOT Hz OR STROUHAL NO.
angles = [45 60 70 90 100 105 110 116 120 126 132 136 140 144 148 152];
NPR = str2double([NPRs{nprs}(5) '.' NPRs{nprs}(7)]);
NTR = str2double([Temps{temps}(5) '.' Temps{temps}(7)]);
if NPR == 3.6
   NPR = 3.67;
end
Uj = sqrt((NPR^(0.4/1.4)-1)*2/0.4)*sqrt(1.4*286.7*(NTR*297.15/(1+0.4/2*(sqrt((NPR^(0.4/1.4)-1)*2/0.4))^2)));
Matrix = zeros(16,4096);
mic = 0; ctr = 1;   Dsingl = 20.6502E-3;   Dtwin = 18.47E-3;
colr = [1 0.4 0.125; 0.494 0.184 0.556; 1 0 0; 0 0 0; 1 0 1; 0.2 0 0.8; 0.4 0 0.6; 0.6 0 0.4; 0.8 0 0.2;]; 

for y =1:2
  for j = nprs
  for i = temps
    if y == 1
       drive_root = 'Z:\Thesis\Far_Field_Data\Single_Jet\Results\V3_New_Chamber\';
    else
       drive_root = 'Z:\Thesis\Far_Field_Data\Twin_Jet\Results\';
    end
   angles2 = fliplr(angles);
   if orint == 1
      drive_in = [drive_root 'Major\'];
   else
      drive_in = [drive_root 'Minor\'];
   end
   file_name = ['FFT_Nick_nofilter_' Temps{i} '_' NPRs{j}];
   load([drive_in file_name]);   plotmics = ang_seq;
   for ch = plotmics
     mic = mic + 1;  
     if y == 1                                                                         
%         PHI_yy_amp_avg = PHI_yy_amp_avg;% + 20*log10(61.5/55); %+ 10*log10(Dsingl/Dtwin);      %SINGLE JET
        Matrix = PHI_yy_amp_avg;
     else
        PHI_yy_amp_avg = PHI_yy_amp_avg + 20*log10(61.5/55)+10*log10(Dsingl/Dtwin);     %TWIN JET
        Matrix = PHI_yy_amp_avg;
     end
    PSD = 10.^(Matrix/10);
    PSD2 =(PSD.*(20e-6)^2./50)./2;
    df = [50:50:204800];
    data = zeros(size(PSD,1),size(PSD,2));
    for a=1:size(PSD,1) %NUMERICAL INTEGRATION TO COMPUTE OASPL
      data(a,:) = trapz(df,PSD2(a,:));   
      data(a,:)=10*log10(data(a,:)/(20e-6)^2);
    end
    if strcmp(mode,'St') == 1
       if y == 1
          freq = freq*Dsingl/Uj;
       elseif y == 2
          freq = freq*(Dtwin)/Uj;
       end
       limit = [0.02 4.7];
       ax = 'St';
    else
        limit = [500 10^5];
        ax = 'Hz';
    end           
    figure(1)
     plot(angles2,data(:,1),'color',colr(ctr,:),'LineWidth',2);
     hold on;   box on;     
     grid on;   %title('OASPL','fontweight','bold');
     xlim([45 160]);    ylim([105 135]);
     xlabel('Mic Angle','fontweight','bold');      ylabel('OASPL(dB)','fontweight','bold');
     legend('Single Jet','Twin Jet');
     text = ['NPR ' num2str(NPR) ' TR ' num2str(NTR)];
     title(text);
     testsize=11;      set(gca,'FontSize',testsize);
     h = set(gcf,'Position',[550 300 400 350]);
     
    figure(2)
     semilogx(freq(1:4096/2),PHI_yy_amp_avg(ch,1:4096/2),'color',colr(ctr,:),'LineWidth',2);
     hold on;   box on;     grid on;
     xlim(limit);     ylim([60 130]);
     legend('Single Jet','Twin Jet');
     testsize=11;      set(gca,'FontSize',testsize);
     title(['\psi_o_{} = ' num2str(angles2(ch))]);
     xlabel(ax,'fontweight','bold');      ylabel('dB','fontweight','bold');
%      h = set(gcf,'Position',[850 300 390 350]);
     h = set(gcf,'Position',[350 400 350 330]);
    end
  end
  end
ctr = ctr +1;
end
hold off;
hold off;       set(0,'defaultfigurecolor',[1 1 1]);
toc
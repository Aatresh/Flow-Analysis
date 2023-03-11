%% PLOTS POD MODES FROM SCHLIEREN & SHADOWGRAPH IMAGES
%  FOR SINGLE & TWIN JETS
clearvars;  clc;    set(0,'defaultfigurecolor',[1 1 1]);  code = 'SchPODPlot';
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';
addpath([root1 'Jet_Analysis\Global_Functions\']); addpath([root1 'Jet_Analysis\POD_Codes\']);      
tests = {'NPR_2p5_TR_1p0', 'NPR_2p6_TR_1p0', 'NPR_2p9_TR_1p0','NPR_3p0_TR_1p0',...
         'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_5p0_TR_1p0', 'NPR_6p0_TR_1p0'};

condition = tests(3);   config = 'C';    nozzle = 'Minor';    NF = 'D';  

%  SELECTING DRIVE BASED ON FRAME RATE 
[OutputStruct] = GF_DriveSelect(config,nozzle,code); 
%  NORMALISATOIN FACTORS
Deq = OutputStruct.dt(2);        dt = OutputStruct.dt(1);
%% COLOR SCHEMA
%  ALL NPRS
colr = [0 0.4470 0.7410;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840;0.301 0.745 0.933] ;
%  EXCLUDING SCREECH CASES
% colr = [0 0.4470 0.7410;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840;0.301 0.745 0.933] ;
%  SCREECH CASES
% colr = [0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];
%  ENERGY AT SPECIFIC LOCATIONS
colr1 = [0 0 1;1 0 0;0 0 0;0 0.75 0.75;0.75 0 0.75;0 0.5 0;0.75 0.75 0];
%
for n = 1:length(condition)
  drive_in  = [OutputStruct.in_root '\' condition{n}(9:14) '\' condition{n}(1:7) '\'];
  drive_out = [OutputStruct.out_root '\' condition{n}(9:14) '\' condition{n}(1:7) '\'];
  [~,Uj,NPR,NTR] = GF_Velocity(condition{n});
  Norm_Fac.root1 = root1; Norm_Fac.typ = config(5:end);  Norm_Fac.op_cond = condition{n};
  Norm_Fac.code = code;   Norm_Fac.nozzle  = OutputStruct.nozzle;  Norm_Fac.NF = NF;  Norm_Fac.De = Deq;   
%  LOADING POD FILES
  load([drive_in 'modeEnergy-Rho']);   load([drive_in 'PODModes-Rho']);   load([drive_in 'Y']);
  load([drive_in 'tempCoeff-Rho']);    load([drive_in 'X']); 
  [Xn,Yn,lim_x,lim_y,xname,yname,lenScales,fig_siz]= GF_AxisDefnSch(config,nozzle,NF,X,Y);
%  NAMING CONVENTION
  if NPR == 3 || NPR == 4 || NPR == 5
     condName = ['NPR ' num2str(NPR) '.0, TR ' num2str(NTR) '.0'];
  else
     condName = ['NPR ' num2str(NPR) ', TR ' num2str(NTR) '.0'];       
  end
  if strcmp(nozzle,'Minor') || strcmp(nozzle,'Circular')
     ptClr = 'auto';
  else
     ptClr = colr(n,:); 
  end
%  MODAL ENERGY DISTRIBUTION & CUMULATIVE SUM OF ALL MODAL ENERGIES
  edisX = (lambdaRho/sum(lambdaRho))*100;  lambdaSum = zeros(size(lambdaRho,1),1);
  for ctr3 = 1:length(lambdaRho)
    lambdaSum(ctr3) = sum(edisX(1:ctr3));
  end
  prompt = [newline '-> Show Mode Metrics?(y/n) - '];    plt_mtrcs = input(prompt,'s');
  if strcmp(plt_mtrcs,'y')
     figure(1);    subplot(211);    plot(edisX,'color',colr(n,:),'LineWidth',1.2,'DisplayName',condName);
     set(gca,'FontSize',13);        ylabel('$(\lambda_n / \Sigma\lambda_n)$','Interpreter','latex');   
     grid on;   box off;  hold on;  xlabel('$Mode(n)$','Interpreter','latex');          ax = gca;
     title('$Mode \thinspace Energy \thinspace Distribution$','Interpreter','latex');   ax.TickLabelInterpreter = 'latex';
     subplot(212);                  plot(lambdaSum,'color',colr(n,:),'LineWidth',2,'DisplayName',condName), 
     set(gca,'FontSize',13);        ylabel('$(\Sigma\lambda_n)$','Interpreter','latex'),ax = gca;
     grid on,   box off;  hold on,  xlabel('$n$','Interpreter','latex');     ax.TickLabelInterpreter = 'latex';          
     title('$Cumulative \thinspace Energy$','Interpreter','latex');          set(gcf,'Position',[35 220 550 625]); 
     title_lm = ['$Mode \thinspace Energy - ',nozzle,'$'];
     figure(2);    loglog(lambdaRho,'--o','Color',colr(n,:),'LineWidth',1.2,'DisplayName',condName,'MarkerFaceColor',ptClr,'MarkerEdgeColor',colr(n,:));
     set(gca,'FontSize',13);        ylabel('$\lambda_{\Delta\rho}$','Interpreter','latex'),   ax = gca;
     grid on,   box off,  hold on;  xlabel('$n$','Interpreter','latex'),     xlim([0 100]),   ylim([4.5*10^9 1.3*10^11]); 
     title(title_lm,'Interpreter','latex');         set(gcf,'Position',[70 220 520 420]);     ax.TickLabelInterpreter = 'latex'; 
     prompt = [newline 'SHOW LEGENDS?(y/n) - '];    choice = input(prompt,'s');
     if strcmp(choice,'y')
        figure(1); subplot(211);legend('show'), legend boxoff; figure(1);subplot(212);legend('show'), legend boxoff;
        figure(2);legend('show'), legend boxoff;
     end
  end
%% PLOTTING POD MODES
    prompt = [newline '-> Plot POD Modes?(y/n) - '];    plt_mods = input(prompt,'s');
    if strcmp(plt_mods,'y')
       prompt = [newline '--> Enter no. of Modes to plot: '];  mode = input(prompt);
       for ctr = 1:mode
         if strcmp(nozzle,'Minor') || strcmp(nozzle,'Circular')
            name = ['$\phi_{\Delta\rho y}:Mode-' num2str(ctr) '$'];
         else
            name = ['$\phi_{\Delta\rho z}:Mode-' num2str(ctr) '$'];
         end
         figure;    pcolor(Xn,Yn,phiRho(:,:,ctr));    shading interp;    caxis([-0.5 0.5]);    colormap bluewhitered;  
         colorbar,  box off;   axis equal;            xlim(lim_x),       ylim(lim_y);          set(gca,'FontSize',13)            
         xlabel(xname,'Interpreter','latex');         ylabel(yname,'Interpreter','latex');     ax = gca;          colorbar off;  
         ax.TickLabelInterpreter = 'latex';           title(name,'Interpreter','latex');       set(gcf,'Position',[160 190 590 420]);
%        cd(drive_out);    saveas(figure(y+m*5),['Mode_' num2str(y) '.png']);   cd(drive_code);
%        set(gcf,'Position',[425 468 540 273]);
%        set(gcf,'Position',[198 501 756 313]);
%        set(gcf,'Position',[173 498 938 399]); 
%        set(gcf,'Position',[550 350 857 396]);
       end
    end
 %% PLOTTING MAGNITUDE ALONG A LINE
    prompt = [newline '-> Show POD Modes along a line(y/n): '];  choice = input(prompt,'s');  
    if strcmp(choice,'y')
       figure(500);    set(gcf,'Position',[110 270 700 340]);    ctr = 1;
       while(strcmp(choice,'y'))
         prompt = [newline '--> Enter Y-Location: '];    yLoc  = input(prompt);
      	 prompt = [newline '--> Mode no.: '];            modeNo = input(prompt);
         [~,ctr1] = find(Yn<yLoc);  pltNname = [yname(2:4) ' = ',num2str(round(Yn(ctr1(end)+1),2))];
         modeVal = smoothdata(phiRho(ctr1(end)+1,:,modeNo),2,'gaussian',2);
         plot(Xn,modeVal,':','color',colr1(ctr,:),'LineWidth',1.2,'DisplayName',pltNname);   
         hold on;   grid on;   box off;    ylim([-1 1]);    xlim(lim_x);   ctr = ctr+1; 
         prompt = [newline '-> Plot another Location(y/n): '];    choice = input(prompt,'s');
       end
       xlabel(xname,'Interpreter','latex');ylabel('$Mode \thinspace Energy \thinspace Intensity$','Interpreter','latex');
       set(gca,'FontSize',12,'fontname','times'); legend('show'), legend boxoff;      
     end
%% PLOTTING PSD OF FIRST 5 MODES
    prompt = [newline 'PLOT MODE PSDS?(y/n) - '];    choice = input(prompt,'s');
    if strcmp(choice,'y')
       PSD = [];   fs = 1/dt; 
%  COMPUTING PSD's FOR THE FIRST 20 MODES
       for ctr2 = 1:20
         [PSD(ctr2,:),f] = cpsd(tempRho(ctr2,:),tempRho(ctr2,:),[],[],[],fs);
       end;          St = f*Deq/(Uj*1000);         figure(12)
       set(gcf,'Position',[460 50 1040 375]);
       for ctr2 = 1:5
         subplot(121); loglog(St,PSD(ctr2,:),'LineWidth',1.5,'DisplayName',['Mode ',num2str(ctr2)]);hold on;
         subplot(122); loglog(f,PSD(ctr2,:),'LineWidth',1.5,'DisplayName',['Mode ',num2str(ctr2)]); hold on;
       end
       low = min(min(min(PSD)));    high = max(max(max(PSD)));      subplot(121);
       box off; grid on; legend('show','Location','northwest'); legend boxoff;  xlabel('St');   ylabel('PSD_{aa}');
       xlim([0 St(end)]); ylim([low/5 high*5]); set(gca,'FontSize',12,'fontname','times');
       subplot(122); box off;  grid on;    legend('show','Location','northwest'); legend boxoff;  xlabel('Hz');   ylabel('PSD_{aa}');
       xlim([0 f(end)]); ylim([low/5 high*5]);    set(gca,'FontSize',12,'fontname','times');
%% POD PHASE PORTRAIT
    m1 = 1;     m2 = 2;     inc = 1;  Rm = mean(sqrt(tempRho(m1,:).^2+tempRho(m2,:).^2));   step = 1;
    figure(13);   subplot(131)
    plot(tempRho(m1,1:step:end)/Rm,tempRho(m2,1:step:end)/Rm,'o','MarkerSize',3,'MarkerFaceColor',colr(1,:),'MarkerEdgeColor',colr(1,:));
    grid on;    axis equal, set(gca,'FontSize',12,'FontName','times'),    xlim([-1.4 1.4]),   ylim([-1.4 1.4]);
    xlabel(['a_',num2str(m1)],'FontWeight','bold'), ylabel(['a_',num2str(m2)],'FontWeight','bold'), hold on;
    viscircles([0 0],1); box off; 
    subplot(133); plot(tempRho(m2,1:step:end)/Rm,tempRho(m2+inc,1:step:end)/Rm,'s','MarkerSize',6,'MarkerFaceColor',colr(2,:),'MarkerEdgeColor',colr(2,:));
    grid on;box off; axis equal, set(gca,'FontSize',12,'FontName','times'),    xlim([-1.4 1.4]),   ylim([-1.4 1.4]);
    xlabel(['a_',num2str(m2)],'FontWeight','bold'), ylabel(['a_',num2str(m2+inc)],'FontWeight','bold')
    subplot(132); plot(tempRho(m1,1:step:end)/Rm,tempRho(m2+inc,1:step:end)/Rm,'d','MarkerSize',6,'MarkerFaceColor',colr(4,:),'MarkerEdgeColor',colr(4,:));
    grid on,box off,axis equal, set(gca,'FontSize',12,'FontName','times'),    xlim([-1.4 1.4]),   ylim([-1.4 1.4]);
    xlabel(['a_',num2str(m1)],'FontWeight','bold'), ylabel(['a_',num2str(m2+inc)],'FontWeight','bold'), set(gcf,'Position',[70 315 1400 425]);
    end
%% PLOT RECONSTRUCTED IMAGES
%     if strcmp(plt,'y')
%        figure(14)
%        set(gcf,'Position',[113 144 934 796]);   j = 5;
%        for ctr = 1:4
%          subplot(2,2,ctr)
%          pcolor(X_sc,Y_sc,U_rec(:,:,ctr+j)), shading interp, colormap gray;
%          xlim([0 14]), ylim([-6 6]), axis equal, xlabel('X/D_e'); ylabel(yname);
%          xticks([0 2 4 6 8 10 12 14]);
%          set(gca,'FontSize',12,'FontName','times');
%        end
% %%  ANIMATE IMAGES
%        figure(15)
%        set(gcf,'Position',[136 549 805 290]);
%        for ctr = 1:100
%          pcolor(X_sc,Y_sc,squeeze(U_rec(:,:,ctr))); shading interp, axis equal, colormap gray;
%          xlabel('X/D_{eq}');   ylabel(yname);     xlim([0 14]);    ylim([-2 2]);   %caxis([-0.5 0.2]);  
%          set(gca,'FontSize',12,'FontName','times');
%          drawnow
%          frame = getframe(15);
%          im{ctr} = frame2im(frame);                          %#ok
%        end
% %  WRITING TO GIF 
%        prompt = [newline '>> SAVE GIF?(y/n) - '];    gifwrite_img = input(prompt,'s'); 
%        if strcmp(gifwrite_img,'y')
%           cd(drive_out);       gifname = [condition{1} '.gif'];
%           for idx = 1:size(im,2)
%             [A,map] = rgb2ind(im{idx},256);
%             if idx == 1
%                imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.1);
%             else
%                imwrite(A,map,gifname,'WriteMode','append','DelayTime',0.1);
%             end
%           end
%        end
%     end
 end

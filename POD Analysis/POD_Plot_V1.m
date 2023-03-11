%% PLOT POD MODES AND RECONSTRUCTUED VELOCITY IMAGES
clc;    clearvars; fclose all;  set(0,'defaultfigurecolor',[1 1 1]);   code = 'PIV_PODPlot';
root = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\'; 
cd([root 'Jet_Analysis\Global_Functions\']); addpath([root 'Jet_Analysis\POD_Codes\']);       
tests = {'NPR_2p5_TR_1p0','NPR_2p9_TR_1p0','NPR_3p0_TR_1p0', 'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0','NPR_5p0_TR_1p0'};

condition = [tests(2)];     nozzle = {'Minor'};       jet_typ = 'S';        textsize = 13;   

%% COLORMAP 
% ALL NPRS
% colr = [0 0.4470 0.7410;0 0 0;1 0 0;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.301 0.745 0.933] ;
% EXCLUDING SCREECH CASES
% colr = [0 0.4470 0.7410;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.301 0.745 0.933] ;
% SCREECH CASES
% colr = [0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];
% colr = [0 0 0;1 0 0];
% JASA 
 colr = [0 0.75 0.75;1 0 0;0 0 1; 0.3 0.3 0.3;0.4660 0.6740 0.1880];
% SIZES
% set(gcf,'Position',[853 386 523 558]);
%% WINDOW SIZE - JASA
% CUMILATIVE & ENERGY DISTRIBUTN. - set(gcf,'Position',[665 500 520 320]);
% xlim([1 10]); ylim([0.1 10.5]); xticks([1 2 4 6 8 10]);
% POD MODE - ylim([-1.5 1.5]); xticks([0.5 1 2 3 4 5]); set(gcf,'Position',[149 332 617 611]);
% RECONSTRUCTED VELOCITY (CONVEC, PHASED, STREAMLINES) - ylim([-1.2 1.2]); set(gcf,'Position',[164 395 643 308]);
% caxis([0 0.01]); xticks([0.5 1 2 3 4 5]);
% STREAMLINES - set(gcf,'Position',[79 372 467 399]);
%% PARAMETERS TO PLOT
mode_versn = 1;   Lgnd = 'on';  lambdax_all = [];  lambday_all = []; Avg_U = []; Avg_V = [];load vorticity_map.mat;
%  LIST OF CONFIGURATIONS
%  COMMON    - Major, Minor         for epsilon in ylabel use - char(949)
%  TWIN ONLY - Minor_J2, Sym_Lin
switch jet_typ
  case 'S'
    drive_root = [root 'POD_Data\PIV - Single_Rectan\'];
    Deq = 20.6502;  diff = 5.63;
  case 'T'
    drive_root = [root '\POD_Data\Twin_Jet\'];
    Deq = 18.47;   diff = 5;
end

for m = 1:length(nozzle)
for n = 1:length(condition)
  [~,Uj,NPR,TR] = GF_Velocity(condition{n});
  drive_in = [drive_root nozzle{m} '\' condition{n}(9:14) '\' condition{n}(1:7) '\'];
% LOADING DATA
  load([drive_in 'Modal_Ampl_X']);  load([drive_in 'mode_energy_X']);  load([drive_in 'Recon_var_X']);  load([drive_in 'POD_Modes_X']);
  load([drive_in 'Modal_Ampl_Y']);  load([drive_in 'mode_energy_Y']);  load([drive_in 'Recon_var_Y']);  load([drive_in 'POD_Modes_Y']);
  load([drive_in 'X']);  load([drive_in 'Y']);    load([drive_in 'AvgVx']);  load([drive_in 'AvgVy']);
  load('custom_map1.mat');   
  X1 = X + abs(X(1))+ diff;         
  X1 = X1/Deq;     Y1 = Y/Deq;
  if strcmp(jet_typ,'S')
     [jet_edge,~] = Find_Center_SJ_POD(U_bar,Y,Uj);     centr = jet_edge(1);
  else
     [jet_edge,~] = Find_Center_TJ_POD(U_bar,Y,Uj,nozzle{m});
     if strcmp(nozzle{m},'Minor') || strcmp(nozzle{m},'Sym_Lin')
        centr = jet_edge(1);
     else
        centr = round((jet_edge(1)+jet_edge(4))/2);
     end
  end
% CENTERING Y AXIS
  Avg_U = U_bar/Uj;     Avg_V = V_bar/Uj;
  Y1 = Y1 - Y1(centr);
  lambdax_all(1:100,n,m) = lambda_x(1:100);                 %#ok
  lambday_all(1:100,n,m) = lambda_y(1:100);                 %#ok
% MODAL ENERGY
  prompt = [newline '>> PLOT MODAL ENERGY?(y/n) - ']; modal_energy = input(prompt,'s');
  if strcmp(modal_energy,'y')
% ENERGY PERCENT PER MODE
     edis_x = (lambda_x/sum(lambda_x))*100;      esum_x = [];
     edis_y = (lambda_y/sum(lambda_y))*100;      esum_y = [];
% MODAL ENERGY PERCENT SUMMATION - SUMMED UPTO THE Nth MODE
     for ctr = 1:size(edis_x,1)
       esum_x(ctr) = sum(edis_x(1:ctr));                    %#ok
       esum_y(ctr) = sum(edis_y(1:ctr));                    %#ok
     end
% CONDITION DISPLAY TEXT
    if NPR == 3 || NPR == 4 || NPR == 5
       N_pr = ['NPR ' num2str(NPR) '.0, TR ' num2str(TR) '.0'];
    else
       N_pr = ['NPR ' num2str(NPR) ', TR ' num2str(TR) '.0'];       
    end
    if strcmp(nozzle{m},'Major')
       f_clr = 'auto';
    else
       f_clr = colr(n,:);
    end
    if length(condition)>1
       name_x = [N_pr ' - ' nozzle{m}];
       name_y = [N_pr ' - ' nozzle{m}];
       title_x1 = ['Stream-wise Modal Energy(\lambda_x) - ' nozzle{m}];
       title_x2 = ['Cumilative Stream-wise Modal Energy(\Sigma\lambda_{nx}) - ' nozzle{m}];
       if strcmp(nozzle{m},'Minor')
          title_y1 = ['Cross-stream Modal Energy(\lambda_y) - ' nozzle{m}];
          title_y2 = ['Cumilative Cross-stream Modal Energy(\Sigma\lambda_{ny}) - ' nozzle{m}];
       else
          title_y1 = ['Cross-stream Modal Energy(\lambda_z) - ' nozzle{m}];
          title_y2 = ['Cumilative Cross-wise Modal Energy(\Sigma\lambda_{nz}) - ' nozzle{m}];
       end      
       figure(1)
       if strcmp(nozzle{m},'Minor') && NPR == 3.67 
          rng = 2:length(lambda_x);
       else
          rng = 1:length(lambda_x);
       end
       loglog(lambda_x(rng),'--s','Color',colr(n,:),'DisplayName',name_x,'MarkerFaceColor',f_clr,'MarkerEdgeColor',colr(n,:)),grid on,box off,hold on;
       loglog(lambda_y(rng),'--o','Color',colr(n,:),'DisplayName',name_y,'MarkerFaceColor',f_clr,'MarkerEdgeColor',colr(n,:)),grid on,box off,hold on;
       ylabel('\bf \lambda_n'),  xlabel('\bf n'), ylim([10^6 10^10]), xlim([0 1000]);
       title(['Mode Energy - ' nozzle{m}]), set(gca,'FontName','times','FontSize',textsize);
       if strcmp(Lgnd,'on')
          c = legend('show'); c.EdgeColor = [1 1 1];        
       end
    else
       name_x = '\bf \lambda_x';
       if strcmp(nozzle{m},'Minor')
          name_y = '\bf \lambda_y';
          title_y1 = ['Cross-stream Modal Energy(\lambda_y) - ' nozzle{m}];
          title_y2 = ['Cumilative Stream-wise Modal Energy(\Sigma\lambda_{ny}) - ' nozzle{m}];
       else
          name_y = '\bf \lambda_z';
          title_y1 = ['Cross-stream Modal Energy(\lambda_z) - ' nozzle{m}];
          title_y2 = ['Cumilative Stream-wise Modal Energy(\Sigma\lambda_{nz}) - ' nozzle{m}];
       end
       title_x1 = ['\lambda_x - ' nozzle{m}];  title_x2 = ['Cumilative Stream-wise Modal Energy(\Sigma\lambda_{nx}) - ' nozzle{m}];
       
       figure(1)
       loglog(lambda_x,'--s','Color',colr(n,:),'DisplayName',name_x,'MarkerFaceColor',f_clr,'MarkerEdgeColor',colr(n,:)),grid on,box off,hold on;
       loglog(lambda_y,'--o','Color',colr(n,:),'DisplayName',name_y,'MarkerFaceColor',f_clr,'MarkerEdgeColor',colr(n,:));
       ylabel('\bf \lambda_n'),  xlabel('\bf n'), ylim([10^7 10^9]), xlim([0 100]);
       title([N_pr ' - ' nozzle{m}]), set(gca,'FontName','times','FontSize',textsize);
       if strcmp(Lgnd,'on')
          c = legend('show'); c.EdgeColor = [1 1 1];        
       end
    end
    figure(2)
    plot(edis_x,'--s','Color',colr(n,:),'LineWidth',2,'DisplayName',name_x,'MarkerFaceColor',f_clr,'MarkerEdgeColor',colr(n,:)),grid on,box off,hold on;
    ylabel('\bf % \lambda_n / \Sigma\lambda_n'),  xlabel('\bf n'), xlim([1 6]), ylim([ 0 5]); title(title_x1);
    set(gca,'FontName','times','FontSize',textsize); 
    if strcmp(Lgnd,'on')
       legend('show'); legend boxoff;   
    end
    figure(3)
    plot(edis_y,'--o','Color',colr(n,:),'LineWidth',2,'DisplayName',name_y,'MarkerFaceColor',f_clr,'MarkerEdgeColor',colr(n,:)),hold on, grid on, box off;
    ylabel('\bf % \lambda_n / \Sigma\lambda_n'),  xlabel('\bf n'),    xlim([1 6]), ylim([ 0 5]); title(title_y1);
    set(gca,'FontName','times','FontSize',textsize);     
    if strcmp(Lgnd,'on')
       c = legend('show'); c.EdgeColor = [1 1 1];   
    end
    figure(4)
    plot(esum_x,'--s','Color',colr(n,:),'LineWidth',2,'DisplayName',name_x,'MarkerFaceColor',f_clr,'MarkerEdgeColor',colr(n,:)), grid on, hold on, box off;
    ylabel('% \Sigma\lambda_n','FontWeight','bold'), xlabel('n','FontWeight','bold'); title(title_x2);
    set(gca,'FontName','times'),set(gca,'FontSize',textsize);     % set(gcf,'Position',[106 511 494 416]);
    if strcmp(Lgnd,'on')
       c = legend('show'); c.EdgeColor = [1 1 1];   
    end
    xlim([1 6]), ylim([0 25]);
    figure(5)
    plot(esum_y,'--o','Color',colr(n,:),'LineWidth',2,'DisplayName',name_y,'MarkerFaceColor',f_clr,'MarkerEdgeColor',colr(n,:)), grid on, hold on, box off;
    ylabel('% \Sigma\lambda_n','FontWeight','bold'), xlabel('n','FontWeight','bold'); title(title_y2);
    set(gca,'FontName','times'),set(gca,'FontSize',textsize);    % set(gcf,'Position',[106 511 494 416]);
    if strcmp(Lgnd,'on')
       c = legend('show'); c.EdgeColor = [1 1 1];   
    end
    xlim([1 6]), ylim([0 25]);
    
  end
  if strcmp(nozzle{m},'Minor')
     yname = 'Y/D_e';
  else
     yname = 'Z/D_e';
  end
  % set(gcf,'Position',[267 531 567 275]);
%% PLOT POD MODES - VERSION 1
  prompt = [newline '>> PLOT POD MODES?(y/n) - ']; pod_modes = input(prompt,'s');
  if strcmp(pod_modes,'y')
     prompt = [newline '>> NUMBER OF MODES TO PLOT - '];  mode = input(prompt);    ctr = 0; 
     if mode_versn == 1
        while ctr<= mode
         for y = 1:2
           name = ['Mode ' num2str(y+ctr)];
           figure(5+ctr+1)
           subplot(2,1,y)
           pcolor(X1,Y1,phi_x(:,:,y+ctr));  shading interp;     axis equal; caxis([-1 1]);    colorbar;     colormap bluewhitered;
           xlabel('X/D_e');    ylabel(yname);       set(gca, 'Ticklength', [0 0])
           grid off;          box off;   lgd = legend('Vx');       legend boxoff;       lgd.FontWeight = 'bold'; 
           set(gca,'FontName','times','FontSize',textsize);          
           xlim([0.29 5]);     title(name);  set(gcf,'Position',[184 47 750 871]);   
           figure(5+ctr+2)
           subplot(2,1,y)
           pcolor(X1,Y1,phi_y(:,:,y+ctr));  shading interp;     axis equal; caxis([-1 1]);    colorbar;     colormap bluewhitered;
           xlabel('X/D_e');    ylabel('Y/D_e');       set(gca, 'Ticklength', [0 0])
           grid off;          box off;   lgd = legend('Vy');       legend boxoff;       lgd.FontWeight = 'bold';   
           set(gca,'FontName','times','FontSize',textsize);   
           xlim([0.3 5]);     title(name);  set(gcf,'Position',[184 47 750 871]);   
         end
         ctr = ctr+2;
        end
        % SIZING FILTERS FOR PUBLICATIONS
        %  set(gcf,'Position',[550 550 550 350]);
        %  set(gcf,'Position',[87 -11 647 1007]);
        %  set(gcf,'Position',[87 271 508 682]);
     else
%% PLOT POD MODES - VERSION 2
        ctr = 0;  ctr1 = 0;
        while ctr<mode
         ctr = ctr+1;
         if mode<8
            ctr1 = ctr;
         else
            ctr1 = ctr1+1;
         end
        name_x = ['\phi_{x' num2str(ctr) '}'];
        if strcmp(nozzle{m},'Minor')
           name_y = ['\phi_{y' num2str(ctr) '}'];         
        else
           name_y = ['\phi_{z' num2str(ctr) '}'];          
        end
        figure(100)
        subplot(2,4,ctr1)
        pcolor(X1,Y1,phi_x(:,:,ctr)), axis equal, shading interp, caxis([-1 1]), colormap bluewhitered;colorbar;
        xlabel('X/D_e'), ylabel(yname), grid off, box off; 
        ylim([-1.5 1.5]), xlim([0.3 5]), title(name_x), set(gca,'FontName','times','FontSize',textsize);   
        set(gcf,'Position',[84 230 2110 733]);
        figure(101)
        subplot(2,4,ctr1)
        pcolor(X1,Y1,phi_y(:,:,ctr)), axis equal, shading interp, caxis([-1 1]), colormap bluewhitered;colorbar;
        xlabel('X/D_e'), ylabel(yname), grid off, box off; 
        ylim([-1.5 1.5]), xlim([0.3 5]), title(name_y), set(gca,'FontName','times','FontSize',textsize);   
        set(gcf,'Position',[84 230 2110 733]);
        end
    end
%% PHASE PORTRAIT 
%   AXIAL VELOCITY
    m1 = 1;     m2 = 2;
    rm = mean(sqrt(aVx(m1,:).^2+aVx(m2,:).^2));       step = 2;
    figure(102)
    subplot(131)
    plot(aVx(m1,1:step:end)/rm,aVx(m2,1:step:end)/rm,'s','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');
    grid on;    axis equal, set(gca,'FontSize',textsize,'FontName','times'),    xlim([-1.4 1.4]),   ylim([-1.4 1.4]);
    xlabel('a_{x1}','FontWeight','bold'), ylabel('a_{x2}','FontWeight','bold')
    subplot(133)
    plot(aVx(m2,1:step:end)/rm,aVx(m2+1,1:step:end)/rm,'s','MarkerSize',6,'MarkerFaceColor',colr(2,:),'MarkerEdgeColor',colr(2,:));
    grid on;    axis equal, set(gca,'FontSize',textsize,'FontName','times'),    xlim([-1.4 1.4]),   ylim([-1.4 1.4]);
    xlabel('a_{x2}','FontWeight','bold'), ylabel('a_{x3}','FontWeight','bold')
    subplot(132)
    plot(aVx(m2+1,1:step:end)/rm,aVx(m1,1:step:end)/rm,'s','MarkerSize',6,'MarkerFaceColor',colr(4,:),'MarkerEdgeColor',colr(4,:));
    grid on,axis equal, set(gca,'FontSize',textsize,'FontName','times'),    xlim([-1.4 1.4]),   ylim([-1.4 1.4]);
    xlabel('a_{x1}','FontWeight','bold'), ylabel('a_{x3}','FontWeight','bold'), set(gcf,'Position',[100 162 1540 391])
%   RADIAL VELOCITY
    rm = mean(sqrt(aVy(1,:).^2+aVy(2,:).^2));   
    figure(103)
    subplot(131)
    plot(aVy(m1,1:step:end)/rm,aVy(m2,1:step:end)/rm,'o','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','k');
    grid on;    axis equal, set(gca,'FontSize',textsize,'FontName','times'),    xlim([-1.4 1.4]),   ylim([-1.4 1.4]);
    xlabel('a_{y1}','FontWeight','bold'), ylabel('a_{y2}','FontWeight','bold')
    subplot(133)
    plot(aVy(m2,1:step:end)/rm,aVy(m2+1,1:step:end)/rm,'o','MarkerSize',6,'MarkerFaceColor',colr(2,:),'MarkerEdgeColor',colr(2,:));
    grid on;    axis equal, set(gca,'FontSize',textsize,'FontName','times'),    xlim([-1.4 1.4]),   ylim([-1.4 1.4]);
    xlabel('a_{y2}','FontWeight','bold'), ylabel('a_{y3}','FontWeight','bold')
    subplot(132)
    plot(aVy(m2+1,1:step:end)/rm,aVy(m1,1:step:end)/rm,'o','MarkerSize',6,'MarkerFaceColor',colr(4,:),'MarkerEdgeColor',colr(4,:));
    grid on,axis equal, set(gca,'FontSize',textsize,'FontName','times'),    xlim([-1.4 1.4]),   ylim([-1.4 1.4]);
    xlabel('a_{y1}','FontWeight','bold'), ylabel('a_{y3}','FontWeight','bold'), set(gcf,'Position',[100 162 1540 391])
  end
%% IMAGE RECONSTRUCTION
  prompt = [newline '>> PLOT RECONSTRUCTED IMAGES?(y/n) - ']; rec = input(prompt,'s');
  if strcmp(rec,'y')
     nImg = size(phi_x,3);                                %%%NUMBER OF IMAGES AND MODES TO PLOT
     if strcmp(nozzle{m},'Minor')||strcmp(nozzle{m},'Major') && NPR == 3.0 || NPR == 2.9
        mode = 1:2;     lims1 = [0 0.015];   lims2 = [0 0.1];
     else
        mode = 1:size(phi_x,3); lims1 = [0 0.03];   lims2 = [0 0.1];  
     end
     Vel_rec = zeros(size(beta_U_rec,1),nImg,2); 
     for ctr1 = 1:nImg
       Vel_rec(:,ctr1,1) = beta_U_rec(:,mode)*aVx(mode,ctr1);
       Vel_rec(:,ctr1,2) = beta_V_rec(:,mode)*aVy(mode,ctr1);
     end
     U_prim1 = reshape(Vel_rec(:,:,1),[size(phi_x,1) size(phi_x,2) nImg]);
     U_prim_R = (mean((U_prim1).^2,3));
     V_prim1 = reshape(Vel_rec(:,:,2),[size(phi_y,1) size(phi_y,2) nImg]);
     V_prim_R = (mean((V_prim1).^2,3));
     TKE_rec1 = ((U_prim_R+2*V_prim_R)*0.5)/Uj^2;
     TKE_rec2 = sqrt(TKE_rec1(:,:,n)*Uj^2)/Uj;
     figure(104);   subplot(211);   title('TKE'); set(gcf,'Position',[1028 209 534 740]);
     pcolor(X1,Y1,TKE_rec1);shading interp,colormap jet,axis equal;colorbar;caxis(lims1); 
     xlim([0.3 5]);  ylim([-2 2]);  xlabel('X/D_e'), ylabel(yname); set(gca,'FontSize',textsize,'FontName','times');
     subplot(212)
     pcolor(X1,Y1,TKE_rec2);shading interp;colormap jet;axis equal;colorbar;caxis(lims2); xlim([0.3 5]);
     figure(105);   subplot(211);   set(gcf,'Position',[1028 209 534 740]);
     pcolor(X1,Y1,rms(U_prim1,3)/Uj);shading interp;colormap jet;axis equal;c = colorbar;caxis([0 0.1]); 
     xlim([0.3 5]);  ylim([-2 2]);  xlabel('X/D_e'), ylabel(yname); set(gca,'FontSize',textsize,'FontName','times');
     title('RMS of Fluctuating Velocity'); title(c,'V"_{rms}/Uj');
     subplot(212)
     pcolor(X1,Y1,rms(V_prim1,3)/Uj);shading interp;colormap jet;axis equal; colorbar;caxis([0 0.1]); 
     xlim([0.3 5]);  ylim([-2 2]);  xlabel('X/D_e'), ylabel(yname); set(gca,'FontSize',textsize,'FontName','times');
     U_inst = U_prim1 + U_bar;          V_inst = V_prim1 + V_bar;
     % PLOTTING INSTANTANEOUS IMAGES
     prompt = [newline '>> NUMBER OF IMAGES TO PLOT - ']; img_no = input(prompt); 
     prompt = [newline '>> IMAGES OR GIF?(1/2) - '];     chc    = input(prompt);
     if chc == 1
        for ctr2 = 1:img_no
          figure(ctr2*200); set(gcf,'Position',[87 420 1139 369]); subplot(121)
          pcolor(X1,Y1,U_inst(:,:,ctr2));shading interp;colormap jet;caxis([0 420]); axis equal; 
          xlim([0.3 5]);  ylim([-2 2]);  xlabel('X/D_e'), ylabel(yname); set(gca,'FontSize',textsize,'FontName','times');
          subplot(122)
          pcolor(X1,Y1,V_inst(:,:,ctr2));shading interp;colormap jet; caxis([-30 30]);axis equal; 
          xlim([0.3 5]);  ylim([-2 2]);  xlabel('X/D_e'), ylabel(yname); set(gca,'FontSize',textsize,'FontName','times');
        end
     else
        for ctr2 = 1:img_no
          figure(106); set(gcf,'Position',[87 420 1139 369]); subplot(121)
          pcolor(X1,Y1,U_inst(:,:,ctr2));shading interp;colormap jet;caxis([0 420]); axis equal; 
          xlim([0.3 5]);  ylim([-2 2]);  xlabel('X/D_e'), ylabel(yname); set(gca,'FontSize',textsize,'FontName','times');
          subplot(122)
          pcolor(X1,Y1,V_inst(:,:,ctr2));shading interp;colormap jet; caxis([-30 30]);axis equal; 
          xlim([0.3 5]);  ylim([-2 2]);  xlabel('X/D_e'), ylabel(yname); set(gca,'FontSize',textsize,'FontName','times');
          pause(0.2)  
          drawnow
          frame = getframe(106);
          im{ctr2} = frame2im(frame);
        end
        cd(drive_in);       gifname = [condition{1} '.gif'];
        for idx = 1:size(im,2)
          [A,map] = rgb2ind(im{idx},256);
          if idx == 1
             imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.1);
          else
             imwrite(A,map,gifname,'WriteMode','append','DelayTime',0.1);
          end
        end
        clear im;
     end
     if strcmp(nozzle{m},'Minor')|| strcmp(nozzle{m},'Major') && NPR == 3.0 || NPR == 2.9   
     %  CONVECTIVE VELOCITY COMPONENTS
         U_conv = U_prim1;      V_conv = V_prim1;   
        figure(107)
%       subplot(211)
        pcolor(X1,Y1,mean(U_conv.^2,3)/Uj^2), shading interp, axis equal, colormap(custom_map1);
        ylim([-1.5 1.5]), xlim([0.3 5.48]), xticks([0.5 1 2 3 4 5]), c = colorbar; hold on;
        title(c,'$\overline{u_c.u_c}/{U_j}^2$','Interpreter','Latex','FontWeight','bold'), xlabel('X/D_e','FontWeight','bold'),ylabel('Y/D_e','FontWeight','bold');
        set(gca,'FontSize',textsize,'FontName','times');      caxis([0 0.015]);
        figure(108)
%       subplot(212)
        pcolor(X1,Y1,mean(V_conv.^2,3)/Uj^2), shading interp, axis equal, colormap(custom_map1);
        ylim([-1.5 1.5]), xlim([0.3 5.48]), xticks([0.5 1 2 3 4 5]), c = colorbar;
        title(c,'$\overline{v_c.v_c}/{U_j}^2$','Interpreter','Latex'), xlabel('X/D_e','FontWeight','bold'),ylabel('Y/D_e','FontWeight','bold');
        set(gca,'FontSize',textsize,'FontName','times');      caxis([0 0.01]);  hold on;
%% PLOTTING AVERAGE CONVECTIVE ENERGY ALONG SHEAR LAYERS
         Ux_convec = mean(U_conv.^2,3)/Uj^2;  Uy_convec = mean(V_conv.^2,3)/Uj^2;     
         figure(110)
         plot(X1,(Ux_convec(jet_edge(2),:)+Ux_convec(jet_edge(3),:))/2,'b','LineWidth',1.5,'DisplayName','V_{c \cdot Axial} - Lip line'), hold on;
         plot(X1,(Uy_convec(jet_edge(2),:)+Uy_convec(jet_edge(3),:))/2,'r','LineWidth',1.5,'DisplayName','V_{c \cdot Radial} - Lip line');
         plot(X1,Ux_convec(jet_edge(1),:),'-.b','LineWidth',1.5,'DisplayName','V_{c \cdot Axial} - Centerline'); title('Convective Velocity Profiles');
         plot(X1,Uy_convec(jet_edge(1),:),'-.r','LineWidth',1.5,'DisplayName','V_{c \cdot Radial} - Centerline');
         xlabel('\bf X/D_e'), ylabel('\bf \it(V_c.V_c)/{U_j}^2'), legend('show'), a = gca; l = legend;
         set(gca,'FontSize',textsize,'FontName','times'); grid on, box off, xlim([0 5.1]), ylim([-0.001 0.02]);
         a.YTickLabel = {'0','0.005','0.01','0.015','0.02'}; l.EdgeColor = [1 1 1]; l.Location = 'northwest'; set(gcf,'Position',[550 500 550 300]); 
     end
  end
%% VORTICITY CALCULATION
  prompt = [newline '>> PLOT VORTICITY?(y/n) - ']; vrtcy = input(prompt,'s');
  if strcmp(vrtcy,'y')
     prompt = [newline '>> ENTER NUMBER OF IMAGES - ']; imgs = input(prompt);
     Curl = zeros(size(U_inst,1),size(U_inst,2),imgs);   [X,Y] = meshgrid(X1,Y1);    
     for ctr = 1:imgs
       Curl(:,:,ctr) = curl(X,Y,U_inst(:,:,ctr),V_inst(:,:,ctr));
       figure(111)
       pcolor(X,Y,Curl(:,:,ctr)*Deq/(1000*Uj)),shading interp,hold on,colormap(vorticity_map);axis equal;xlim([0.3 5.48]), ylim([-2 2]);
       xlabel('X/D_e'), ylabel(yname); set(gca,'FontSize',textsize,'FontName','times'); 
       q = quiver(X,Y,U_inst(:,:,ctr),V_inst(:,:,ctr),'y');   q.Color = [0 0 0]; q.AutoScaleFactor = 3;
       pause(0.2)  
       drawnow
       frame = getframe(111);
       im{ctr} = frame2im(frame);
     end  
     cd(drive_in);       gifname = ['Vorticity_' condition{1} '.gif'];
        for idx = 1:size(im,2)
          [A,map] = rgb2ind(im{idx},256);
          if idx == 1
             imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.1);
          else
             imwrite(A,map,gifname,'WriteMode','append','DelayTime',0.1);
          end
        end
        clear im;       clear Curl;
        % AVERAGE VORTICITY
        commandwindow; prompt = [newline '>> PLOT AVERAGE VORTICITY?(y/n) - ']; chc = input(prompt,'s');
        if strcmp(chc,'y')
           Curl = zeros(size(U_inst,1),size(U_inst,2),size(U_inst,3));
            for ctr = 1:size(U_inst,3)  
              Curl(:,:,ctr) = curl(X,Y,U_inst(:,:,ctr),V_inst(:,:,ctr));
            end
            Avg_Curl = mean(Curl,3);
            figure(112)
            pcolor(X,Y,Avg_Curl*Deq/(1000*Uj)), shading interp,colormap(vorticity_map); hold on;axis equal;xlim([0.3 5.48]), ylim([-2 2]);
            xlabel('X/D_e'), ylabel(yname); set(gca,'FontSize',textsize,'FontName','times'); 
            q = quiver(X,Y,U_inst(:,:,ctr),V_inst(:,:,ctr),'y');   q.Color = [0 0 0]; q.AutoScaleFactor = 5;
        end
  end
%%  PLOT CONTOURS
  prompt = [newline '>> PLOT CONTOUR?(y/n) - ']; contr_plt = input(prompt,'s');
  if strcmp(contr_plt,'y')
     figure(500*n)
     pcolor(X1,Y1,Avg_U), shading interp, colormap jet,axis equal;   
     ylim([-2 2]), xlim([0.3 5.48]), xlabel('X/D_e'), ylabel(yname);xticks([0.5 1 2 3 4 5]);
     caxis([0 1]),cl = colorbar; title(cl,'U/U_j','FontWeight','bold'); set(gca,'FontSize',textsize,'FontName','times');
     figure(201*n)
     pcolor(X1,Y1,Avg_V), shading interp, colormap jet,axis equal;   
     ylim([-2 2]), xlim([0.3 5.48]), xlabel('X/D_e'), ylabel(yname);xticks([0.5 1 2 3 4 5]);
     caxis([-0.05 0.05]),cl = colorbar; title(cl,'U/U_j','FontWeight','bold');set(gca,'FontSize',textsize,'FontName','times');
  end
%%  IMAGE BINNING
  prompt = [newline '>> SHOW BINNED PROFILES(y/n) - ']; chc = input(prompt,'s');
  if strcmp(chc,'y') && strcmp(nozzle{m},'Minor')
     if NPR == 2.9 || NPR == 3.0
        [Phs1_U,Phs1_V,Phs2_U,Phs2_V,Phs1_Um,Phs1_Vm,Phs2_Um,Phs2_Vm] = Img_Binr(U_inst,V_inst,jet_edge);
%  PLOTTING PHASE 1
        figure(113); title('Avg Flapping Phase 1')
        pcolor(X,Y,Phs1_Um/Uj), shading interp, colormap jet, axis equal, xlim([0.3 5.48]), ylim([-2 2]);caxis([0 1]);
        xlabel('X/D_e');ylabel(yname), cl = colorbar; title(cl,'U/U_j','FontWeight','bold'); set(gca,'FontSize',textsize,'FontName','times');
%  PLOTTING PHASE 2
        figure(114); title('Avg Flapping Phase 2')
        pcolor(X,Y,Phs2_Um/Uj), shading interp, colormap jet, axis equal, xlim([0.3 5.48]), ylim([-2 2]);caxis([0 1]);
        xlabel('X/D_e');ylabel(yname), cl = colorbar; title(cl,'U/U_j','FontWeight','bold'); set(gca,'FontSize',textsize,'FontName','times');
        commandwindow;
%  VORTICITY CALCULATION
        prompt = [newline '>> PLOT VORTICITY METRICS?(y/n) - ']; chc = input(prompt,'s');
        if strcmp(chc,'y')
           Curl_Phs1 = curl(X,Y,Phs1_U(:,:,1),Phs1_V(:,:,1));     Curl_Phs2 = curl(X,Y,Phs2_U(:,:,1),Phs2_V(:,:,1));
           figure(115); pcolor(X,Y,Curl_Phs1*Deq/(1000*Uj)), shading interp, colormap(vorticity_map);
           xlabel('X/D_e'),ylabel(yname); hold on; title('Normalized Vorticity(\omega\cdotD_e/{U_j}) Phase - 1');
           set(gca,'FontSize',textsize,'FontName','times');
           figure(116); pcolor(X,Y,Curl_Phs2*Deq/(1000*Uj)), shading interp, colormap(vorticity_map);
           title('Normalized Vorticity(\omega\cdotD_e/{U_j}) Phase - 2'); xlabel('X/D_e'),ylabel(yname); hold on;
           set(gca,'FontSize',textsize,'FontName','times');
           prompt = [newline '>> PLOT VORTICITY VECTORS?(y/n) - '];   vec = input(prompt,'s');
           if strcmp(vec,'y')
              figure(115); q1 = quiver(X,Y,Phs1_U(:,:,1),Phs1_V(:,:,1),'y'); q1.Color = [0 0 0]; q1.AutoScaleFactor = 5; 
              figure(116); q2 = quiver(X,Y,Phs2_U(:,:,1),Phs2_V(:,:,1),'y'); q2.Color = [0 0 0]; q2.AutoScaleFactor = 5; 
           end
        end
     else
        disp('FOR SCREECH CONDITION ONLY');
     end
  end
%%  PLOT STREAMLINES
  prompt = [newline '>> PLOT STREAMLINES(y/n) - ']; chc = input(prompt,'s');
  if strcmp(chc,'y')
     [~,b] = find(X1<=5); [a,~] = find(Y1 > -1.2 & Y1 <= 1.2); X_s = X(a,b); Y_s = Y(a,b);
     U1 = Phs1_U(a,b,1); U2 = Phs2_U(a,b,1);  V1 = Phs1_V(a,b,1);  V2 = Phs2_V(a,b,1);
     figure(117); pcolor(X_s,Y_s,U1/Uj); shading interp, colormap jet, axis equal,caxis([0 1]); hold on;
     c = colorbar; title(c,'U/U_j','FontSize',textsize,'FontName','times','FontWeight','bold');
     hslice = streamslice(X_s,Y_s,[],U1,V1,[],X_s,Y_s,[],10); set(hslice,'Color',[0.6 0.6 0.6]);
     title('Streamlines Phase - 1'); xlabel('X/D_e'),ylabel(yname); hold on;
     set(gca,'FontSize',textsize,'FontName','times'); xlim([0.3 5]); ylim([-1.2 1.2]);
     figure(118); pcolor(X_s,Y_s,U2/Uj); shading interp, colormap jet, axis equal,caxis([0 1]); hold on;
     hslice = streamslice(X_s,Y_s,[],U2,V2,[],X_s,Y_s,[],10); set(hslice,'Color','w');
     title('Streamlines Phase - 2'); xlabel('X/D_e'),ylabel(yname); hold on;
     set(gca,'FontSize',textsize,'FontName','times'); xlim([0.3 5]); ylim([-1.2 1.2]);

  end
  %% PLOTTING THE SHOCK REFLECTION LINES ON CONVECTIVE VELOCITY PLOTS 
     if strcmp(nozzle{m},'Minor') && strcmp(rec,'y')
        if NPR == 2.9
           Sh_lc = [0.4358, 1.382, 2.329, 3.275, 4.157, 5.005];
        elseif NPR == 3
           Sh_lc = [0.4343, 1.437, 2.503, 3.538,4.443,5.413,6.189,6.933,7.741];
        end
        S = zeros(size(Y));
        for q = 1:length(Sh_lc)
            S(:) = Sh_lc(q);   figure(107);   plot(S,Y,'--k','LineWidth',1.5);  hold on;
            figure(108);   plot(S,Y,'--w','LineWidth',1.5);  hold on;
        end
     end

end
end
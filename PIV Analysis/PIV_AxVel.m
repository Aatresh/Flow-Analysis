 %% PLOTS THE CENTERLINE & RADIAL DISTRIBUTION OF AXIAL VELOCITIES
% clc;   clearvars;   fclose all;   set(0,'defaultfigurecolor',[1 1 1]);   
code = 'PIV_Utilization'; 
root1 = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\'; 
cd([root1 'Jet_Analysis\Global_Functions\']); addpath([root1 'Jet_Analysis\PIV_Stiching_Codes\V2\']);
colr = [0 0 0;0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.301 0.745 0.933] ;
root2 = 'X:\OneDrive - University of Cincinnati\Shared_Data_Sets\Stanford_Shared_Data\PIV\Stanford Data\Data Set - 2\';
tests = {'NPR_2p4_TR_1p0', 'NPR_2p5_TR_1p0', 'NPR_2p6_TR_1p0', 'NPR_2p9_TR_1p0', 'NPR_3p0_TR_1p0', 'NPR_3p6_TR_1p0', 'NPR_4p0_TR_1p0', 'NPR_5p0_TR_1p0',...
         'NPR_2p5_TR_1p1', 'NPR_3p0_TR_1p1', 'NPR_3p5_TR_1p1', 'NPR_3p6_TR_1p1', 'NPR_4p0_TR_1p1', 'NPR_4p5_TR_1p1',...
         'NPR_2p5_TR_2p0', 'NPR_3p0_TR_2p0', 'NPR_3p5_TR_2p0', 'NPR_3p6_TR_2p0', 'NPR_4p0_TR_2p0', 'NPR_4p5_TR_2p0',...
         'NPR_2p5_TR_2p6', 'NPR_3p0_TR_2p6', 'NPR_3p5_TR_2p6', 'NPR_3p6_TR_2p6', 'NPR_4p0_TR_2p6','NPR_4p5_TR_2p6'};
configTable;  a_amb = 347.0998;   load blckToGold.mat;  

condition = [tests(5)];   nozzle = 'Major';   config = 'S2';   NF = 'D';   writeTec = 'n';

%%  CONTOUR LIMITS
% Average Velocity & Text Size
limU = 1.17;         cUnorm = [0 limU];   textSize = 13;
% Turbulence & Reynold's Stres
cTKE = [8e-4 0.03];  cRxx = [0 0.03];     cRyy = [0 0.015];    cRxy = [-0.01 0.01]; 
%  DRIVE SELECTION 
[OutputStruct] = GF_DriveSelect(config,nozzle,code);  Deq = OutputStruct.dt;  nozzle = OutputStruct.nozzle; 
for n = 1:length(condition)
%  NORMALIZATION FACTORS & AXIS DEFINITION FOR LABELS
  NormFac.root1 = root1; NormFac.root2 = root2;  NormFac.typ = config(5:end);    NormFac.nozzle = nozzle;
  NormFac.NF = NF;       NormFac.De = Deq;       NormFac.op_cond = condition{n}; NormFac.code = code;
  driveIn = [OutputStruct.in_root condition{n}(9:14) '\' condition{n}(1:7) '\'];
%  FILE EXISTANCE CHECK
  if exist(driveIn)==0                                                              %#ok
     disp('->> Configuration needs to be computed <<-');
     break;
  end
%  PLOT TITLE
  [Mj,Uj,NPR,NTR] = GF_Velocity(condition{n});     name = ['NPR \thinspace ',num2str(NPR)];
%%  FILE NAMES
  [file_U,file_V,file_TKE,file_Rxx,file_Ryy,file_Rxy] = FileNamer(condition{n},nozzle);
%%  LOADING DATA FROM DAT FILES
  fid1 = fopen([driveIn file_U]);      fid2 = fopen([driveIn file_V]);     fid3 = fopen([driveIn file_TKE]); 
  fid4 = fopen([driveIn file_Rxx]);    fid5 = fopen([driveIn file_Ryy]);   fid6 = fopen([driveIn file_Rxy]); 
  data_U   = textscan(fid1,'%n%n%n');  data_V = textscan(fid2,'%n%n%n');   data_TKE = textscan(fid3,'%n%n%n');
  data_Rxx = textscan(fid4,'%n%n%n');  data_Ryy = textscan(fid5,'%n%n%n'); data_Rxy = textscan(fid6,'%n%n%n');
%  REARRANGING MATRICES FOR SORTING
  Dat1 = cat(2,data_U{1},data_U{2},data_U{3});         Dat2 = cat(2,data_V{1},data_V{2},data_V{3});
  Dat3 = cat(2,data_TKE{1},data_TKE{2},data_TKE{3});   Dat4 = cat(2,data_Rxx{1},data_Rxx{2},data_Rxx{3});
  Dat5 = cat(2,data_Ryy{1},data_Ryy{2},data_Ryy{3});   Dat6 = cat(2,data_Rxy{1},data_Rxy{2},data_Rxy{3});
%  MATRIX SORTING & RESIZING
  [X,Y,uBar] = GF_TrixSort(Dat1);   [~,~,vBar] = GF_TrixSort(Dat2);  [~,~,tkE] = GF_TrixSort(Dat3);
  [~,~,Rxx]  = GF_TrixSort(Dat4);   [~,~,Ryy]  = GF_TrixSort(Dat5);  [~,~,Rxy] = GF_TrixSort(Dat6);
%%  DIMENSION CHECK
  higherOrderStat.tkE = tkE;   higherOrderStat.Rxx = Rxx;   higherOrderStat.Ryy = Ryy;   higherOrderStat.Rxy = Rxy;
  [X,Y,uBar,vBar,higherOrderStat] = DimCheck(X,Y,Uj,uBar,vBar,higherOrderStat);
  tkE = higherOrderStat.tkE;   Rxx = higherOrderStat.Rxx;   Ryy = higherOrderStat.Ryy;   Rxy = higherOrderStat.Rxy;
  TKE2 = sqrt(abs(tkE*Uj^2))/Uj;
%%  AXIS DEFINITIONS & LIMITS
  [X,Y,xname,yname,lim_x,lim_y,Postn,sub] = GF_AxisNormTags(X,Y,NormFac,config);    %TKE = smoothdata(TKE,'gaussian',5);
  clear data_U data_V data_TKE;
%  FINDING THE CENTER OF THE JET & NORMALIZING WITH AMBIENT SPEED OF SOUND
 [jetEdge,Y] = GF_CenterFinder(uBar,Y,config,nozzle,'Ux');
 uNorm = (uBar)/Uj;   
%% CONTOUR - ACTUAL VELOCITY 
  figure(10*n);  
  aX = subplot(sub(1));  pcolor(X,Y,uBar);  shading interp;  caxis([0 Uj*limU]);  colormap(aX,blckToGold);      axis equal;  xlim(lim_x);  ylim(lim_y);  box off;          
  aX.TickLabelInterpreter = 'latex';        set(gca,'FontSize',textSize);         xlabel(xname,'Interpreter','latex');       ylabel(yname,'Interpreter','latex');
  title(['$|U|_{actual} - ',name,'$'],'Interpreter','latex'); c = colorbar;       c.TickLabelInterpreter = 'latex';          title(c,'${U_j}(m/s)$','Interpreter','latex');     
%% CONTOUR - RADIAL VELOCITY 
  aX = subplot(sub(2));  pcolor(X,Y,vBar);  shading interp;  caxis([-30 30]);     colormap(aX,'bluewhitered');  axis equal;  xlim(lim_x);  ylim(lim_y);  box off;        
  aX.TickLabelInterpreter = 'latex';        set(gca,'FontSize',textSize);         xlabel(xname,'Interpreter','latex');       ylabel(yname,'Interpreter','latex');
  title(['$|V|_{actual} - ',name,'$'],'Interpreter','latex'); c = colorbar;       c.TickLabelInterpreter = 'latex';          title(c,'$V_j(m/s)$','Interpreter','latex');
  set(gcf,'Position',Postn);
%% CONTOUR - VELOCITY NORMALIZED BY THEORETICAL EXIT VELOCITY
  figure(20*n);  
  aX = subplot(sub(1));  pcolor(X,Y,uNorm); shading interp;  caxis(cUnorm);  colormap(blckToGold);     axis equal;  xlim(lim_x);  ylim(lim_y);  box off;
  aX.TickLabelInterpreter = 'latex';        set(gca,'FontSize',textSize);    xlabel(xname,'Interpreter','latex');   ylabel(yname,'Interpreter','latex');
  title(['$|U|_{norm} - ',name,'$'],'Interpreter','latex');   c = colorbar;  c.TickLabelInterpreter = 'latex';      title(c,'$U/U_j$','Interpreter','latex');     
%% CONTOUR - TURBULENCE DISTRIBUTION  
  aX = subplot(sub(2));  pcolor(X,Y,tkE);   shading interp;  caxis(cTKE);    colormap(blckToGold);     axis equal;  xlim(lim_x);  ylim(lim_y);  box off;   
  aX.TickLabelInterpreter = 'latex';        set(gca,'FontSize',textSize);    xlabel(xname,'Interpreter','latex');   ylabel(yname,'Interpreter','latex');
  title(['$Turbulence - ',name,'$'],'Interpreter','latex');  c = colorbar;   c.TickLabelInterpreter = 'latex';      title(c,'$TKE/{U_j}^2$','Interpreter','latex');
  set(gcf,'Position',Postn);
%% CONTOUR - REYNOLD'S STRESSES(Rxx) DISTRIBUTION  
  figure(30*n);  
  aX = subplot(sub(1));  pcolor(X,Y,Rxx);   shading interp;  caxis(cRxx);    colormap(blckToGold);     axis equal;  xlim(lim_x);  ylim(lim_y);  box off;
  aX.TickLabelInterpreter = 'latex';        set(gca,'FontSize',textSize);    xlabel(xname,'Interpreter','latex');   ylabel(yname,'Interpreter','latex');
  title(['$Reynold''s Stress(R_{xx}) - ',name,'$'],'Interpreter','latex');   c = colorbar;                          c.TickLabelInterpreter = 'latex';
  title(c,'$R_{xx}/{U_j}^2$','Interpreter','latex');                         
%% CONTOUR - REYNOLD'S STRESSES(Ryy) DISTRIBUTION  
  aX = subplot(sub(2));  pcolor(X,Y,Ryy);   shading interp;   caxis(cRyy);   colormap(blckToGold);     axis equal;  xlim(lim_x);  ylim(lim_y);  box off;
  aX.TickLabelInterpreter = 'latex';        set(gca,'FontSize',textSize);    xlabel(xname,'Interpreter','latex');   ylabel(yname,'Interpreter','latex');
  title(['$Reynold''s Stress(R_{yy}) - ',name,'$'],'Interpreter','latex');   c = colorbar;                          c.TickLabelInterpreter = 'latex';
  title(c,'$R_{yy}/{U_j}^2$','Interpreter','latex');                         set(gcf,'Position',Postn);             
%% CONTOUR - REYNOLD'S STRESSES(Rxy) DISTRIBUTION  
  figure(50*n);aX = gca; pcolor(X,Y,Rxy);   shading interp;   caxis(cRxy);   colormap bluewhitered;    axis equal;  xlim(lim_x);  ylim(lim_y);  box off;
  aX.TickLabelInterpreter = 'latex';        set(gca,'FontSize',textSize);    xlabel(xname,'Interpreter','latex');   ylabel(yname,'Interpreter','latex');
  title(['$Reynold''s Stress(R_{xy}) - ',name,'$'],'Interpreter','latex');   c = colorbar;                          c.TickLabelInterpreter = 'latex';
  title(c,'$R_{xy}/{U_j}^2$','Interpreter','latex');                         %c.Title.Position = [0 100 0];
  if strcmp(config,'TS3')
     set(gcf,'Position',[25 80 1460 420]);
  else
     set(gcf,'Position',[105 305 886 382]);
  end
%% CONTOUR - LOCAL MACH NO.
To = 297.15;   gma = 1.4;   R = 287.058;
localMach = abs(sqrt( (uBar.^2)./ (gma*R*To - (uBar.^2)*((gma-1)/2) )));
figure(60*n);   aX = gca;  pcolor(X,Y,localMach); shading interp;   caxis([0 1.5]);  colormap(blckToGold);     axis equal;  xlim(lim_x);  ylim(lim_y);  box off;   hold on;
aX.TickLabelInterpreter = 'latex';   set(gca,'FontSize',textSize);    xlabel(xname,'Interpreter','latex');     ylabel(yname,'Interpreter','latex'); 
title(['$\bf M_{local} - ',name,'$'],'Interpreter','latex');   c = colorbar;  c.TickLabelInterpreter = 'latex';      title(c,'$U/U_j$','Interpreter','latex');    
set(gca,'Layer','top');  contour(X,Y,localMach,[1 1],':k','LineWidth',1.2);    set(gcf,'Position',[15 200 920 410]);
%%  PLOT - CENTERLINE LOCAL MACH NO.
figure(70);   aX = gca;    hold on;       aX.TickLabelInterpreter = 'latex';
if strcmp(nozzle,'Major') && strcmp(config(1),'T')
   plot(X,localMach(jetEdge(2),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor','k',      'DisplayName',[name ':Jet - 1']);
   plot(X,localMach(jetEdge(3),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor',colr(n,:),'DisplayName',[name ':Jet - 2']);
else
   plot(X,localMach(jetEdge(1),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor',colr(n,:),'DisplayName',['$',name,'$']);
end;  grid on;     xlim(lim_x);    ylim([0.6 1.6]);      c = legend;       xlabel(xname,'Interpreter','latex');   ylabel('$U/U_j$','Interpreter','latex'); 
set(gcf,'Position',[50 300 1050 480]);    set(gca,'FontSize',textSize);    c.EdgeColor = [1 1 1];                 c.Interpreter = 'latex';                                               

%% CHECK ALTERNATE POINTS FOR CENTERLINE
%% WRITE TO TECPLOT FORMAT
if strcmp(writeTec,'y')
   tecWrite(X,Y,flipud(uBar),driveIn,nozzle,condition{n},root1);
end
%% PLOT - NORMALIZED CENTERLINE VELOCITY
  figure(4);   aX = gca;    hold on;       aX.TickLabelInterpreter = 'latex';
  if strcmp(nozzle,'Major') && strcmp(config(1),'T')
     plot(X,uNorm(jetEdge(2),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor','k',      'DisplayName',[name ':Jet - 1']);
     plot(X,uNorm(jetEdge(3),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor',colr(n,:),'DisplayName',[name ':Jet - 2']);
  else
     plot(X,uNorm(jetEdge(1),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor',colr(n,:),'DisplayName',['$',name,'$']);
  end;  grid on;     xlim(lim_x);    ylim([0.5 1.3]);      c = legend;       xlabel(xname,'Interpreter','latex');   ylabel('$U/U_j$','Interpreter','latex'); 
  set(gcf,'Position',[50 300 1050 480]);    set(gca,'FontSize',textSize);    c.EdgeColor = [1 1 1];                 c.Interpreter = 'latex';                                               
%% PLOT - CENTERLINE VELOCITY(m/s)
  figure(5);   aX = gca;    hold on;       aX.TickLabelInterpreter = 'latex';
  if strcmp(nozzle,'Major') && strcmp(config(1),'T')
     plot(X,uBar(jetEdge(2),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor','k','DisplayName',[name ':Jet - 1']);
     plot(X,uBar(jetEdge(3),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor',colr(n,:),'DisplayName',[name ':Jet - 2']);
  else
     plot(X,uBar(jetEdge(1),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor',colr(n,:),'DisplayName',['$',name,'$']);
  end; grid on;     xlim(lim_x);    ylim([100 500]);      c = legend;        xlabel(xname,'Interpreter','latex');   ylabel('$U(m/s)$','Interpreter','latex'); 
  set(gcf,'Position',[50 300 1050 480]);    set(gca,'FontSize',textSize);    c.EdgeColor = [1 1 1];                 c.Interpreter = 'latex'; 
%% PLOT - LIPLINE VELOCITIES
  figure(6);   aX = gca;    hold on;       aX.TickLabelInterpreter = 'latex';
  if strcmp(nozzle,'Major') && strcmp(config(1),'T') && strcmp(config,'TR3')~=1
     plot(X,uNorm(jetEdge(4),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor','k','DisplayName',[name ':Jet - 1']);
     plot(X,uNorm(jetEdge(6),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor',colr(n,:),'DisplayName',[name ':Jet - 2']);
  else
     plot(X,uNorm(jetEdge(3),:),'-o','Color',colr(n,:),'LineWidth',1.2,'MarkerSize',3,'MarkerFaceColor',colr(n,:),'MarkerEdgeColor',colr(n,:),'DisplayName',['$',name,'$']);
  end; grid on;     xlim(lim_x);    ylim([0.3 1.01]);     c = legend;        xlabel(xname,'Interpreter','latex');   ylabel('$U/U_j$','Interpreter','latex'); 
  set(gcf,'Position',[50 300 1050 480]);    set(gca,'FontSize',textSize);    c.EdgeColor = [1 1 1];                 c.Interpreter = 'latex';
%% PLOT - RADIAL DISTRIBUTION of STATISTICS
   locs = [0.5, 1, 2.5, 5, 8];  mrkr = {'-s';'-v';'-o';'-^';'-p';'-h';'-*'};   figure(80*n);  set(gcf,'Position',[85 300 1270 490]);
   for ctr = 1:length(locs)
     [~,a1] = find(X<locs(ctr)); 
     if isempty(a1)
        continue;
     end;    axLoc = ['$Axial \thinspace Loc. - ',num2str(locs(ctr)),xname(4:end-1),'$'];
     % Radial distribution of axial velocity
     subplot(121);     plot(uNorm(:,a1(end)),Y,'Color',colr(ctr,:),'LineWidth',1.2,'DisplayName',axLoc);  
     xlim([0 1.2]);    ylim(lim_y);             grid on;      box off;         hold on;
     c = legend;       c.EdgeColor = [1 1 1];   set(gca,'FontSize',textSize);  c.Interpreter = 'latex';
     title('$Radial \thinspace Profile - Velocity(U_{norm,U_j})$','Interpreter','latex');  
     xlabel('$U/U_j$','Interpreter','latex');   ylabel(yname,'Interpreter','latex');
     % Radial distribution of turbulence
     subplot(122);     plot(smoothdata(tkE(:,a1(end)),'gaussian',1.5),Y,'Color',colr(ctr,:),'LineWidth',1.2,'DisplayName',axLoc);
     xlim([0 0.04]);   ylim(lim_y);             grid on;      box off;         hold on;
     c = legend;       c.EdgeColor = [1 1 1];   set(gca,'FontSize',textSize);  c.Interpreter = 'latex';
     title('$Radial \thinspace Profiles - Turbulence(TKE)$','Interpreter','latex');  
     xlabel('$TKE/{U_j}^2$','Interpreter','latex'); ylabel(yname,'Interpreter','latex'); 
   end; clear a1;  aX = subplot(121);  aX.TickLabelInterpreter = 'latex';    aX = subplot(122);  aX.TickLabelInterpreter = 'latex';
   % Saving jet edge location
   save([driveIn 'jetEdge'],'jetEdge');
end
%% FILE NAMING CONVENTION FOR INPUT DAT FILES
function [U_name,V_name,TKE_name,Rxx_name,Ryy_name,Rxy_name] = FileNamer(testCondtn,nozzle)
   U_name    = [testCondtn '_' nozzle '_AvgVx.DAT']; 
   V_name    = [testCondtn '_' nozzle '_AvgVy.DAT'];
   TKE_name  = [testCondtn '_' nozzle '_TurbKineticE.DAT'];
   Rxx_name  = [testCondtn '_' nozzle '_ReyStresXX.DAT'];
   Ryy_name  = [testCondtn '_' nozzle '_ReyStresYY.DAT'];
   Rxy_name  = [testCondtn '_' nozzle '_ReyStresXY.DAT'];
end
%% CHECK IF DIMENSIONS OF ALL THE DATA MATRICES ARE EQUAL
function [X,Y,uBar,vBar,higherOrderStat] = DimCheck(X,Y,Uj,uBar,vBar,higherOrderStat)
% Row & Column counts for each statistic
Rows = zeros(6,1);   Cols = zeros(6,1);
% Read all row sizes
Rows(1) = size(uBar,1);                         Rows(2) = size(vBar,1); 
Rows(3) = size(higherOrderStat.tkE,1);          Rows(4) = size(higherOrderStat.Rxx,1);          
Rows(5) = size(higherOrderStat.Ryy,1);          Rows(6) = size(higherOrderStat.Rxy,1);
% Read all column sizes
Cols(1) = size(uBar,2);                         Cols(2) = size(vBar,2); 
Cols(3) = size(higherOrderStat.tkE,2);          Cols(4) = size(higherOrderStat.Rxx,2);
Cols(5) = size(higherOrderStat.Ryy,2);          Cols(6) = size(higherOrderStat.Rxy,2);
% Finding Min. and Max. number of rows & columns
rowMax = max(Rows);  colMax = max(Cols);
rowMin = min(Rows);  colMin = min(Cols);
% Normalizinf Uj sq to plot DaVis results
if max(max(higherOrderStat.tkE)) > 1
   higherOrderStat.tkE = higherOrderStat.tkE/Uj^2;    higherOrderStat.Rxx = higherOrderStat.Rxx/Uj^2;
   higherOrderStat.Ryy = higherOrderStat.Ryy/Uj^2;    higherOrderStat.Rxy = higherOrderStat.Rxy/Uj^2;
end
if isequal(rowMin/rowMax,colMin/colMax)
   return;
else
% Assigning Min value to Row & Column to data sets
X = X(1:colMin);      Y = Y(1:rowMin);
uBar = uBar(1:rowMin,1:colMin);                 vBar = vBar(1:rowMin,1:colMin);
higherOrderStat.tkE = higherOrderStat.tkE(1:rowMin,1:colMin);
higherOrderStat.Rxx = higherOrderStat.Rxx(1:rowMin,1:colMin);
higherOrderStat.Ryy = higherOrderStat.Ryy(1:rowMin,1:colMin);
higherOrderStat.Rxy = higherOrderStat.Rxy(1:rowMin,1:colMin);
end
end
%% CONVERT TO TECPLOT FORMAT 
function tecWrite(X,Y,uBar,driveOut,nozzle,condition,root1)
addpath([root1 'Jet_Analysis\Global_Functions\']);      cd(driveOut);    
[X1,Y1] = meshgrid(X,Y);
% Holding Struct to write to Tecplot format
holdMat = [];    holdMat.Nvar = 4;              % Number of variables to store
holdMat.varnames = {'x','y','z','Vx'};          % Variable Names
holdMat.surfaces(1).zonename = 'Jet-1';         % Contour Zone
holdMat.surfaces(1).x = X1;                     % Axial Co-ordinate
holdMat.surfaces(1).y = Y1;                     % Radial Co-ordinate
% holdMat.surfaces(1).order = 2;                % To plot XZ plot
holdMat.surfaces(1).v(1,:,:) = uBar;            % Contour Values(uBar,vBar,etc.)
contourName = ['Jet-1_' nozzle '_' condition '.plt'];
% Function to export Tecplot format
mat2tecplot(holdMat,contourName);
end
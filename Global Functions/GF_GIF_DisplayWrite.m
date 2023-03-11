function GF_GIF_DisplayWrite(FigMatrix,gifInput)
%  DESCRIPTION:FUNCTION TO DISPLAY & SAVE .GIF FILES
%  Accepts the image matrix as part of an input Struct file & displays the .gif file 
%  in an image window & saves it to a specifed output folder
%  Parameters of Struct(iNpStruct) file:
%  nozzle   - Orientation of nozzle   | NF        - axis normalization factor
%  Xn,Yn    - Axes                    | TYP       - Type of data used for .gif
%  bckLim   - Contour limits          | Mode      - SPOD mode(if included)
%  theta    - Phase angle(if included)| drive_out - Output folder
%  gif_name - Name of output .gif     | frameDt   - Time between frames for recording
%  pausTimer- Time between frames for display
nozzle = gifInput.nozzle;  Xn = gifInput.Xn;   Yn = gifInput.Yn;    pausTimer = 0.01;   
if strcmp(gifInput.NF,'H')
   if strcmp(nozzle,'Minor') || strcmp(nozzle,'Sym_Lin') 
      yName = '$Y/h$';    
   else
      yName = '$Z/h$';    
   end;    xName = '$X/h$';
else
   if strcmp(nozzle,'Minor') || strcmp(nozzle,'Circular')
      yName = '$Y/D_e$';  
   else
      yName = '$Z/D_e$';   
   end;  xName = '$X/D_e$';
end;     limX = gifInput.lim_x; limY = gifInput.lim_y;
%% GIF GENERATION
%  SCHLIEREN .gif
if strcmp(gifInput.TYP,'SC')
   figNo = 1000; figure(figNo);  set(gcf,'Position',gifInput.Sz); frameDt = 1/15;
   disp([newline,'->> Frame Rate - ',num2str(1/frameDt),' fps']);
   if strcmp(gifInput.config,'TS2') %|| strcmp(INP_STRUCT.config,'TR2')
      clim = [10 400];
   else
      if isempty(gifInput.bckLim) ~= 1
         clim = gifInput.bckLim;
      else 
         clim = [10 900];
      end
   end
  for ctr = 1:100
      pcolor(Xn,Yn,squeeze(FigMatrix(ctr,:,:)));   shading interp,   axis equal;    colormap gray; 
      xlabel(xName,'Interpreter','latex');         ylabel(yName,'Interpreter','latex'); 
      set(gca,'FontSize',13);    caxis(clim);      xlim(limX),      ylim(limY);   ax = gca;
      ax.TickLabelInterpreter = 'latex';           xticks([0 2 4 6 8 10 12]);              
      pause(pausTimer);   drawnow;   curntFrame = getframe(figNo);
      gifHold{ctr} = frame2im(curntFrame);                                 %#ok
  end
%  POD RECONSTRUCTION .gif
elseif strcmp(gifInput.TYP,'POD')
   figNo = 1000;    figure(figNo);     set(gcf,'Position',gifInput.Sz);    frameDt = 1/2;
   disp([newline,'->> Frame Rate - ',num2str(1/frameDt),' fps']);           load schJet.mat;      %#ok
   if strcmp(gifInput.config,'TS2') %|| strcmp(INP_STRUCT.config,'TR2')
      clim = [10 400];
   else
      if isempty(gifInput.bckLim) ~= 1
         clim = gifInput.bckLim;
      else 
         clim = [10 900];
      end
   end
  for ctr = 1:100
      pcolor(Xn,Yn,squeeze(FigMatrix(ctr,:,:)));   shading interp,   axis equal;    colormap(gray); 
      xlabel(xName,'Interpreter','latex');         ylabel(yName,'Interpreter','latex'); 
      set(gca,'FontSize',13);    caxis(clim);      xlim(limX),      ylim(limY);   ax = gca;
      ax.TickLabelInterpreter = 'latex';           xticks([0 2 4 6 8]);              
      pause(pausTimer);   drawnow;   curntFrame = getframe(figNo);
      gifHold{ctr} = frame2im(curntFrame);                                 %#ok
  end
%  SPOD .gif
elseif strcmp(gifInput.TYP,'SPOD')
  figNo = 1000; figure(figNo);  nt = gifInput.nt;     time = gifInput.time;   Mode = gifInput.Mode;   L = gifInput.L;   frameDt = 0.05;
  freqNo   = gifInput.freqNo;  freq = gifInput.Freq; gifMap = gifInput.map;     clim = gifInput.clim;   St = gifInput.St;
  for ctr = 1:nt
      pcolor(Xn,Yn,real(squeeze(FigMatrix(freqNo,:,:,Mode)*exp(2i*pi*freq(freqNo)*time(ctr)))));
      shading interp,  axis equal,  colormap(gifMap),  caxis(clim);  set(gcf,'Position',gifInput.Sz);
      xlim(limX),     ylim(limY); xlabel(xName,'Interpreter','latex');  ylabel(yName,'Interpreter','latex'); yticks([-4 -3 -2 -1 0 1 2 3 4]);
      title(['$St \approx ' num2str(round(St(freqNo),2)) '; Mode = ',num2str(Mode) '; \lambda = ' num2str(L(freqNo,Mode),'%.4g$')],'Interpreter','latex');
      set(gca,'FontSize',12,'FontName','times');    hold on, ax = gca; ax.TickLabelInterpreter = 'latex'; 
      pause(pausTimer);  drawnow;     curntFrame = getframe(figure(figNo));
      gifHold{ctr} = frame2im(curntFrame);  hold off;                     %#ok
  end
%  DMD .gif
elseif strcmp(gifInput.TYP,'DMD')
  figNo = 1000;  figure(figNo);  mode = gifInput.mode;  theta = gifInput.theta;  frameDt = 0.05;
  for ctr = 1:length(theta)
      pcolor(Xn,Yn,squeeze(real(FigMatrix(mode,:,:,1)*exp(1i*theta(ctr)))));
      shading interp;  colormap(gifInput.dmd_Colormap);  axis equal;   caxis([-0.004 0.004])
      xlabel(xName,'Interpreter','latex');  ylabel(yName,'Interpreter','latex'); title(gifInput.gifTitle,'Interpreter','latex'); 
      set(gca,'FontSize',12);  xlim(limX);  ylim(limY); set(gca,'TickLabelInterpreter','latex');
      pause(pausTimer);  drawnow;    curntFrame = getframe(figure(figNo));
      gifHold{ctr} = frame2im(curntFrame);  hold off                      %#ok    
  end
%  Schlieren Phase .gif
elseif strcmp(gifInput.TYP,'SchFourier')
       if strcmp(gifInput.config,'C')
          limX = [-1 6.5];     limY = [-2.2 2.2];   xFil = [0 -1 -1 0];       yFil = [0.5 0.82 -0.82 -0.5];
       elseif strcmp(gifInput.config,'S') && strcmp(gifInput.nozzle,'Major')
          limX = [-0.5 6.5];   limY = [-2.2 2.2];   xFil = [0 -0.5 -0.5 0];   yFil = [0.63 0.79 -0.79 -0.63]; 
       elseif strcmp(gifInput.config,'S') && strcmp(gifInput.nozzle,'Minor')
          limX = [-0.5 6.5];   limY = [-2.2 2.2];   xFil = [0 -0.5 -0.5 0];   yFil = [0.35 0.53 -0.53 -0.35];
       end
    figNo = 1000;   figure(figNo);   nt = gifInput.nt;  time = gifInput.time;  set(gcf,'Position',gifInput.Sz); 
    freqNo = gifInput.freqNo;       frameDt = 0.05;     freq = gifInput.freq;  clim = gifInput.clim;
    for ctr = 1:nt
        pcolor(Xn,Yn,-angle(squeeze(FigMatrix(:,:,freqNo)*exp(2i*pi*freq(freqNo)*time(ctr)))));   shading interp;   axis equal;
        xlabel(xName,'Interpreter','latex');    xlim(limX);   caxis(clim);    colormap(gifInput.cMap);
        ylabel(yName,'Interpreter','latex');    ylim(limY);   c = colorbar;   c.Ticks = [-pi,-pi/2,0,pi/2,pi];
        c.TickLabels = {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'};         c.TickLabelInterpreter = 'latex';
        title(['$\Delta \phi \thinspace Freq. =',num2str(freq(freqNo-1)),'Hz$'],'Interpreter','latex');
        title(c,'$\bf \Delta\phi$','Interpreter','latex');   set(gca,'FontSize',12,'TickLabelInterpreter','latex');   hold on;   
        set(gca,'Layer','top');   fill(xFil,yFil,'k');       ax = gca;   ax.TickLength = [0 0];   box off;
        pause(pausTimer);   drawnow;   curntFrame = getframe(figure(figNo));
        gifHold{ctr} = frame2im(curntFrame);   hold off;                  %#ok
    end
%  Schlieren Convection Separation .gif
elseif strcmp(gifInput.TYP,'CONVC-SEP')
    figNo = 1000;   figure(figNo);   set(gcf,'Position',gifInput.Sz);  nt = gifInput.nt;  frameDt = gifInput.frameDt;
    for ctr = 1:nt
    % Instantaneous
        subplot(311);  pcolor(Xn,Yn,real(squeeze(FigMatrix.Fluc(:,:,ctr))));      shading interp;  ax1 = gca;  axis equal;  hold on;
        caxis(gifInput.cLI); colormap(ax1,'gray');  xlim(limX);   ylim(limY);  set(gca,'FontSize',13,'TickLabelInterpreter','latex');
        xlabel(xName,'Interpreter','latex');   ylabel(yName,'Interpreter','latex');  title(gifInput.Ftitle,'Interpreter','latex');
    % Downstream Components
        subplot(312);   pcolor(Xn,Yn,real(squeeze(FigMatrix.DwnStrm(:,:,ctr))));  shading interp;  ax2 = gca;  axis equal;  hold on;
        caxis(gifInput.cLC);  colormap(ax2,gifInput.cMap);  xlim(limX);   ylim(limY);  set(gca,'FontSize',13,'TickLabelInterpreter','latex');
        xlabel(xName,'Interpreter','latex');   ylabel(yName,'Interpreter','latex');  title(gifInput.Dtitle,'Interpreter','latex');
    % Upstream Components
        subplot(313);   pcolor(Xn,Yn,real(squeeze(FigMatrix.UpStrm(:,:,ctr))));   shading interp;  ax3 = gca;  axis equal;  hold on;
        caxis(gifInput.cLC);  colormap(ax3,gifInput.cMap);  xlim(limX);   ylim(limY);  set(gca,'FontSize',13,'TickLabelInterpreter','latex');
        xlabel(xName,'Interpreter','latex');   ylabel(yName,'Interpreter','latex');  title(gifInput.Utitle,'Interpreter','latex');
        pause(pausTimer);  drawnow;    curntFrame = getframe(figure(figNo));   
        gifHold{ctr} = frame2im(curntFrame);  hold off;                  %#ok     
    end
end
%%  GIF SAVING TO DISK
prompt2 = [newline '->> Save .gif?(y/n) - '];    chc = input(prompt2,'s');   gifName = gifInput.gifName;
if strcmp(chc,'y')
   cd(gifInput.driveOut); 
   for idx = 1:size(gifHold,2)
     [gifVals,gifMap] = rgb2ind(gifHold{idx},256);
     if idx == 1
        imwrite(gifVals,gifMap,[gifName,'.gif'],'gif','LoopCount',Inf,'DelayTime',frameDt);
     else
        imwrite(gifVals,gifMap,[gifName,'.gif'],'WriteMode','append','DelayTime',frameDt);
     end
   end
end
%% CHANGING DIRECTORY BACK TO GLOBAL FUNCTIONS
cd('X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\Jet_Analysis\Global_Functions\');


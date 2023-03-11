function GF_MovieMaker(movieInput,FigMatrix)
%% DESCRIPTION: Generates a .MP4 video by taking a series of images as input
cd(movieInput.driveOut);      addpath(movieInput.driveCode);
% Frame Rate                  Video Object
fRate = 20;                   vw = VideoWriter(movieInput.videoName,'MPEG-4');
vw.FrameRate = fRate;         open(vw);
if isfield(movieInput,'sX') && isfield(movieInput,'sY')
   sourceLocX = movieInput.sX; sourceLocY = movieInput.sY;
   plotSource = movieInput.source;
end
%% Schlieren Phase
%  Plot figure & figure formatting
if strcmp(movieInput.code,'SchFourier')
   figNo = 1000;               figure(figNo);                 set(gcf,'Position',movieInput.figSiz);    cLim  = movieInput.cLim;   
   time  = movieInput.time;    pausTimer = 1/fRate;           StNo = movieInput.StNo;               %St    = movieInput.St;   
   Xn    = movieInput.Xn;      Yn        = movieInput.Yn;     nt   = movieInput.nt;                 xName = movieInput.xName;
   yName = movieInput.yName;   freq  = movieInput.freq;       
   if strcmp(movieInput.config,'C') 
      if strcmp(movieInput.style,'hf')
         limY = [0 2.2];      txtOneA = 1.1;       txtOneB = -1.3;     txtTwoA = 1.7;     txtTwoB = -1.7;
      else
         limY = [-2.2 2.2];   txtOneA = 1.1;       txtOneB = -3.3;     txtTwoA = 1.7;     txtTwoB = -3.7;  
      end
      limX = [-1 6.5];     xFil = [0 -1 -1 0];     yFil = [0.5 0.82 -0.82 -0.5];
   elseif strcmp(movieInput.config,'S') 
      limX = [-0.5 6.5];   limY = [-2.2 2.2];      txtOneA = 1.5;      txtOneB = -3.3;    txtTwoA = 2;       txtTwoB = -3.7;
      if strcmp(movieInput.nozzle,'Major')
         xFil = [0 -0.5 -0.5 0];   yFil = [0.63 0.79 -0.79 -0.63];
      elseif strcmp(movieInput.nozzle,'Minor')
         xFil = [0 -0.5 -0.5 0];   yFil = [0.35 0.53 -0.53 -0.35];
      end
   end
   for ctr = 1:nt
       pcolor(Xn,Yn,-angle(squeeze(FigMatrix(:,:,StNo)*exp(2i*pi*freq(StNo)*time(ctr)))));   shading interp;   axis equal;       aX = gca;
       xlabel(xName,'Interpreter','latex');    xlim(limX);     caxis(cLim);     colormap(movieInput.schJet);   c = colorbar;     aX.YTick = [-2 -1 0 1 2];      
       ylabel(yName,'Interpreter','latex');    ylim(limY);     c.Ticks = [-pi,-pi/2,0,pi/2,pi];                c.TickLabels = {'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'};
       c.TickLabelInterpreter = 'latex';       set(gca,'FontSize',13,'TickLabelInterpreter','latex');          hold on; 
       set(gca,'Layer','top');                 fill(xFil,yFil,'k');   ax = gca;   ax.TickLength = [0 0];       box off;   
       title(movieInput.fTitle,'Interpreter','latex');   
       text(txtOneA,txtOneB,movieInput.descrpText1, 'Color','k','FontSize',12,'Interpreter','latex');
       text(txtTwoA,txtTwoB,movieInput.descrpText2, 'Color','k','FontSize',12,'Interpreter','latex'); 
       if isfield(movieInput,'sX') && isfield(movieInput,'sY')
          sourcePoint(sourceLocX,sourceLocY,plotSource);
       end
       pause(pausTimer);                      drawnow;         curntFrame = getframe(figure(figNo));      
       writeVideo(vw,curntFrame);             hold off;                   
   end
end
% Close video file & end movie
close(vw);                  cd(movieInput.driveCode);
end
function sourcePoint(sourceLocX,sourceLocY,plotSource)
if strcmp(plotSource,'on')
   for ctr = 1:length(sourceLocX)
       set(gca,'Layer','top');   plot(sourceLocX(ctr),sourceLocY(ctr),'h','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','k');
       ax = gca;   ax.TickLength = [0 0];
   end
end
end
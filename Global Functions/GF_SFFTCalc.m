function [outputStruct] = GF_SFFTCalc(funcInput)
% DESCRIPTION: SPATIAL FOURIER TRANSFORM ANALYSIS   
% Function computes the spatial Fourier tranform of input data to convert
% spatial information into wavenumber spectra. The input matrix can be
% either of the following:
% 1. Coherent Fluctuation matrix computed from POD(Schlieren or PIV)
% 2. Energy & Phase matrix obtained from Spatial Fourier analysis at a specific frequency
% 3. Energy & Phase matrix obtained from SPOD analysis at a specific frequency
% ---------------------------------------------------------------------------
% Inputs:  Struct with Spatial Axes, Spatial Energy distribution; Peak Frequency(Hz); Length Scales(mm)
% Outputs: Spatial wavenumber vs Radial Energy distribution contour; Wave Nos.; FFT Spectrum
%---------------------------------------------------------------------------
    Xn = funcInput.Xn;   xName = funcInput.xName;   spatMat = funcInput.inpMat;        aAmb = 345.501931;        NF = funcInput.NF;             
    Yn = funcInput.Yn;   yName = funcInput.yName;   lenScales = funcInput.lenScales;   normFac = lenScales(1);   nozHt = lenScales(2);
    scrFreq = funcInput.scrFreq;  
%  Defining nozzle lipline & axis labels based on Configuration & NF check
%  Final axial limits to control spatial domain    
    prompt = [newline,'==> Final Axial  Location(0-',num2str(Xn(end)),'): '];    xEnd = input(prompt);  
    prompt = [newline,'==> Limit Radial Location(+/-',num2str(Yn(end)),'): '];   yLim = input(prompt);
    [~,xStrt] = find(Xn==0);       [~,xEnd] = find(Xn<xEnd);        xEnd = xEnd(end)+1;   schJet = [];     
    [yLimUp,~] = find(Yn<yLim);    [yLimDw,~] = find(Yn>-1*yLim);   radSpan = yLimDw(1)-1:yLimUp(end)+1;       
    yNew = Yn(radSpan);      fD = pwd;    load([fD,'\','schJet.mat']);     %#ok                   
%  Making the number of axial locations odd numbered 
    if rem(size(xStrt:xEnd,2),2)
       axSpan = xStrt:xEnd;          % Odd numbered
    else
       axSpan = xStrt:xEnd - 1;      % Changing Even numbered span to Odd
    end;   xEnd = Xn(axSpan(end));   % Re-assigning final axial valu to variable
%  Re-defining input matrix to match axial limits
    spatMat = spatMat(radSpan,axSpan);
%  Defining the limits of the axial domain(in mm), the number of points in the  domain & the domain resolution
    axDom = Xn(axSpan)*normFac;   domainPts = size(axDom,2);   dx = (axDom(3) - axDom(2))/1000;   xNew = Xn(axSpan);  
%  Wavenumber Range computation - calls another function
%  Computing the wavenumber range, finding the negative & positive wavenumber locations & rearranging to be zero centered
    kAct = GF_FourierFreq(domainPts,2*pi/dx);   kNos = fftshift(kAct);   kRangeD = kNos*nozHt/1000;   zLoc = find(kRangeD==0);
%  Computing spatial FFT along the streamwise direction for the entire domain & scaling abs with max([0 1] scale)
    inpSig = spatMat - mean(spatMat,2);    sigFFT = fft(inpSig,axSpan(end),2);   baseFFT = abs(sigFFT);     %#ok
%  Plotting the Base Spatial Fourier Result
%     fig1 = contourPlot(baseFFT,kNos,yNew,'$k$',yName,'auto',schJet,'$Base \thinspace Spatial \thinspace FFT$','$|\hat{\psi}|$');   axis tight;
%  Defining wave number axis label based on nozzle configuration
   if lenScales(1) == nozHt && strcmp(NF,'D')
      kName = '$kD_e$';
   else
      kName = '$kh_j$';   
   end
%  Transforming(shift & flip) & plotting the Spatial Fourier result    
    sigFFT1 = fftshift(sigFFT,2);   sigFFT1 = fliplr(sigFFT1);   sigFFT1 = abs(sigFFT1);   sigFFT2 = sigFFT1/(max(max(sigFFT1)));  sigFFT2 = smoothdata(sigFFT2,'gaussian',10);   
    fig2 = contourPlot(sigFFT2,kRangeD,yNew,kName,yName,[-18 18],schJet,['$f_s=',num2str(scrFreq),'Hz$'],'$|\hat{\psi}|$');   
%  Plotting log10 of abs with scaling shifted to positive values & [0 1] range
    sigNew = log10(sigFFT1);   sMin = min(min(sigNew));   
    if sMin == -Inf
       sMin = 0;
    end
    sigNew = sigNew + abs(sMin);   sigNew = sigNew/max(max(sigNew));                                        %#ok
%     fig3 = contourPlot(sigNew,kRangeD,yNew,kName,yName,[-20 20],schJet,'$Spatial \thinspace FFT:log_{10} \thinspace Scaled$','$log_{10}(|\hat{\psi}|)$');
%  Separating & plotting upstream & downstream components
    S1 = sigFFT;   S1(:,1:zLoc)   = 0;   dwStrmRes = S1;   clear S1;   
    S1 = sigFFT;   S1(:,zLoc:end) = 0;   upStrmRes = S1;   clear S1;
%     fig4 = contourPlot(abs(fliplr(fftshift(dwStrmRes,2))),kRangeD,yNew,kName,yName,[-20 20],schJet,'$Downstream \thinspace Component$','$|\hat{\psi}|$');
%     fig5 = contourPlot(abs(fliplr(fftshift(upStrmRes,2))),kRangeD,yNew,kName,yName,[-20 20],schJet,'$Upstream \thinspace Component$'  ,'$|\hat{\psi}|$');
%  Inverting the seperated components & plotting 
    dwStrm = ifft(dwStrmRes,size(dwStrmRes,2),2);     upStrm = ifft(upStrmRes,size(upStrmRes,2),2);
    fig6 = contourPlot(real(dwStrm),xNew,yNew,xName,yName,[0 xEnd],schJet,'$Downstream \thinspace Component$','$R_e(inv(\hat{\psi_d}))$');   
    figure(fig6);   axis equal;   xlim([0 xEnd]);   set(gcf,'Position',[70 400 715 400]);     
    fig7 = contourPlot(real(upStrm),xNew,yNew,xName,yName,[0 xEnd],schJet,'$Upstream \thinspace Component$'  ,'$R_e(inv(\hat{\psi_u}))$');
    figure(fig7);   axis equal;   xlim([0 xEnd]);   set(gcf,'Position',[70 400 715 400]);     
%  Computing sonic wave number(kA) - representing the speed of sound - for upstream & downstream
    kAdw = zeros(size(Yn));   kAdw(:) = ((2*pi*scrFreq)/aAmb)*nozHt/1000;   kAup = kAdw*-1;
%  Plotting kA line on contours
%     figure(fig1);   hold on;  set(gca,'Layer','top');   plot(kAup,Yn,'-w','LineWidth',1.2);   plot(kAdw,Yn,'-w','LineWidth',1.2);
    figure(fig2);   hold on;  set(gca,'Layer','top');   plot(kAup,Yn,'-w','LineWidth',1.2);   plot(kAdw,Yn,'-w','LineWidth',1.2);
%     figure(fig3);   hold on;  set(gca,'Layer','top');   plot(kAup,Yn,'-w','LineWidth',1.2);   plot(kAdw,Yn,'-w','LineWidth',1.2);
%     figure(fig4);   hold on;  set(gca,'Layer','top');   plot(kAdw,Yn,'-w','LineWidth',1.2);
%     figure(fig5);   hold on;  set(gca,'Layer','top');   plot(kAup,Yn,'-w','LineWidth',1.2);   
%  Outputs - fftRes is the zer-centered FFT result
    outputStruct.dwStrmRes = dwStrmRes;   outputStruct.kD   = kRangeD;   outputStruct.baseFFT = sigFFT;   outputStruct.xNew = xNew;      
    outputStruct.upStrmRes = upStrmRes;   outputStruct.kAct = kAct;      outputStruct.figNo   = fig2;     outputStruct.yNew = yNew;
    
%%  Function to plot the contours
    function figNo = contourPlot(inpMat,xAxis,yAxis,xName,yName,xLim,contourMap,contourTitle,cBarTitle)
        GF_FigurePlot(inpMat,xAxis,yAxis);    xlabel(xName,'Interpreter','latex');   axis normal;   ax = gca;    set(gcf,'Position',[35 300 900 550]);    
        ax.TickLabelInterpreter = 'latex';    ylabel(yName,'Interpreter','latex');   xlim(xLim);    ax.FontSize = 12;   curntFig = gcf;   c = colorbar;
        c.TickLabelInterpreter = 'latex';     cL = title(c,cBarTitle,'Interpreter','latex');   cL.Rotation = 90;   cL.Position = [55 170 0];  colormap(contourMap);  
        title(contourTitle,'Interpreter','latex');   cL.FontSize = 13;   figNo = curntFig.Number;
    end
end
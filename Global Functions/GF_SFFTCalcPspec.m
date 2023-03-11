function [outputStruct] = GF_SFFTCalcPspec(funcInput)
% DESCRIPTION: SPATIAL FOURIER TRANSFORM ANALYSIS USING psectrum  
% Function computes the spatial Fourier tranform of input data to convert spatial information into
% wavenumber spectra using the spectrum function. The input matrix can be either of the following:
% 1. Coherent Fluctuation matrix computed from POD(Schlieren or PIV)
% 2. Energy & Phase matrix obtained from Spatial Fourier analysis at a specific frequency
% 3. Energy & Phase matrix obtained from SPOD analysis at a specific frequency
% ---------------------------------------------------------------------------
% Inputs:  Struct with Spatial Axes, Spatial Energy distribution; Peak Frequency(Hz); Length Scales(mm)
% Outputs: Spatial wavenumber vs Radial Energy distribution contour; Wave Nos.; FFT Spectrum
%---------------------------------------------------------------------------

    Xn = funcInput.Xn;   xName = funcInput.xName;   spatMat   = funcInput.inpMat;      NF    = funcInput.NF;               
    Yn = funcInput.Yn;   yName = funcInput.yName;   lenScales = funcInput.lenScales;   nozHt = lenScales(2);      
    scrFreq = funcInput.scrFreq;                    aAmb = 345.501931; 
%  Defining wave number axis label based on nozzle configuration
    if lenScales(1) == nozHt && strcmp(NF,'D')
       kName = '$kD_e$';
    else
       kName = '$kh_j$';   
    end
%  Defining nozzle lipline & axis labels based on Configuration & NF check
%  Final axial limits to control spatial domain    
    prompt = [newline,'==> Final Axial  Location(0-',num2str(Xn(end)),'): '];    xEnd = input(prompt);  
    prompt = [newline,'==> Limit Radial Location(+/-',num2str(Yn(end)),'): '];   yLim = input(prompt);
    [~,xStrt] = find(Xn==0);       [~,xEnd] = find(Xn<xEnd);        xEnd = xEnd(end)+1;     
    [yLimUp,~] = find(Yn<yLim);    [yLimDw,~] = find(Yn>-1*yLim);   radSpan = yLimDw(1)-1:yLimUp(end)+1;       
    yNew = Yn(radSpan);            schJet = [];    fD = pwd;    load([fD,'\','schJet.mat']);                   %#ok                   
%  Making the number of axial locations odd numbered 
    if rem(size(xStrt:xEnd,2),2)
       axSpan = xStrt:xEnd;          % Odd numbered
    else
       axSpan = xStrt:xEnd - 1;      % Changing Even numbered span to Odd
    end 
 % Re-assigning final axial value to variable & defining SFF matrix
    xEnd = Xn(axSpan(end));   colNos = 1:size(axSpan,2);    pSpec = [];
%  Re-sizing matrix to match axial limits & defining Spatial Frequency: pixels/m
    spatMat = spatMat(radSpan,axSpan);   pixLen = size(axSpan,2);   spaLen = xEnd*lenScales(1);   spaFreq = pixLen/spaLen;   
%  Computing the spatial FFT using the pspectrum function
    for ctr = 1:size(radSpan,2)
        [P,~] = pspectrum(spatMat(ctr,:),spaFreq,'Leakage',1,'TwoSided',true); 
        pSpec(ctr,:) = P';                                                                                     %#ok
    end
%  Spatial resolution, wavenumbers based on length of pSpec & correction factor for correct 
%  wavenumber scaling & fftshift to center 0 in spatial frequency
    newPixLen = size(pSpec,2);     dxNew = spaLen/newPixLen;        kPspec = GF_FourierFreq(newPixLen,2*pi/dxNew);   
    corctnFac = pixLen/newPixLen;  kPspec  = fftshift(kPspec);      kPspecH = kPspec*nozHt*corctnFac;    
%  pSpectrum scaling to log & plotting the result
    pSpecNorm = fliplr(pSpec/max(max(pSpec)));   figTitle = ['$f_s=',num2str(scrFreq),'Hz$'];   cBTitle = '$|\hat{\psi}|$';
    figNo1  = contourPlot(pSpecNorm,kPspecH,yNew,kName,yName,[-18 18],schJet,figTitle,cBTitle);                  
    zLoc = find(kPspecH==0);
%  Separating & plotting upstream & downstream components
%  Downstream Component
    dwStrmRes = pSpec;   dwStrmRes(:,zLoc:end) = 0;      dwStrmRes = ifftshift(dwStrmRes);   figTitle = '$Downstream \thinspace Components$'; 
    dwStrm    = ifft(dwStrmRes,size(dwStrmRes,2),2);     dwStrm = ifftshift(dwStrm,1);       cBTitle = '$R_e(inv(\hat{\psi_d}))$';
%     figNo2    = contourPlot(real(dwStrm(:,colNos)),Xn(axSpan),yNew,xName,yName,[0 xEnd],schJet,figTitle,cBTitle);  
%     figure(figNo2);   axis equal;   xlim([0 xEnd]);   set(gcf,'Position',[70 400 715 400]);
%  Upstream Component
    upStrmRes = pSpec;   upStrmRes(:,1:zLoc)   = 0;      upStrmRes = ifftshift(upStrmRes);   figTitle = '$Upstream \thinspace Components$';      
    upStrm    = ifft(upStrmRes,size(upStrmRes,2),2);     upStrm = ifftshift(upStrm,1);       cBTitle = '$R_e(inv(\hat{\psi_u}))$';
%     figNo3    = contourPlot(real(upStrm(:,colNos)),Xn(axSpan),yNew,xName,yName,[0 xEnd],schJet,figTitle,cBTitle);  
%     figure(figNo3);   axis equal;   xlim([0 xEnd]);   set(gcf,'Position',[70 400 715 400]);
%  Computing sonic wave number(kA) - representing the speed of sound - for upstream & downstream & plotting them
    kAdw = zeros(size(Yn));   kAdw(:) = ((2*pi*scrFreq)/aAmb)*nozHt;   kAup = kAdw*-1;
    figure(figNo1);   hold on;  set(gca,'Layer','top');   plot(kAup,Yn,'-w','LineWidth',1.2);   plot(kAdw,Yn,'-w','LineWidth',1.2);
%  Outputs
    outputStruct.dwStrmRes = dwStrmRes;   outputStruct.kD   = kPspecH;   outputStruct.spSpec = pSpecNorm;      
    outputStruct.upStrmRes = upStrmRes;   outputStruct.kAct = kPspec;    outputStruct.figNo   = figNo1;
    outputStruct.xNew      = Xn(axSpan);  outputStruct.yNew = yNew;
%%  Function to plot the contours
    function figNo = contourPlot(inpMat,xAxis,yAxis,xName,yName,xLim,contourMap,contourTitle,cBarTitle)
        GF_FigurePlot(inpMat,xAxis,yAxis);    xlabel(xName,'Interpreter','latex');   axis normal;   ax = gca;    set(gcf,'Position',[35 300 900 550]);    
        ax.TickLabelInterpreter = 'latex';    ylabel(yName,'Interpreter','latex');   xlim(xLim);    ax.FontSize = 12;   curntFig = gcf;   c = colorbar;
        c.TickLabelInterpreter = 'latex';     cL = title(c,cBarTitle,'Interpreter','latex');   cL.Rotation = 90;   cL.Position = [55 170 0];  colormap(contourMap);  
        title(contourTitle,'Interpreter','latex');   cL.FontSize = 13;   figNo = curntFig.Number;
    end
end
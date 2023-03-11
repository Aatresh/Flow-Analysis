function [kVal] = GF_SFFTModePlot(funcInput)
% DESCRIPTION: FUNCTION RECONSTRUCTS THE SPATIAL WAVENUMBER MODES
% Function computes the spatial modes based on the selected wavenunmber from the upstream &
% downsteram matrices. Wavenumber selection is done directly from the plot through a mouse click
% ---------------------------------------------------------------------------
% Inputs:  Struct with Spatial Frequency; Spatial Axes; Component modes; Spatial frequency buffer
% Outputs: Plots showing the spatial distribution
%---------------------------------------------------------------------------
    X = funcInput.X;     upMat = funcInput.upStrmRes;   NF = funcInput.NF;        lenScales = funcInput.lenScales;      
    Y = funcInput.Y;     dwMat = funcInput.dwStrmRes;   normFac = lenScales(1);   nozHt = lenScales(2);   aAmb = 345.501931;
    kD = funcInput.kD;   baseFFT = funcInput.baseFFT;   load("schJet.mat");       freqVal = funcInput.fPeak;                 %#ok   
    Uj = funcInput.Uj;   figNo = funcInput.figNo;               
%  Selecting input from the spatial transform plot & locating the correct wave number
    figure(figNo);   [kInp,~] = ginput(1);   [~,posTn] = find(kD<kInp);   kVal = kD(posTn(end)+1);   clear posTn;
    disp([newline,'==> Selected Wave Number: ',num2str(kVal)]);           uProp = (2*pi*freqVal)/(abs(kVal)/nozHt);
    disp([newline,'==> Propagation velocity: ',num2str(uProp),'m/s']);    
%  Transforming the wavenumber matrix to correlate the wave number location with the base spatial spectra
    kT = ifftshift(fliplr(kD));   posTn = find(kT==kVal);   
%  Adding a buffer before & after selected wavenumber for windowing
    if strcmp(NF,'D')
       bufr = 3;
    else
       bufr = 4;
    end
    while posTn+bufr > size(kT,2)
        bufr = bufr-1;
    end
    kRng = posTn-bufr:posTn+bufr;
%  Zeroing & inverting based on the direction of propagation 
    if kVal < 0
       inpSpec = upMat;   cmp = 'u';   disp([newline,'==> Velocity fraction(U): ',num2str(uProp/aAmb),'a0']);
    else
       inpSpec = dwMat;   cmp = 'd';   disp([newline,'==> Velocity fraction(D): ',num2str(uProp/Uj),'Uj']);
    end
    holdMat = zeros(size(inpSpec));   holdMat(:,kRng) = inpSpec(:,kRng);   winFun = kaiser(size(kRng,2),10)';
    holdMat(:,kRng) = holdMat(:,kRng).*winFun;                             outMat = ifft(holdMat,size(holdMat,2),2);  
%  Plotting the result
    GF_FigurePlot(real(outMat),X,Y);   ax = gca;     ax.TickLabelInterpreter = 'latex';   colormap bluewhitered;
    xlabel(funcInput.xName,'Interpreter','latex');   ax.FontSize = 12;   cB = colorbar;   cB.TickLabelInterpreter = 'latex';
    ylabel(funcInput.yName,'Interpreter','latex');   cL = title(cB,['$R_e[I_{fft}(\hat{\psi}_',cmp,')]$'],'Interpreter','latex'); 
    cL.Rotation = 90;   cL.Position = [55 125 0];    title(['$k_',cmp,':',num2str(kVal),'$'],'Interpreter','latex');    
    set(gcf,'Position',[75 340 750 410]);
%% DEPREICATED
%  Normalizing with peak values
%     holdMat = reshape(real(outMat),size(outMat,1)*size(outMat,2),1);   [pR,~] = find(holdMat>0);   [nR,~] = find(holdMat<0);
%     pMax = max(holdMat);   holdMat(pR) = holdMat(pR)/pMax;             nMin = min(holdMat);   holdMat(nR) = holdMat(nR)/abs(nMin);
%     outMat = reshape(holdMat,size(outMat,1),size(outMat,2));
      

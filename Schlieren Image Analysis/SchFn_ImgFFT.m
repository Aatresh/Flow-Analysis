function [imgPSD,imgFFT,imgPhs,imgSpec,intgEng] = SchFn_ImgFFT(inpMat,blockSize,nOvlp,dt,baseRef)
%% Description: Function computes temporal FFT
%  Function accepts the image matrix consisting of fluctuating components
%  and computes the Fast Fourier Transform in the temporal direction based
%  on the specified block size.
%  Inputs : Image matrix, fft blocksize, block overlap sampling frequency, Base Reference
%  Outputs: fft matrix, PSD matrix, Phase matrix, Integral energy
%  distribution based on BaseReference value
%% Freq. resolution & Output matrix initialization
freqRes = (1/dt)/blockSize;         fs = 1/dt;
nOvlp = blockSize * nOvlp;                             % multiplying by factor to get no. of overlap images 
nBlks = (size(inpMat,3)- nOvlp)/(blockSize - nOvlp);   % number of blocks 
imgPSD = zeros(size(inpMat,1),size(inpMat,2),blockSize);      
imgPhs = zeros(size(inpMat,1),size(inpMat,2),blockSize);   
imgFFT = zeros(size(inpMat,1),size(inpMat,2),blockSize);   
%% FFT Computation, averaging & spectrum calculation
for iBlk = 1:nBlks
  offset = min((iBlk-1)*(blockSize - nOvlp)+blockSize,size(inpMat,3))-blockSize;
  nImgs = (1:blockSize) + offset;                      % adding offset to get correct overlap
  block = inpMat(:,:,nImgs(1):nImgs(end));
  res_fft = fft(block,blockSize,3);
  res_fft(:,:,2:end-1) = 2*res_fft(:,:,2:end-1);       % converting to single sided spectrum
  imgFFT = imgFFT + res_fft;                           % Adding fft results
  imgPSD = abs(res_fft).^2 + imgPSD;                   % Computing Power Spectral Density & block addition
  imgPhs = res_fft;                                    % Frequency specific phase 
end
imgPSD = imgPSD/(fs*blockSize*nBlks);                  % Normalizing to Pa/Hz^2 & averaging
imgPSD(:,:,2:end-1) = 2*imgPSD(:,:,2:end-1);
imgFFT = imgFFT/iBlk;                                  % Averaging by number of blocks
%% Converting to SPL
imgSpec = 10*log10(imgPSD*freqRes/(baseRef^2));
%% Integral Energy
intgEng = trapz(freqRes,imgPSD./2,3);
intgEng = 10*log10(intgEng/baseRef^2);
end

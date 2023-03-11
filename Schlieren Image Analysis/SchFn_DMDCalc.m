function [outputStruct] = SchFn_DMDCalc(Master_U,nModes,dt,matDims)
% DESCRIPTIN: FUNCTION COMPUTES THE DMD MODES FOR A GIVEN INPUT MATRIX
%  Input matrix is a reconstructed 2D-matrix with each image represented as a column & the number of
%  such columns equaling the number of snapshots.
% ---------------------------------------------------------------------------
%  Inputs : Fluctuating image matrix; No. of DMD Modes; Temporal resolution; Image matrix dimensions
%  Outputs: DMD Modes; Associated DMD Frequencies; Mode Amplitudes; Mode Energies; Growth Rate
% ---------------------------------------------------------------------------

%  Output Struct
    outputStruct = [];                  
%  Reshaping to 2D matrix with each snapshot/image as column
%     Master_U = (reshape(Master_U,dims(1),prod(dims(2:end)))).';    
%  Splitting the images into two set
%  S1 - Set without last  snapshot;   S2 - Set without first snapshot
    S1 = Master_U(:,1:end-1);             S2 = Master_U(:,2:end);    clear Master_U
%  SVD Computation on first image set 
    [leftSinglr, diagMat, rightSinglr] = svd(S1,'econ');
%  POD Mode energies
    outputStruct.PODModEne = diag(diagMat).^2/max(diag(diagMat).^2);    
%  Reducing all matrix ranks to required number of modes
    leftSinglr = leftSinglr(:,1:nModes);  rightSinglr = rightSinglr(:,1:nModes);   diagMat = diagMat(1:nModes,1:nModes);
%  Computing Companion/Convulution matrix
    A_tilde = leftSinglr'*S2*rightSinglr/diagMat;           clear leftSinglr;
%  (For debugging), we can compare if A_tilde=A, for r=max(r):
%   A = S2*pinv(Snap1);
%  Extract Eigen values & Eigen vectors from A~
    [eigenVecs, eigenVals] = eig(A_tilde);                  clear A_tilde;
%  Computing the DMD Modes & extracting mode energy(discrete time Eigen values)
    allModes = S2*rightSinglr*inv(diagMat)*eigenVecs;                                           %#ok
  % dmdModes = S2*rightSinglr*eigenVecs/(diagMat);    - Alternative
    allLambda =  diag(eigenVals);         outputStruct.allLambda = allLambda;
%  Computing Mode amplitudes & 
    allModeAmp = allModes\S1(:,1);        outputStruct.allModeAmp = allModeAmp;
%  Computuing the continuous time Eigen values
    allOmega = log(allLambda)/dt;         outputStruct.allOmega  = allOmega;
%  Frequencies associated with modes
    allFreq = imag(allOmega)/(2*pi);      [rng,~] = find(allFreq>0);
%  Filtering modes & energies associated to positive frequencies
    dmdFreq = allFreq(rng);               dmdModes = allModes(:,rng);   
    dmdOmega = allOmega(rng);             dmdModeAmp = allModeAmp(rng);
    dmdModeAmp = abs(dmdModeAmp).^2/max(abs(dmdModeAmp).^2);   dmdLambda = allLambda(rng);
%  Mode growth Rates & normalization based on max & min values
    allGrwthRate = log(abs(allLambda));   dmdGrwthRate = allGrwthRate(rng);
    maxGR = max(dmdGrwthRate);            minGR = min(dmdGrwthRate);
    [posVal,~] = find(dmdGrwthRate>0);    [negVal,~] = find(dmdGrwthRate<0);
    dmdGrwthRate(posVal) = dmdGrwthRate(posVal)/maxGR;
    dmdGrwthRate(negVal) = dmdGrwthRate(negVal)/abs(minGR);
%  Rearranging & output matrices based on updated number of snapshots
    matDims(3) = size(dmdModes,2);        dmdModes = reshape(dmdModes,matDims);
    outputStruct.dmdModes = dmdModes;     outputStruct.dmdFreq = dmdFreq;
    outputStruct.dmdModeAmp = dmdModeAmp; outputStruct.dmdLambda = dmdLambda;
    outputStruct.dmdOmega = dmdOmega;     outputStruct.dmdGrwthRate = dmdGrwthRate;
end
% DEPRICATED
%  COMPUTE THE FREQUENCIES ASSOCIATED WITH THE MODES
%     Nyquist_cutoff = 1/(2*dt);     Output_Struct.Freq = (angle(log(Output_Struct.EigenVal))/pi)*Nyquist_cutoff;
%     Output_Struct.Freq = imag((1/2*dt)*log(Output_Struct.EigenVal));

function f = GF_FourierFreq(sigPts,acqFreq,aliasFreq)
%% DESCRIPTION:FUNCTION TO RETURN VECTOR OF FREQUENCIES(SPATIAL OR TEMPORAL)
%  Computes the frequency/wavenumber range depending on the number of points
%  in the domain(temporal{blockSize}/spatial{domainSize}) & the resolution
%  (Acquisition frequency/Inverse of Spatial Resolution). Also computes 
%  the frequency/wavenumber range for aliased signals{Acquisition frequency < Nyquist limit}
%  Cite: Peter Mao (2022). fftfreq (https://www.mathworks.com/matlabcentral/fileexchange/67026-fftfreq), 
%  MATLAB Central File Exchange. Retrieved December 14, 2022.
%  Inputs  : Signal Length; Inverse Resolution(1/dt or 1/dx), Aliased Resolution(optional)
%  Outputs : Frequency(Temporal)/Wavenumbers(Spatial) 

    if nargin < 3
       aliasFreq = acqFreq;
    end
%  Nyquist cut-off Frequency definition
    minFreq = -acqFreq/2;
%  Frequency resolution definition & negative frequency starting value
    df = acqFreq/sigPts;            f0 = -minFreq;
%  Nyqyuist cut-off for Aliased acquisition frequency & negative frequency starting value
    fMinAlias = -aliasFreq/2;       f0Alias = -fMinAlias;
%  Frequency range definition - trua & aliased
    fRange = mod(linspace(0,2*f0-df,sigPts)+f0, 2*f0) - f0;
    fRangeAlias = mod(fRange+f0Alias,2*f0Alias) - f0Alias;
    f = fRangeAlias;
end
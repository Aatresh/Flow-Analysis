# Flow-Analysis
## Description:
The code directory consists of scripts written in MATLAB used to post-process temporal & spatial data obtained from various expaeriments conducted at the *Heated Jet Noise Rig (HJNR)* located in the *Gas Dynamics & Propulsion Lab (GDPL)* at the University of Cincinnati. The scripts are categorized according to the type of experiment conducted on supersonic jets. Also included are some common functions utilized by all the scripts for data catergorization, location tracking, normalization parameters, length scale assignment & other parameters to aid in data representation & interpretation. The sub directories are as follows:

### *Acoustic Analysis*: 
Scripts used for processing pressure data from far field & near field experiments. Rhis includes extracting temporal frequency distribution, Spectrograms to study temporal variation in frequencies of supersonic jets, phase - coherence analysis to understand the oscillation modes of the jet & visualization of frequency specific noise radiation patterns obtained from near field experimental data.
  
### *Global Functions*: 
Common function scripts that used to conduct data organization & categorizations, computation of Spatial Fourier decomposition & mode reconstruction, computation of critical flow properties such as isentropic Mach no. & axis normalization with a specific lenght scale. 
  
### *PIV Analysis*: 
Scripts that deal with data reconstruction, vector filetering, flow property extraction (such as Turbulence Kinetic Energy, Reynold's Stresses & Vorticity) from Particle Image Velocimetry (PIV) vector data.
  
### *POD Analysis*: 
Script files that deal with implementation of Proper Orthogonal Decomposition (POD) technique of spatial energy decomposition on PIV as well as Schlieren images.
  
### *Schlieren Image Analysis*: 
Script files that deal with processing high-speed Schlieren images  including Spatial & Temporal Fourier decomposition techniques, sequential energy integral of temporal frequency energy distribution & custom implementation of [Spectral Proper Orthogonal Decomposition (SPOD) - developed by Oliver Schmidt](https://www.mathworks.com/matlabcentral/fileexchange/65683-spectral-proper-orthogonal-decomposition-spod) as well as Dynamic Mode Decomposition Techniques.

# Flow-Analysis
## Description:
The code directory consists of scripts written in MATLAB used to post-process temporal & spatial data obtained from various expaeriments conducted at the *Heated Jet Noise Rig (HJNR)* located in the *Gas Dynamics & Propulsion Lab (GDPL)* at the University of Cincinnati. The scripts are categorized according to the type of experiment conducted on supersonic jets. Also included are some common functions utilized by all the scripts for data catergorization, location tracking, normalization parameters, length scale assignment & other parameters to aid in data representation & interpretation. The sub directories are as follows:

### *Acoustic Analysis*: 
Scripts used for processing pressure data from far field & near field experiments. Rhis includes extracting temporal frequency distribution, Spectrograms to study temporal variation in frequencies of supersonic jets, phase - coherence analysis to understand the oscillation modes of the jet.  
  
### *Global Functions*: 
Common function scripts that used to conduct data organization & categorizations, computation of Spatial Fourier decomposition & mode reconstruction, computation of critical flow properties such as isentropic Mach no. & axis normalization with a specific lenght scale. 
  
### *PIV Analysis: Script files that deal with organizing, filetering, processing & data extraction from Particle Image Velocimetry (PIV) vector data.
  
  =>  POD Analysis: Script files that deal with implementation of Proper Orthogonal Decomposition (POD) technique of spatial energy decomposition on PIV as well as Schlieren images
  
  =>  Schlieren Image Analysis: Script files that deal with processing high-speed Schlieren images which includes spatial & temporal decomposition of images to extract frequency specific energy distribution in supersonic flows.

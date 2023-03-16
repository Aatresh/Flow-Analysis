# Flow-Analysis
## Description:
The code directory consists of scripts written in MATLAB used to post-process temporal & spatial data obtained from various expaeriments conducted at the *Heated Jet Noise Rig (HJNR)* located in the *Gas Dynamics & Propulsion Lab (GDPL)* at the University of Cincinnati. The scripts are categorized according to the type of experiment conducted on supersonic jets. Also included are some common functions utilized by all the scripts for data catergorization, location tracking, normalization parameters, length scale assignment & other parameters to aid in data representation & interpretation.

  =>  Acoustic Analysis: Contains scripts for processing pressure data from far field & near field experiments and extracting temporal frequency distribution of supersonic jet noise.
  
  =>  Global Functions: Contains common function scripts that are accesed by all the codes. This includes scripts that specify data location directory & critical geometric parameters for scaling & length scale normalization.
  
  =>  PIV Analysis: Script files that deal with organizing, filetering, processing & data extraction from Particle Image Velocimetry (PIV) vector data.
  
  =>  POD Analysis: Script files that deal with implementation of Proper Orthogonal Decomposition (POD) technique of spatial energy decomposition on PIV as well as Schlieren images
  
  =>  Schlieren Image Analysis: Script files that deal with processing high-speed Schlieren images which includes spatial & temporal decomposition of images to extract frequency specific energy distribution in supersonic flows.

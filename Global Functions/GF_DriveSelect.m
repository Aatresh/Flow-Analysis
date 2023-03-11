function [OutputStruct] = GF_DriveSelect(config,nozzle,code)
%% DESCRIPTION: SELECTS INPUT DRIVE BASED ON TYPE OF TEST & CONFIGURATION
%  Syntax:
%  GF_Drive_Select(nozzleConfiguration,nozzleOrientation,callingCode)
%  nozzleConfiguration: Type of nozzle. Eg: S - Single Rectangular
%  nozzleOrientation  : Major/Minor/Sym_Lin etc.
%  callingCode        : Unique code assigned to specific analysis type
%  The function assigns an input/output path based on the type of analysis
%  the calling code is performing. Outputs are sent in a Struct file &
%  include data paths & other case specific metrics. Eg: dt in case of
%  Schlieren

% Initializing nozzle for Circular case
if strcmp(config(1),'C')
   nozzle = 'Circular';
end
%%  Output Struct & other variables
OutputStruct = [];   bckgrnd_root = [];   dt = [];
%%  Root Folders
dataRoot    = 'X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\';
schRoot     = 'G:\Data\Schlieren\';           % Schlieren Raw Data
pivDatRoot  = 'G:\Data\Vector_Files\';        % PIV DAT files 
pivMatRoot  = 'G:\Data\Vector_Matrices\';     % Unfiletered PIV Matrices
pivResRoot  = [dataRoot 'PIV_Data\'];         % Root for processed PIV Data
pivPODRoot  = [dataRoot 'POD_Data\'];         % Root for processed POD Data
ffAcousRoot = [dataRoot 'Far_Field_Data\'];   % Root for Far Field Data
%%  FINDING CORRECT OUTPUT PARAMETERS
if strcmp(code(1),'S')
%%  Schlieren Configration Tables
% Configuration Legend: Active Use
% C      - CircularMedium --------------------------> Acquisition Rate: 45,000  Hz
% S      - SingleRectagualr ------------------------> Acquisition Rate: 41,000  Hz
% S2/SS2 - SingleRectagular(Schlieren/Shadowgraph) -> Acquisition Rate: 204,800 Hz
% TR     - TwinRectagualr --------------------------> Acquisition Rate: 41,000  Hz
% TR2    - TwinRectangular -------------------------> Acquisition Rate: 204,800 Hz
% TS     - TwinSquare ------------------------------> Acquisition Rate: 41,000  Hz
% TS1    - TwinSquare ------------------------------> Acquisition Rate: 112,000 Hz
% TS2    - TwinSquare ------------------------------> Acquisition Rate: 204,800 Hz
%---------------------------------------------------------------------------------------
% Legacy: 
% TRV0 - OldNozzle Twin Rectagular -----------------> Acquisition Rate: 45,000Hz
%---------------------------------------------------------------------------------------
% Configuration Table for Schlieren
nozlConfig    = {'C';'S';'S2';'SS2';'TR';'TR2';'TS';'TS1';'TS2';'TRV0'};
frameRate     = [45E3;41E3;204.8E3;204.8E3;41E3;204.8E3;41E3;112E3;204.8E3;45E3];
frameTxt      = {'45kHz';'41kHz';'204kHz';'204kHz_s';'41kHz';'204kHz';'41kHz';'112kHz';'204kHz';'45kHz'};
folderTxt     = {'SC - Circular_Med';'SC - Single_Rectan';'SC - Single_Rectan';'SC - Single_Rectan';'SC - Twin_Rectan';...
                 'SC - Twin_Rectan';'SC - Twin_Square';'SC - Twin_Square';'SC - Twin_Square';'SC - Twin_Rectan_V0'};
frameDt       = 1./frameRate;
equivDiam     = [20.6502;20.6502;20.6502;20.6502;18.4658;18.4658;18.4658;18.4658;18.4658;18.4658;];
schTestConfig = table(nozlConfig,frameRate,frameDt,frameTxt,folderTxt,equivDiam);
%---------------------------------------------------------------------------------------------------
% Finding the right configuration details
inDxSch       = strcmp(schTestConfig.nozlConfig,config);       % Location of the input configuration
frameTag      = cell2mat(schTestConfig.frameTxt(inDxSch));     % Corresponding frame rate text
folderTag     = cell2mat(schTestConfig.folderTxt(inDxSch));    % Corresponding folder name
%---------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------
% __________________________________ Schlieren POD Computation _____________________________________
   if strcmp(code,'Sch_PODCalc')
      in_root  = [schRoot folderTag '\DAT Versions\' nozzle '\' frameTag '\'];
      out_root = [dataRoot 'POD_Data\' folderTag '\' nozzle '\' frameTag '\'];
      dt       = schTestConfig.frameDt(inDxSch);
      if strcmp(config,'TRV0')
         in_root  = [schRoot folderTag '\AVI Versions\' nozzle '\' frameTag '\'];
      end
      if isfolder(in_root) ~= 1
         disp([newline '->> Check nozzle configuration']); return;
      end
   % Configurations with Background Folders
      bckgYes = {'C','S','S2','TR','TR2','TS','TS1','TS2'};  bckIndx = strcmp(config,bckgYes);
      if find(bckIndx,1) ~= 0
         bckgrnd_root = [schRoot folderTag '\DAT Versions\' nozzle '\' frameTag '\Bckgrnd\'];
      else
         bckgrnd_root = [];
      end
%---------------------------------------------------------------------------------------------------
% _________________________________  Schlieren Fourier Analysis  ___________________________________
   elseif strcmp(code,'SchFourier')
      in_root  = [schRoot folderTag '\DAT Versions\' nozzle '\' frameTag '\'];
      out_root = [dataRoot 'Schlieren\' folderTag '\' nozzle '\' frameTag '\'];
      dt       = schTestConfig.frameDt(inDxSch);
      if strcmp(config,'TRV0')
         in_root  = [schRoot folderTag '\AVI Versions\' nozzle '\' frameTag '\'];
      end
      if isfolder(in_root) ~= 1
         disp([newline '->> Check nozzle configuration']); return;
      end
      % Configurations with Background Folders
      bckgYes = {'C','S','S2','TR','TR2','TS','TS1','TS2'};  bckIndx = strcmp(config,bckgYes);
      if find(bckIndx,1) ~= 0
         bckgrnd_root = [schRoot folderTag '\DAT Versions\' nozzle '\' frameTag '\Bckgrnd\'];
      else
         bckgrnd_root = [];
      end
%---------------------------------------------------------------------------------------------------
% _______________ Plotting Schlieren POD Modes: dt has both 1/Fs & equivalent diameter _____________
   elseif strcmp(code,'SchPODPlot')
      in_root = [dataRoot 'POD_Data\' folderTag '\' nozzle '\' frameTag '\']; out_root = in_root;
      dt      = [schTestConfig.frameDt(inDxSch),schTestConfig.equivDiam(inDxSch)];
      if isfolder(in_root) ~= 1
         disp([newline '->> Check nozzle configuration']); return;
      end
%---------------------------------------------------------------------------------------------------
% __________________________________ Schlieren Spectral POD ________________________________________
   elseif strcmp(code,'SPOD')
      in_root  = [schRoot folderTag '\DAT Versions\' nozzle '\' frameTag '\'];
      out_root = [dataRoot 'SPOD\' folderTag '\' nozzle '\' frameTag '\'];
      dt       = schTestConfig.frameDt(inDxSch);
      if strcmp(config,'TRV0')
         in_root  = [schRoot folderTag '\AVI Versions\' nozzle '\' frameTag '\'];
      end
      if isfolder(in_root) ~= 1
         disp([newline '->> Check nozzle configuration']); return;
      end
      % Configurations with Background Folders
      bckgYes = {'C','S','TR','TR2','TS','TS2'};  bckIndx = strcmp(config,bckgYes);
      if find(bckIndx,1) ~= 0
         bckgrnd_root = [schRoot folderTag '\DAT Versions\' nozzle '\' frameTag '\Bckgrnd\'];
      else
         bckgrnd_root = [];
      end
%---------------------------------------------------------------------------------------------------
% __________________________________ Schlieren DMD _________________________________________________
   elseif strcmp(code,'SchDMD')
      in_root  = [schRoot folderTag '\DAT Versions\' nozzle '\' frameTag '\'];
      out_root = [dataRoot 'DMD_Data\' folderTag '\' nozzle '\' frameTag '\'];
      dt       = schTestConfig.frameDt(inDxSch);
      if strcmp(config,'TRV0')
         in_root  = [schRoot folderTag '\AVI Versions\' nozzle '\' frameTag '\'];
      end
      if isfolder(in_root) ~= 1
         disp([newline '->> Check nozzle configuration']); return;
      end
      % Configurations with Background Folders
      bckgYes = {'C','S','TR','TR2','TS','TS2'};  bckIndx = strcmp(config,bckgYes);
      if find(bckIndx,1) ~= 0
         bckgrnd_root = [schRoot folderTag '\DAT Versions\' nozzle '\' frameTag '\Bckgrnd\'];
      else
         bckgrnd_root = [];
      end
%---------------------------------------------------------------------------------------------------
% _______________________ Simultaneous Acoustics accquirred with Schlieren _________________________
   elseif strcmp(code,'Sch_SimlAcs')
      in_root  = [dataRoot 'Schil_Simult_Acoustics\' folderTag '\'  nozzle '\' frameTag '\'];   
      out_root = [dataRoot 'Schil_Simult_Acoustics\' folderTag '\Results\'  nozzle '\' frameTag '\'];  
      dt       = 0.0184658;  % Equivalent diameter for Strouhal number calculation
   end
elseif strcmp(code(1),'P')
%%  PIV Configuration Tables
% Configuration Naming
% Syntax: {Configuration}-'L'{LensFcLen}{ImageState}{FinalWindow}{ProessingType}{Resolution*}. 
%   Configuration : Type of nozzle. Eg: C - Circular
%   LensFcLen     : Focal length of lens used. Eg: 50 - 50mm
%   ImageState    : Individual(T:TwinCam) or Merged(M)
%   FinalWindow   : 8(8x8 pxWindow) or 16(16x16 pxWindow)
%   ProcessingType: Davis(D) or Self Processed(no end letter)
%   Resolution    : *For Davis results only - full resolution statistics
%    Eg.1: C-L50M8D - Circular, 50mm FOV, merged, 8x8 window, Davis Results
%    Note: By default all Circular nozzles are taken with 50mm('L50' not included in 
%          configuration name)
% --------------------------------------------------------------------------------------------------
% Configuration Legend: Circular Medium
% C1  - CM8  : Computed{Merged}      |  C6  - CM16  : Computed{Merged}
% C2  - CM8D : Davis{Merged}         |  C7  - CM16D : Davis{Merged}
% C3  - CM8DF: Davis(FullRes){Merged}|  C8  - CM16DF: Davis(FullRes){Merged}
% C4  - CT8  : Computed{TwinCam}     |  C9  - CT16  : Computed{TwinCam}
% C5  - CT8D : Davis{TwinCam}        |  C10 - CT16D : Davis{TwinCam}
%---------------------------------------------------------------------------------------------------
% Circular Medium configurations table
   circTypes     = {'C1';'C2';'C3';'C4';'C5';'C6';'C7';'C8';'C9';'C10'};
   circCamState  = {'Merged';'Merged';'Merged';'TwinCam';'TwinCam';'Merged';'Merged';'Merged';'TwinCam';'TwinCam'};
   circProcsWin  = {'8by8';'8by8';'8by8';'8by8';'8by8';'16by16';'16by16';'16by16';'16by16';'16by16'};
   circCompTyp   = {'';'Davis';'Davis\Full_Res';'';'Davis';'';'Davis';'Davis\Full_Res';'';'Davis'};
   circPivConfig = table(circTypes,circCamState,circProcsWin,circCompTyp);
%---------------------------------------------------------------------------------------------------
% Configuration Legend: Circular Large
% CL1  - CM8  : Computed{Merged}      |  CL6  - CM16  : Computed{Merged}
% CL2  - CM8D : Davis{Merged}         |  CL7  - CM16D : Davis{Merged}
% CL3  - CM8DF: Davis(FullRes){Merged}|  CL8  - CM16DF: Davis(FullRes){Merged}
% CL4  - CT8  : Computed{TwinCam}     |  CL9  - CT16  : Computed{TwinCam}
% CL5  - CT8D : Davis{TwinCam}        |  CL10 - CT16D : Davis{TwinCam}
%---------------------------------------------------------------------------------------------------
% Circular Large configurations table
   circLrgTypes     = {'CL1';'CL2';'CL3';'CL4';'CL5';'CL6';'CL7';'CL8';'CL9';'CL10'};
   circLrgCamState  = {'Merged';'Merged';'Merged';'TwinCam';'TwinCam';'Merged';'Merged';'Merged';'TwinCam';'TwinCam'};
   circLrgProcsWin  = {'8by8';'8by8';'8by8';'8by8';'8by8';'16by16';'16by16';'16by16';'16by16';'16by16'};
   circLrgCompTyp   = {'';'Davis';'Davis\Full_Res';'';'Davis';'';'Davis';'Davis\Full_Res';'';'Davis'};
   circLrgPivConfig = table(circLrgTypes,circLrgCamState,circLrgProcsWin,circLrgCompTyp);
%---------------------------------------------------------------------------------------------------
% Configuration Legend: Twin Rectangular
% TR1  - TR-L50T16 : Computed{TwinCam} |  TR6  - TR-L50M16D  : Davis{Merged}
% TR2  - TR-L50M16 : Computed{Merged}  |  TR7  - TR-L105S8D  : Davis{SinglCam}
% TR3  - TR-L105S8 : Computed(SinglCam}|  TR8  - TR-L105S6D  : Davis{SinglCam}
% TR4  - TR-L105S6 : Computed{SinglCam}|  TR9  - TR-L105S8J2 : Computed{SinglCam}
% TR5  - TR-L50T16D: Davis{TwinCam}    |  TR10 - TR-L105S8J2D: Davis{SinglCam}
%---------------------------------------------------------------------------------------------------
% Twin Rectangular configurations table
   twRecTypes     = {'TR1';'TR2';'TR3';'TR4';'TR5';'TR6';'TR7';'TR8';'TR9';'TR10'};
   twRecCamState  = {'TwinCam';'Merged';'SinglCam';'SinglCam';'TwinCam';'Merged';...
                    'SinglCam';'SinglCam';'SinglCam';'SinglCam'};
   twRecProcsWin  = {'16by16';'16by16';'8by8';'6by6';'16by16';'16by16';'8by8';'6by6';'8by8';'8by8'};
   twReCompTyp    = {'';'';'';'';'Davis';'Davis';'Davis';'Davis';'';'Davis'};
   twRecLensTyp   = {'L-50';'L-50';'L-105';'L-105';'L-50';'L-50';'L-105';'L-105';'L-105';'L-105'};
   twRecPivConfig = table(twRecTypes,twRecCamState,twRecProcsWin,twReCompTyp,twRecLensTyp);
%---------------------------------------------------------------------------------------------------
% Configuration Legend: Twin Square
% TS1  - TS-L50T16 : Computed{TwinCam} 
% TS2  - TS-L50M16 : Computed{Merged}
% TS3  - TS-L50Q16 : Computed{QuadCam}
% TS3  - TS-L50M16D: Davis{Merged}
%---------------------------------------------------------------------------------------------------
% Twin Square configurations table
   twSquTypes     = {'TS1';'TS2';'TS3';'TS4'};
   twSquCamState  = {'TwinCam';'Merged';'QuadCam';'Merged'};
   twSquProcsWin  = {'16by16';'16by16';'16by16';'16by16'};
   twSquCompTyp   = {'';'';'';'Davis'};
   twSquLensTyp   = {'L-50';'L-50';'L-50';'L-50'};
   twSquPivConfig = table(twSquTypes,twSquCamState,twSquProcsWin,twSquCompTyp,twSquLensTyp);
%---------------------------------------------------------------------------------------------------
% Configuration Legend: Single Rectangular
% S1   - S-L50Q16  : Computed{QuadCam}
% S2   - S-L50T16  : Computed{TwinCam}
% S3   - S-L50M16  : Computed{Merged}
% S4   - S-L50M16D : Davis{Merged}
% S5   - S-L50M16DF: Davis(Full_Res){Merged}
%---------------------------------------------------------------------------------------------------
% Single Rectangalar configuration table
   siRecTypes     = {'S1';'S2';'S3';'S4';'S5'}; 
   siRecCamState  = {'QuadCam';'TwinCam';'Merged';'Merged';'Merged'};
   siRecProcsWin  = {'16by16';'16by16';'16by16';'16by16';'16by16'};
   siRecCompTyp   = {'';'';'';'Davis';'Davis\Full_Res'};
   siRecLensTyp   = {'L-50';'L-50';'L-50';'L-50';'L-50'};
   siRecPivConfig = table(siRecTypes,siRecCamState,siRecProcsWin,siRecCompTyp,siRecLensTyp);
%---------------------------------------------------------------------------------------------------
% -------------------------------- All nozzle configurations ---------------------------------------
   nozlConfig    = {'Circular_Med';'Circular_Lrg';'Single_Rectan';'Twin_Rectan';'Twin_Square'};
%---------------------------------------------------------------------------------------------------
% Finding the right configuration details
   if strcmp(config(1),'C') && length(config) == 2
      inDxPiv   = strcmp(circPivConfig.circTypes,config);        % Location of input configuration 
      camState  = cell2mat(circPivConfig.circCamState(inDxPiv)); % Merged or TwinCam camera types
      procWindw = cell2mat(circPivConfig.circProcsWin(inDxPiv)); % Final processing window
      calcType  = cell2mat(circPivConfig.circCompTyp(inDxPiv));  % Computed/Davis result
      nozlType  = nozlConfig{1};   lensTyp = 'L-50';             % Nozzle config & lens type
      Deq       = 20.6502;                                       % Nozzle exit diameter
   elseif strcmp(config(1),'C') && length(config) > 2
      inDxPiv   = strcmp(circLrgPivConfig.circLrgTypes,config);        % Location of input configuration 
      camState  = cell2mat(circLrgPivConfig.circLrgCamState(inDxPiv)); % Merged or TwinCam camera types
      procWindw = cell2mat(circLrgPivConfig.circLrgProcsWin(inDxPiv)); % Final processing window
      calcType  = cell2mat(circLrgPivConfig.circLrgCompTyp(inDxPiv));  % Computed/Davis result
      nozlType  = nozlConfig{2};   lensTyp = 'L-50';                   % Nozzle config & lens type
      Deq       = 27.559;                                              % Nozzle exit diameter
   elseif strcmp(config(1),'S')
      inDxPiv   = strcmp(siRecPivConfig.siRecTypes,config);        % Location of input configuration 
      camState  = cell2mat(siRecPivConfig.siRecCamState(inDxPiv)); % Merged or TwinCam camera types
      procWindw = cell2mat(siRecPivConfig.siRecProcsWin(inDxPiv)); % Final processing window
      calcType  = cell2mat(siRecPivConfig.siRecCompTyp(inDxPiv));  % Computed/Davis result
      lensTyp   = cell2mat(siRecPivConfig.siRecLensTyp(inDxPiv));  % Lens type
      nozlType  = nozlConfig{3};                                   % Nozzle config
      Deq      = 20.6502;                                          % Nozzle equivalent diameter
   elseif strcmp(config(1:2),'TR')
      inDxPiv   = strcmp(twRecPivConfig.twRecTypes,config);        % Location of input configuration 
      camState  = cell2mat(twRecPivConfig.twRecCamState(inDxPiv)); % Merged or TwinCam camera types
      procWindw = cell2mat(twRecPivConfig.twRecProcsWin(inDxPiv)); % Final processing window
      calcType  = cell2mat(twRecPivConfig.twReCompTyp(inDxPiv));   % Computed/Davis result
      lensTyp   = cell2mat(twRecPivConfig.twRecLensTyp(inDxPiv));  % Lens type
      nozlType  = nozlConfig{4};                                   % Nozzle config
      Deq      = 18.4658;                                          % Nozzle equivalent diameter
   elseif strcmp(config(1:2),'TS')
      inDxPiv   = strcmp(twSquPivConfig.twSquTypes,config);        % Location of input configuration 
      camState  = cell2mat(twSquPivConfig.twSquCamState(inDxPiv)); % Merged or TwinCam camera types
      procWindw = cell2mat(twSquPivConfig.twSquProcsWin(inDxPiv)); % Final processing window
      calcType  = cell2mat(twSquPivConfig.twSquCompTyp(inDxPiv));   % Computed/Davis result
      lensTyp   = cell2mat(twSquPivConfig.twSquLensTyp(inDxPiv));  % Lens type
      nozlType  = nozlConfig{5};                                   % Nozzle config
      Deq      = 18.4658;                                          % Nozzle equivalent diameter
   end
%---------------------------------------------------------------------------------------------------
% ____________________________ PIV - Saving .DAT files as .mat files _______________________________
   if strcmp(code,'PIV_VectSave')
      prompt   = [newline,'->> Enter Position(D/U): ']; pos = [input(prompt,'s'),'\'];
      in_root  = [pivDatRoot nozlType '\' lensTyp '\' pos camState '\' procWindw '\' nozzle '\'];
      out_root = [pivMatRoot nozlType '\' lensTyp '\' pos camState '\' procWindw '\' nozzle '\'];
      if strcmp(camState,'TwinCam') || strcmp(camState,'QuadCam') 
         camNos = 1:2;
      else
         camNos = 1;
      end;  dt = camNos;  bckgrnd_root = pos; % Camera no. tied to dt; Position -> bckgrnd_root
%---------------------------------------------------------------------------------------------------
% ____________________________ PIV - Compute Statistics from Vector Fields _________________________
   elseif strcmp(code,'PIV_VectCalc')
      prompt   = [newline,'->> Enter Position(D/U): ']; pos = [input(prompt,'s'),'\'];
      in_root  = [pivMatRoot nozlType '\' lensTyp '\' pos camState '\' procWindw '\' nozzle '\'];
      out_root = [pivResRoot nozlType '\' lensTyp '\' calcType '\' camState '\' procWindw '\' nozzle '\'];
      if strcmp(camState,'TwinCam') || strcmp(camState,'QuadCam') 
         camNos = 1:2;
      else
         camNos = 1;
      end;  dt = camNos;  bckgrnd_root = pos; % Camera no. tied to dt; Position -> bckgrnd_root
%---------------------------------------------------------------------------------------------------
% ________________________ PIV - Image Stitching/Axis Shifting/Plotting Results ____________________
   elseif strcmp(code,'PIV_Stitch') || strcmp(code,'PIV_AxisShift') || strcmp(code,'PIV_Utilization')
      in_root = [pivResRoot nozlType '\' lensTyp '\' calcType '\' camState '\' procWindw '\' nozzle '\'];
      out_root = in_root;       dt = Deq;     % Assigning equivalent diameter to dt
%---------------------------------------------------------------------------------------------------
% _______________________________ PIV - Computation of POD Modes ___________________________________
   elseif strcmp(code,'PIV_PODCalc')
      prompt   = [newline,'->> Enter Position(D/U): ']; pos = [input(prompt,'s'),'\'];
      in_root  = [pivMatRoot nozlType '\' lensTyp '\' pos camState '\' procWindw '\' nozzle '\'];
      out_root = [pivPODRoot 'PIV - ' nozlType '\' lensTyp '\' calcType '\' camState '\' procWindw '\' nozzle '\'];
      bckgrnd_root = [pivResRoot nozlType '\' lensTyp '\' calcType '\' camState '\' procWindw '\' nozzle '\'];
      dt = Deq;  % Assigning equivalent diameter to dt; PIV folder -> bckgrnd_root for Velocity filters
%---------------------------------------------------------------------------------------------------
% ______________________________ PIV - Plot POD Modes computed from PIV ____________________________
   elseif strcmp(code,'PIV_PODPlot')
      in_root = [pivPODRoot 'PIV - ' nozlType '\' lensTyp '\' calcType '\' camState '\' procWindw '\' nozzle '\'];
      out_root = in_root;       dt = Deq;     % Assigning equivalent diameter to dt
   end
elseif strcmp(code(1),'F')
%%  Far Field Configuration Tables
   nozlTypes    = {'C';'S';'TR';'TS'};                     % Type of nozzle 
   ffDistances  = [32.52;44.715;44.715;44];                % Distance of Far field microphone array
   nozlEqDiam   = [20.6502;20.6502;18.4658;18.5648];       % Nozzle equivalent diameter (mm)
   ffConfigSet  = table(nozlTypes,ffDistances,nozlEqDiam);
% -------------------------------- All nozzle configurations ---------------------------------------
   nozlConfig    = {'Circular_Med';'Single_Rectan';'Twin_Rectan';'Twin_Square'};
%---------------------------------------------------------------------------------------------------
 % Finding the right configuration details
   if strcmp(config,'C')
      inDx     = strcmp(ffConfigSet.nozlTypes,config);
      Deq      = ffConfigSet.nozlEqDiam(inDx)/1000;
      ffDistnc = ffConfigSet.ffDistances(inDx);
      ffAngles = [60 70 90 100 105 110 116 120 124 128 132 136 140 144 148 152];
      nozlType = nozlConfig{inDx};          resultVer = '';
   elseif strcmp(config,'S')
      inDx     = strcmp(ffConfigSet.nozlTypes,config);
      Deq      = ffConfigSet.nozlEqDiam(inDx)/1000;
      ffDistnc = ffConfigSet.ffDistances(inDx); 
      ffAngles = [45 60 70 90 100 105 110 116 120 126 132 136 140 144 148 152];
      nozlType = nozlConfig{inDx};          
      versions = {'V2_New_Chamber';'V3_New_Chamber';'V4_New_Chamber'};
      prompt   = [newline,'->> Enter Result Version(V2(1),V3(2) or V4(3) - '];  chc = input(prompt);
      resultVer = versions{chc};
   elseif strcmp(config(1),'T')
      inDx     = strcmp(ffConfigSet.nozlTypes,config);
      Deq      = ffConfigSet.nozlEqDiam(inDx)/1000;
      ffDistnc = ffConfigSet.ffDistances(inDx); 
      ffAngles = [45 60 70 90 100 105 110 116 120 126 132 136 140 144 148 152];
      nozlType = nozlConfig{inDx};          resultVer = '';
   end
%---------------------------------------------------------------------------------------------------
% _______________________ Far Field Result Computation from .dat files _____________________________
   if strcmp(code,'FF_Comp')
      in_root  = [ffAcousRoot nozlType '\' nozzle '\'];
      out_root = [ffAcousRoot nozlType '\Results\' resultVer '\' nozzle '\'];
      dt = Deq;                                      % Assigning Equivalent Diameter to dt    
      bckgrnd_root = [];
%---------------------------------------------------------------------------------------------------
% ________________________________ Far Field Result Plotting _______________________________________
   elseif strcmp(code,'FF_Plot')
      in_root  = [ffAcousRoot nozlType '\Results\' resultVer '\' nozzle '\'];
      out_root = in_root;
      dt.Deq   = Deq;       dt.ffDistnc = ffDistnc;  % Assigning Equivalent Diameter & far field
                                                     % microphone distance to dt
      bckgrnd_root = ffAngles;                       % Assigning Far Field angles to bckgrnd_root
   end
end
%% Defining Output Struct file with relavent details
OutputStruct.in_root  = in_root;      OutputStruct.bckgrnd_root = bckgrnd_root;
OutputStruct.out_root = out_root;     OutputStruct.dt           = dt;
OutputStruct.nozzle   = nozzle;      
end
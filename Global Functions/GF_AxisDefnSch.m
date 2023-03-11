function [Xn,Yn,limX,limY,xName,yName,lenScales,figSize] = GF_AxisDefnSch(config,nozzle,NF,X,Y)
% DESCRIPTION: Define axis normalization & naming
% Function defines the axial limits & normalization factor 
% based on nozzle configuration for Schlieren data
%--------------------------------------------------------------------------------------------------
% Assigning nozzle to Circular
if strcmp(config,'C')
   nozzle = 'Circular';
end
% Length scale parameter that contains normalizing length scale for St &
% height of the nozzle
lenScales = zeros(1,2);
% ------------------------------- Nozzle configuration table ---------------------------------------
nozlConfig = {'C';    'S';    'S2';  'SS2';  'TR';   'TR2';  'TS';   'TS1';  'TS2';  'TRV0'};
equivDiam  = [20.6502;20.6502;20.6502;20.6502;18.4658;18.4658;18.4658;18.4658;18.4658;18.4658];
nozlHeight = [20.6502;12.954 ;12.954 ;12.954 ;12.19  ;12.19  ;16.61  ;16.61  ;16.61  ;12.19];
xLimitUpH  = [6.5    ;12     ;4.5    ;4.5    ;12     ;4.95   ;9.2    ;6.5    ;3.65   ;12   ];
xLimitUpD  = [6.5    ;6.5    ;2.82   ;2.82   ;8      ;3.28   ;8      ;6      ;3.3    ;8    ];
yLimitH    = [2.2    ;3.5    ;1.4    ;1.4    ;4      ;1.4    ;3      ;1      ;1      ;4    ];
yLimitD    = [2.2    ;2.2    ;0.85   ;0.85   ;2.8    ;0.932  ;3.5    ;0.95   ;0.9    ;2.8  ];
xLimitDwn  =  -0.05;
axisNameD  = {'$Y/D_e$' ;'$Z/D_e$';'$Y/D_e$'};
axisNameH  = {'$Y/D_e$' ;'$Z/h_j$';'$Y/h_j$'};
nozlType   = {'Circular';'Major'  ;'Minor'};
%---------------------------------------------------------------------------------------------------
% Finding the right configuration details
inDx1         = strcmp(config,nozlConfig);               % Location of the input configuration
lenScales(1)  = equivDiam(inDx1)/1000;                   % Normalizing length(m) used for Strouhal no.
lenScales(2)  = nozlHeight(inDx1)/1000;                  % Height of the smaller nozzle dimension
figSize       = [160 380 620 405];                       % Size of the output figure
inDx2         = strcmp(nozzle,nozlType);                 % Index for checking yName
% Normalizing factor condition(equivalent diameter/nozzle height)
if strcmp(NF,'D')
   limX  = [xLimitDwn xLimitUpD(inDx1)];          limY  = [-yLimitD(inDx1) yLimitD(inDx1)];
   xName = '$X/D_e$';                             yName = cell2mat(axisNameD(inDx2));
   Xn    = X/equivDiam(inDx1);                    Yn    = Y/equivDiam(inDx1);
else
   limX  = [xLimitDwn xLimitUpH(inDx1)];          limY  = [-yLimitH(inDx1) yLimitH(inDx1)];
   xName = '$X/h_j$';                               yName = cell2mat(axisNameH(inDx2));
   Xn    = X/nozlHeight(inDx1);                   Yn    = Y/nozlHeight(inDx1);
   if strcmp(nozzle,'Circular')
      disp([newline '->> Used De as Normalizing length for Circular' newline]);
      xName = '$X/D_e$';
   end
end
end


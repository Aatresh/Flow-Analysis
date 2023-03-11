function [U,V,img_inf,bad_U,bad_V,U_imgset,Vpos_imgset,Vneg_imgset ] = GF_PIV_ImageSort(InputStruct)
% DESCRIPTION:FUNCTION TO SORT IMAGES FOR PIV CALCULATION VERSION - 2
% Comapres the images to a set of limit parameters for velocities &
% filters images that do not meet that criteria

U = InputStruct.U;   V = InputStruct.V;   camera_number = InputStruct.n_cam;
Upos_C1 = InputStruct.Upos_C1;   Upos_C2 = InputStruct.Upos_C2;  
Vpos_C1 = InputStruct.Vpos_C1;   Vneg_C1 = InputStruct.Vneg_C1;
Vpos_C2 = InputStruct.Vpos_C2;   Vneg_C2 = InputStruct.Vneg_C2;
%  COUNTER VARIABLE & INFINITE IMAGE SET
a = 1;   img_inf = [];  
%  SETTING CUTOFF BASED ON FINAL WINDOW SIZE
configPool = {'C1';'C2';'C4';'TR3';'TR4';'TR9'}; inDx = strcmp(InputStruct.config,configPool);
if isempty(find(inDx ~= 0, 1))
   Top_Cutoff = 1;      Bottom_Cutoff = -0.8;
else
   Top_Cutoff = 2.5;    Bottom_Cutoff = -2.5;
end
%% MEAN VELOCITY BEFORE SORTING
  U_bar = mean(U,3);           V_bar = mean(V,3);
  U_prim = zeros(size(U));     V_prim = zeros(size(U));
  U_temp = zeros(1,size(U,2)); V_temp = zeros(1,size(V,2));     
%%  IDENTIFYING IMAGES WITH INFINITE VECTORS - LAYER 1
  for countr2 = 1:size(U,3)                                     
    U_prim(:,:,countr2) = U(:,:,countr2) - U_bar;
    V_prim(:,:,countr2) = V(:,:,countr2) - V_bar;
    for j = 1:size(U,1)
      U_temp = U_temp + U_prim(j,:,countr2);
      V_temp = V_temp + V_prim(j,:,countr2);
    end
%  REMOVES IMAGES WITH INFINITE VELOCITY VALUES
    if max(U_temp)/10^4 > Top_Cutoff || max(U_temp)/10^4 < 0 || min(U_temp)/10^4 < Bottom_Cutoff
       img_inf(a) = countr2;   a = a+1;               %#ok
    elseif max(V_temp)/10^4 > Top_Cutoff || max(V_temp)/10^4 < 0 || min(V_temp)/10^4 < Bottom_Cutoff
       img_inf(a) = countr2;   a = a+1;               %#ok
    end
    U_temp = zeros(1,size(U,2));       V_temp = zeros(1,size(V,2));
  end
%  REMOVING IDENTIFIED IMAGES - LAYER 1
  U(:,:,img_inf) = [];    V(:,:,img_inf) = [];  disp(['->> Images with Infinite Vectors: ',num2str(size(img_inf,2))]);
  U_bar = mean(U,3);         
%% IDENTIFYING IMAGES WITH BAD VECTORS U VEL - LAYER 2
  U_imgset = zeros(1,size(U,3));         bad_U = 0;
  for countr2 = 1:size(U,3)                    
    U_imgset(1,countr2) = max(max(U(:,:,countr2)));
  end
  Upeak = max(max(U_bar));    U_imgset = U_imgset./Upeak;
  j = 1;
%  U VELOCITY LIMIT ON FILTERING IMAGES
  for countr2 = 1:length(U_imgset)
    if camera_number == 1
       if U_imgset(countr2) > Upos_C1
          bad_U(j) = countr2;                         %#ok
          j = j+1;
       end
    else
       if U_imgset(countr2) > Upos_C2
          bad_U(j) = countr2;                         %#ok
          j = j+1;
       end
    end
  end
%  CHECK IF THERE ARE NON-ZERO IMAGES
  if bad_U == 0
     disp('!!--- Increase/Decrease U cut-ff Limit ---!!');
     return;
  else
     disp([newline '->> Images cut after U - filtering: ',num2str(length(bad_U))]);
  end
%  REMOVING IDENTIFIED IMAGES - LAYER 2
  U(:,:,bad_U) = [];      V(:,:,bad_U) = [];
%% IDENTIFYING IMAGES WITH BAD VECTORS V VEL - LAYER 2
  Vpos_imgset = zeros(1,size(U,3));    Vneg_imgset = zeros(1,size(V,3));
  V_bar = mean(V,3);            bad_V = 0;
  for countr2 = 1:size(V,3)
    Vpos_imgset(1,countr2) = max(max(V(:,:,countr2)));          %  POSITIVE V VALUES     
    Vneg_imgset(1,countr2) = min(min(V(:,:,countr2)));          %  NEGATIVE V VALUES
  end
  Vpeak_pos = max(max(V_bar));      Vpeak_neg = min(min(V_bar));
  Vpos_imgset = Vpos_imgset./Vpeak_pos;           Vneg_imgset = Vneg_imgset./Vpeak_neg;
  j = 1;                                             
%  V VELOCITY LIMIT ON FILTERING IMAGES
  for countr2 = 1:length(Vpos_imgset)
    if camera_number == 1
       if Vpos_imgset(countr2) > Vpos_C1 || Vneg_imgset(countr2) > Vneg_C1
          bad_V(j) = countr2;            j = j+1;     %#ok
       end
    else
       if Vpos_imgset(countr2) > Vpos_C2 || Vneg_imgset(countr2) > Vneg_C2
          bad_V(j) = countr2;            j = j+1;     %#ok
       end     
    end
  end
  if bad_V == 0
     disp('!!--- Increase/Decrease V cut-ff Limit ---!!');
     return;
  else
     disp([newline '->> Images cut after V - filtering: ',num2str(length(bad_V))]);
  end
%  REMOVING IDENTIFIED IMAGES - LAYER 3
  U(:,:,bad_V) = [];        V(:,:,bad_V) = []; 
end
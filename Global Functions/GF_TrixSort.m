function [ X_full,Y_full,V_full ] = GF_TrixSort(Input)
% DESCRIPTION: FUNCTION TO CHANGE ORDER OF .DAT DAVIS FILES
% Function accepts the 1D DaVis .DAT file & sorts the values into a 2D
% matrix for easy storage & understanding

%  COUNTERS & HOLDING MATRIX
ctr1 = 1;   ctr2 = 0;   hold = zeros(length(Input),3);
%  COUNTER FOR X - CO-ORDINATE
for ctr3 = 1:size(Input,1)
  if Input(ctr3,2,:) ~= Input(ctr3+1,2,:)
     xlen = ctr3;  break;
  end
end
%  COUNTER FOR Y - CO-ORDINATE
ylen = round(size(Input,1)/xlen);
%  REARRANING MATRIX ACCORDING TO X-LOCATION
for ctr3 = 1:xlen
  for ctr4 = ctr3:xlen:size(Input,1)
    hold(ctr1,:,:) = Input(ctr4,:,:);
    ctr1 = ctr1+1;
  end
end
%%  SEPARATING VELOCITY COMPONENTS 
x = hold(:,1);   y = hold(:,2);   z = hold(:,3);
%  INITIALIZING LOCATION & VELOCITY MATRICES
X_full = zeros(1,xlen);   Y_full = zeros(ylen,1);   V_full = zeros(ylen,xlen);
%  SPLITTING INTO LOCATION AND VELOCITY MATRICES
for ctr5 = 1:xlen
  for ctr6 = 1:ylen
    V_full(ctr6,ctr5) = z(ctr6+ctr2);       X_full(1,ctr5) = x(ctr5*ctr6);
  end
 ctr2 = ctr2 + ylen;
end
for ctr4 = 1:ylen
  Y_full(ctr4,1) = y(ctr4);
end
end


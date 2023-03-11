function [jet_edge,Y] = GF_CenterFinder(U,Y,config,nozzle,typ)
%  DESCRIPTION: Function accepts the average velocity or tubulence profile & 
%  finds the jet centerline based on the radial distrbution at the jet exit
x = 1;            strtPostn = 5;
%  VELOCITY OR TKE LIMIT VALUE BASED ON INPUT CONTOUR
if strcmp(typ,'Ux')
   CheckLim = 250;
else
   CheckLim = 0.013;
end
%  FINDING CENTER FOR MINOR AXIS OR SINGLE JET CONFIGURATION
if strcmp(nozzle,'Minor') || strcmp(config(1),'S') || strcmp(config(1),'C') || strcmp(config,'TR3')  ||...
   strcmp(config,'TR4')   || strcmp(config,'TR5')  || strcmp(config,'TR9')  || strcmp(config,'TR10') ||...
   strcmp(nozzle,'Minor2')
%  DECIDING VELOCITY LIMIT FOR SHEAR LAYER IDENTIFICATION BASED ON MEASUREMENT PLANE
   if strcmp(nozzle,'Sym_Lin')
      CheckLim = 90;
   end
%  NaN CHECK
   while isnan(U(x,strtPostn))
      x = x+1;
   end
%  BOTTOM OF THE JET SHEAR LAYER
   while U(x,strtPostn) < CheckLim
      x = x+1; 
   end;             jetBot = x;
%  TOP OF THE JET SHEAR LAYER    
   while U(x,strtPostn) > CheckLim
         x = x+1;
   end;             jetTop = x-1;     clear x;
   jetCentr = round((jetBot+jetTop)/2);
%  DEFINING CENTER & EDGES
   jet_edge = zeros(3,1);
   jet_edge(1) = jetCentr;   disp('jetEdge(1) - Jet Center');                        %  JET CENTER
   jet_edge(2) = jetBot;     disp('jetEdge(2) - Bottom Shear Layer Center');         %  JET BOTTOM
   jet_edge(3) = jetTop;     disp('jetEdge(3) - Top Shear Layer');                   %  JET TOP
%  SHIFTING Y - AXIS TO MATCH CENTER OF THE JET
   Y = Y - Y(jetCentr); 
%  FINDING CENTER FOR MAJOR AXIS FOR TWIN CONFIGURATION
elseif strcmp(nozzle,'Major') && strcmp(config(1),'T')
%  BOTTOM SHEAR LAYER OF BOTTOM JET
   while U(x,strtPostn) < CheckLim
      x = x+1;
   end;             jetBot_1 = x;       x = size(U,1);
%  TOP SHEAR LAYER OF TOP JET
   while U(x,strtPostn) < CheckLim
      x = x-1;
   end;             jetTop_2 = x;       
%  IMAGE CENTER
   imgCenter = round((jetBot_1 + jetTop_2)/2);
%  TOP SHEAR LAYER OF BOTTOM JET
  x = imgCenter;
  while U(x,strtPostn) < CheckLim
     x = x-1;
  end;              jetTop_1 = x;       x = imgCenter;
%  BOTTOM SHEAR LAYER OF TOP JET
  while U(x,strtPostn) < CheckLim
     x = x+1;
  end;           jetBot_2 = x;
%  CENTERS OF EACH JET
  jetCentr_1 = round((jetBot_1+jetTop_1)/2);
  jetCentr_2 = round((jetBot_2+jetTop_2)/2); 
%  COMBINING INTO AN ARRAY
  jet_edge = zeros(7,1);
  jet_edge(1) = imgCenter;      disp('jetEdge(1) - Image Center');                    %  IMAGE CENTER
  jet_edge(2) = jetCentr_1;     disp('jetEdge(2) - Jet - 1 Center(Bottom)');          %  CENTER OF BOTTOM JET
  jet_edge(3) = jetCentr_2;     disp('jetEdge(3) - Jet - 2 Center(Top)');             %  CENTER OF TOP JET
  jet_edge(4) = jetBot_1;       disp('jetEdge(4) - Jet - 1 Bottom Shear Layer');      %  BOTTOM SHEAR LAYER OF BOTTOM JET
  jet_edge(5) = jetTop_1;       disp('jetEdge(5) - Jet - 1 Top Shear Layer');         %  TOP SHEAR LAYER OF BOTTOM JET
  jet_edge(6) = jetBot_2;       disp('jetEdge(6) - Jet - 2 Bottom Shear Layer');      %  BOTTOM SHEAR LAYER OF TOP JET
  jet_edge(7) = jetTop_2;       disp('jetEdge(7) - Jet - 2 Top Shear Layer');         %  TOP SHEAR LAYER OF TOP JET
%  SHIFTING Y - AXIS TO MATCH CENTER OF THE JET
   Y = Y - Y(imgCenter); 
end
end
 

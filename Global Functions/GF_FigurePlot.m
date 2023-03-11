function GF_FigurePlot(varargin)
%  DESCRIPTION - Plots pcolor from a matrix with or without a defined axes
%  Input format: GF_FigurePlot(plotVariable,axis1,axis2)

if nargin < 2       % Only plotVariable
   figure;   pcolor(varargin{1});      axis equal;    shading interp;   colorbar;
   lim_x = [0  size(varargin{1},2)];   xlim(lim_x);   colormap jet;   
   lim_y = [0  size(varargin{1},1)];   ylim(lim_y);   box off;      
   caxis([min(min(varargin{1})) max(max(varargin{1}))]); 
elseif nargin > 2   % plotVariable with axes
   figure;   pcolor(varargin{2},varargin{3},varargin{1});   axis equal;     shading interp; 
   caxis([min(min(varargin{1})) max(max(varargin{1}))]);    colormap jet;   colorbar;
   lim_x = [min(varargin{2}) max(varargin{2})];   xlim(lim_x);   box off;  
   lim_y = [min(varargin{3}) max(varargin{3})];   ylim(lim_y);  
else
   disp([newline '=> Check axis size/format']);
end


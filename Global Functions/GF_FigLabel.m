function GF_FigLabel(figNo,xName,yName,fSize)
%  DESCRIPTION: Labels figures in LaTEX style with specified axis labels & font size
%  Font labels input with proper formatting except $$'s. They are added in the function
%  Input Format: func(figure no.,xName,yName,<font size>); <optional>
%  --------------------------------------------------------------------------------------%

% FontSize check
    if nargin < 4
       textSize = 12;
    else
       textSize = fSize;
    end
%  Check for Formatting labels & append if necessary
   if strcmp(xName(1),'$') ~= 1
      xName = ['$',xName,'$'];
   end
   if strcmp(yName(1),'$') ~= 1
       yName = ['$',yName,'$'];
   end
%  Setting the parameters
    figure(figNo);   xlabel(xName,'Interpreter','latex');   ylabel(yName,'Interpreter','latex');
    ax = gca;   ax.FontSize = textSize;   ax.TickLabelInterpreter = 'latex';
end
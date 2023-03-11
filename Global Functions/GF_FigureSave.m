function GF_FigureSave(figName,driveSave,figNo)
%  DESCRIPTION: FUNCTION TO SAVE PCOLOR/CONTOUR DATA TO DISK
%   Takes figure no. & saves them as Matlab Figures, TIFF & EPS files
prompt = [newline '->> Save Figure?(y/n) - '];    chc = input(prompt,'s'); 
if strcmp(chc,'y')
   cd(driveSave);   gcf = figure(figNo);
% .eps Format
   figEPS = [figName '.eps'];   saveas(gcf,figEPS,'epsc');
% .mat Figure
   figMAT = [figName '.fig'];
   if strcmp(figName(1:11),'SPL-Contour') == 0
      saveas(gcf,figMAT);
   end
% .tif(uncompressed) Format
   figTIF = [figName '.tif'];   print(figure(figNo), '-dtiffn', '-r300', figTIF);   
%% CHANGING DIRECTORY BACK TO GLOBAL FUNCTIONS
   cd('X:\OneDrive - University of Cincinnati\Working_Directory\Thesis\Jet_Analysis\Global_Functions\');
end
end


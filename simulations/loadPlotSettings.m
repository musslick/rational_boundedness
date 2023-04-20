
% plot settings

fontSize_title = 14;
fontSize_gca = 14;
fontSize_xlabel = 14;
fontSize_ylabel = 14;
fontSize_legend = 14;

fontName = 'Helvetica';

lineWidth = 3;

theColorMap = copper;

colors = [253 120 21; ... % orange
              31 104 172; ... % blue
              44 155 37; ... % green
              0     0   0  ; ... % black
            142 142 142; ... % grey 
            255 255 255; ... % white 
            55 146 171;] / 255; % cyan
        
cContrast1 = 1;
cContrast2 = 2;
cContrast3 = 3;
cSingle = 4;
cWeak = 5;
cWhite = 6;
cCyan = 7;

% for export
plotSettings.fontSize_title = fontSize_title;
plotSettings.fontSize_gca = fontSize_gca;
plotSettings.fontSize_xlabel = fontSize_xlabel;
plotSettings.fontSize_ylabel = fontSize_ylabel;
plotSettings.fontSize_legend = fontSize_legend;
plotSettings.fontName = fontName;
plotSettings.lineWidth = lineWidth;
plotSettings.colors = colors;
plotSettings.cContrast1 = cContrast1;
plotSettings.cContrast2 = cContrast2;
plotSettings.cContrast3 = cContrast3;
plotSettings.cSingle = cSingle;
plotSettings.cWeak = cWeak;
plotSettings.cWhite = cWhite;
%% latexify.m
function latexify(varargin)
%LATEXIFY Sets text interpreter to latex for all figures or current figure.
%   Optional name-value pairs also allow custom font and figure sizing.
%
% Author:
%   Tiger Hou

%%

defaultPlotSize = [20,15];
defaultFontSize = 16;
p = inputParser;
addOptional( p,'Initialize',1);
addParameter(p,'plotSize',defaultPlotSize);
addParameter(p,'fontSize',defaultFontSize);
parse(p,varargin{:});

if p.Results.Initialize == 0
    %% modify global figure settings
    set(groot,'defaulttextinterpreter','latex');  
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

else
    %% modify current figure and axes    
    set(findall(gcf,'-property','FontSize'),'FontSize',p.Results.fontSize)
    set(gca, 'Color', [1 1 1] )
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', p.Results.plotSize)
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0.2 1.2 p.Results.plotSize])
    set(gcf, 'PaperPosition', [0.2 1.2 p.Results.plotSize])
    
end

end %latexify.m


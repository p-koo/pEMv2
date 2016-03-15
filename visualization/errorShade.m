function H = errorShade(x, lE,uE, col, patchSaturation)
% Adapted from Rob Campbell code, at:
% http://www.mathworks.com/matlabcentral/fileexchange/26311-shadederrorbar/content/shadedErrorBar.m
hold on
faceAlpha=patchSaturation;
patchColor=col;
set(gcf,'renderer','openGL')

% 
% %Calculate the y values at which we will place the error bars
% if size(errBar,2) > 1
%     uE = y + errBar(:,1);
%     lE = y - errBar(:,2);
% else
%     uE = y + errBar;
%     lE = y - errBar;
% end


%Make the cordinats for the patch
yP = [ lE ; flipud(uE) ];
xP = [ x ; flipud(x) ];

invalid = isnan(xP) | isnan(yP) | isinf(xP) | isinf(yP);
yP(invalid) = [];
xP(invalid) = [];


H.patch = patch(xP, yP, 1, ...
    'Facecolor', patchColor,...
    'Edgecolor', 'none',...
    'Facealpha', faceAlpha);

%Make nice edges around the patch.
% H.edge(1) = plot(x, lE, '-', 'Color', edgeColor);
% H.edge(2) = plot(x, uE, '-', 'Color', edgeColor);

%The main line is now covered by the patch object and was plotted first to
%extract the RGB value of the main plot line. I am not aware of an easy way
%to change the order of plot elements on the graph so we'll just remove it
%and put it back (yuk!)
% H.mainLine = plot(x, y, 'Color', col);






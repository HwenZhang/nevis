function [X,Y] = ginput_draw(c,n)
% replicates ginput(n), but drawing line segments between points with colour c

if nargin<2, n = 1000; end

if ishold, holding = 1; else holding = 0; end
for i = 1:n
    [x,y] = ginput(1);
    if isempty(x),
        return;
    else
        hold on;
        h = plot(x,y,'.'); set(h,'color',c);
        if i==1, X = x; Y = y;
        else X = [X; x]; Y = [Y; y];
        h = plot([X(end-1) X(end)],[Y(end-1) Y(end)],'-'); set(h,'color',c);
        end
        if holding, else hold off; end
    end
    
end
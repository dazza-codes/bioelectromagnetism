function [lines] = crosshair_set(x,y,percent)
%
% crosshair_set - place a crosshair on current plot at x,y position
%
% [lines] = crosshair_set(x,y,percent)
%
% x,y are the location of the cross (default at 0,0).  The x value
% specifies the intercept of a line with the x axis, the line is
% parallel with the y axis.  The y value specifies the intercept of a
% line with the yaxis, the line is parallel with the xaxis.
%
% lines.x,lines.y are handles to the crosshair lines
%
% The length of the lines is a 'percent' of the axis range; this value
% ranges from 0 to 1, the default is 1
%

if ~exist('x','var'), x = 0; end
if ~exist('y','var'), y = 0; end
if ~exist('percent','var'), percent = 1; end

if isempty(x), x = 0; end
if isempty(y), y = 0; end
if isempty(percent), percent = 1; end

x_rng = get(gca,'Xlim') .* percent;
y_rng = get(gca,'Ylim') .* percent;

% should check whether x,y are within x_rng, y_rng

lines.x = line([x,x],[y_rng(1),y_rng(2)]);
set(lines.x,'Color','k');

lines.y = line([x_rng(1) x_rng(2)],[y y]);
set(lines.y,'Color','k');

return

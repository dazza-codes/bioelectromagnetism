function [xdata,ydata] = get_lines

% Get Line Data from Current Figure

% Lines are referenced as axis children, among other
% axis children; so first get all axis children
sibs = get(gca,'Children');

% Now search axis children for any line types.
% Because the columns of the y data matrix in a plot
% command seem to be reversed in the axis children, 
% count down from max sibs to the first sib.
lines = 0;
xdata = [];
ydata = [];

i = max(size(sibs));
while i >= 1,
  if strcmp(get(sibs(i),'Type'),'line')
    % Found a line child, but check its size
    getline = 1;
    if ~isempty(xdata),
      if isequal(size(get(sibs(i),'XData').',1),size(xdata,1)),
        getline = 1;
      else
        getline = 0;
      end
    end
    if getline,
      % OK, found a line among the axis children.
      lines = lines + 1;
      datalines(lines) = sibs(i);
      
      % put line data into a column of data.xdata|data.ydata
      xdata(:,lines) = get(sibs(i),'XData').';
      ydata(:,lines) = get(sibs(i),'YData').';
    end
  end
  i = i - 1;
end

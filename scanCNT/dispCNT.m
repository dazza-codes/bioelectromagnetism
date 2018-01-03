% dispCNT(data,srate,labels)
%               
%   data                ->  matrix of waveform data (electrode in columns)
%   srate               ->  sample rate of data
%   labels(optional)    ->  vector of numeric labels or cell array of strings 
%                           for electrodes ie. [1 2 3] OR {'Pz','Cz','Oz'}
%                           number and order of labels must match electrodes 
%                           in data matrix
%
%   -- Note: Works only with Scan 4.1+ data files
function dispCNT(data,srate,labels)

% ---- Display Parameters -----------------------

dispRate = 10; % seconds per screen
chanPage =  4; % channels displayed per page

% ---- GUI layout -------------------------------
scrsz = get(0,'ScreenSize');
figPos = [scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2];
axesPos = round([0.1*figPos(3) 0.1*figPos(4) 0.8*figPos(3) 0.8*figPos(4)]);
sliderWidth = 20;
horizSliderPos = [axesPos(1) axesPos(2)-sliderWidth axesPos(3) sliderWidth];
vertSliderPos  = [axesPos(1)+axesPos(3) axesPos(2) sliderWidth axesPos(4)];

% ---- Control Parameters -----------------------

nChan  = size(data,2);
points = size(data,1);

if (size(data,1) >= dispRate*srate)
    data = data(1:dispRate*srate,:);
end

if (nargin == 2)
    labels = [];
elseif (nargin == 3)
    if (length(labels) ~= nChan)
        labels = [];
    end
end

horizSliderStep = [srate/points  (dispRate*srate)/(points-1)];

% ---- Draw Components --------------------------

hFig = figure('MenuBar','none','Position',figPos);

hAxes = axes('Units','pixels','Position',axesPos,'XAxisLocation','top','YDir','reverse');

hVertSlider  = uicontrol(hFig,'style','slider','Position',vertSliderPos);
hHorizSlider = uicontrol(hFig,'style','slider','Position',horizSliderPos,...
    'sliderstep',horizSliderStep, ...
    'min',1,'max',points,'value',1,...
    'Callback','dispCNTcallback(''hslider'')');

dataMean = mean(data);
data = data - dataMean(ones(size(data,1),1),:);

dataMax = max(data);
data = data ./ dataMax(ones(size(data,1),1),:);

increment = [1:2:size(data,2)*2];
data = data + increment(ones(size(data,1),1),:);

time = [1/srate:1/srate:size(data,1)/srate]';
line(time,data);
set(hAxes,'XLim',[0 dispRate],'YTick',increment,'YTickLabel',labels);

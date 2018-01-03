function clk_beep(n,prompt)
% CLK_BEEP	ring console bell
%   	CLK_BEEP(n [,prompt]) 
%   	rings bell N times
%   	displays PROMPT and waits for user response

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

if nargin < 1
    help clk_beep; return
end
for i=1:n, fprintf('\a'); pause(0.25); end
if nargin > 1
    disp(prompt);
    pause
end

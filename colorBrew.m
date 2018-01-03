function cs = colorBrew(cname,ncolor,ctype)

% colorScale = colorBrew([cname] [,ncolor] [,ctype])
%
% This product includes color specifications and designs developed by
% Cynthia Brewer (http://colorbrewer.org/).  This function returns a
% matlab struct array that contains RGB colors.
%
% All input arguments are optional.  The default is to return all color
% scales.  To restrict the return struct array, specify any of the input
% options, eg:
%
% >> cs = colorBrew('Blues', 9);    % Blues scale with 9 colors
% >> colormap(cs.colorRGB)
% >> colorbar
% >> cs = colorBrew([],[],'seq');   % all sequential color scales
% cs = 
% 1x126 struct array with fields:
%     colorName
%     numColors
%     colorType
%     critVal
%     colorNum
%     colorLet
%     colorRGB
%
% color scale names (cname): 'Accent' 'Blues' 'BrBG' 'BuGn' 'BuPu'
% 'Dark2' 'GnBu' 'Greens' 'Greys' 'OrRd' 'Oranges' 'PRGn' 'Paired' 'Pastel1'
% 'Pastel2' 'PiYG' 'PuBu' 'PuBuGn' 'PuOr' 'PuRd' 'Purples' 'RdBu' 'RdGy'
% 'RdPu' 'RdYlBu' 'RdYlGn' 'Reds' 'Set1' 'Set2' 'Set3' 'Spectral' 'YlGn'
% 'YlGnBu' 'YlOrBr' 'YlOrRd'
%
% number of colors (ncolor): For each color scale name there is a set of
% color scales with anywhere from 3-12 color gradations.  Not all color
% scales support all of these levels, most have 3-9 levels.
%
% color scale types (ctype):
% 'div'  - diverging
% 'qual' - qualitative
% 'seq'  - sequential
%
%
% Copyright (c) 2002 Cynthia Brewer, Mark Harrower, and The Pennsylvania
% State University. All rights reserved.  Redistribution and use in source
% and binary forms, with or without modification, are permitted provided
% that the following conditions are met: 1. Redistributions as source code
% must retain the above copyright notice, this list of conditions and the
% following disclaimer.  2. The end-user documentation included with the
% redistribution, if any, must include the following acknowledgment: This
% product includes color specifications and designs developed by Cynthia
% Brewer (http://colorbrewer.org/).  Alternately, this acknowledgment may
% appear in the software itself, if and wherever such third-party
% acknowledgments normally appear.  4. The name "ColorBrewer" must not be
% used to endorse or promote products derived from this software without
% prior written permission. For written permission, please contact Cynthia
% Brewer at cbrewer@psu.edu.  5. Products derived from this software may not
% be called "ColorBrewer", nor may "ColorBrewer" appear in their name,
% without prior written permission of Cynthia Brewer.  THIS SOFTWARE IS
% PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT
% NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL CYNTHIA BREWER,
% MARK HARROWER, OR THE PENNSYLVANIA STATE UNIVERSITY BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
% OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
% SUCH DAMAGE.
%

cs = colorBrewScales;

%colorNames = unique({cs(:).colorName});
%colorTypes = unique({cs(:).colorType});
%numColors  = unique([cs(:).numColors]);

if exist('ctype','var'),
  if ~isempty(ctype),
    typeIndex = strmatch(ctype,{cs(:).colorType},'exact');
    if isempty(typeIndex),
      msg = sprintf('failed to find type: %s', ctype);
      error(msg)
    else
      % sub-select only 'ctype' entries
      cs = cs(typeIndex);
    end
  end
end

if exist('cname','var'),
  if ~isempty(cname),
    nameIndex = strmatch(cname,{cs(:).colorName},'exact');
    if isempty(nameIndex),
      msg = sprintf('failed to find name: %s', cname);
      error(msg)
    else
      % sub-select only 'name' entries
      cs = cs(nameIndex);
    end
  end
end

if exist('ncolor','var'),
  if ~isempty(ncolor),
    ncolorIndex = find([cs(:).numColors] == ncolor);
    if isempty(ncolorIndex),
      msg = sprintf('failed to find ncolor: %d\n', ncolor);
      msg = [msg, sprintf('possible values are: %d\n', unique([cs(:).numColors]))];
      error(msg)
    else
      % sub-select only 'name' entries
      cs = cs(ncolorIndex);
    end
  end
end

return

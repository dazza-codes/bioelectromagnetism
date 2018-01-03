#!/usr/bin/env python

# Copyright (c) 2002 Cynthia Brewer, Mark Harrower, and The Pennsylvania
# State University. All rights reserved.  Redistribution and use in source
# and binary forms, with or without modification, are permitted provided
# that the following conditions are met: 1. Redistributions as source code
# must retain the above copyright notice, this list of conditions and the
# following disclaimer.  2. The end-user documentation included with the
# redistribution, if any, must include the following acknowledgment: This
# product includes color specifications and designs developed by Cynthia
# Brewer (http://colorbrewer.org/).  Alternately, this acknowledgment may
# appear in the software itself, if and wherever such third-party
# acknowledgments normally appear.  4. The name "ColorBrewer" must not be
# used to endorse or promote products derived from this software without
# prior written permission. For written permission, please contact Cynthia
# Brewer at cbrewer@psu.edu.  5. Products derived from this software may not
# be called "ColorBrewer", nor may "ColorBrewer" appear in their name,
# without prior written permission of Cynthia Brewer.  THIS SOFTWARE IS
# PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL CYNTHIA BREWER,
# MARK HARROWER, OR THE PENNSYLVANIA STATE UNIVERSITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.

import sys
import os
import re


if '-matlab' in sys.argv:
    printMatlab = 1
else:
    printMatlab = 0

if '-verify' in sys.argv:
    verify = 1
else:
    verify = 0


class colorScale:

    def __init__(self):
        self.colorName = ''
        self.numColors = int()
        self.colorType = ''
        self.critVal = float()
        self.colorNum = []
        self.colorLet = []
        self.colorRGB = []
    
    def __str__(self):
        string = "%-18s %2d %-4s" % (self.colorName, self.numColors, self.colorType)
        if self.critVal > 0:
            if (self.critVal % 1) > 0:
                string += " %7.1f" % self.critVal
            else:
                string += " %7d" % int(self.critVal)
        else:
            string += " %7s" % ''
        string += " %8d" % self.colorNum[0]
        string += " %-11s" % self.colorLet[0]
        string += " %3d %3d %3d\n" % tuple(self.colorRGB[0])
        for i in range(1,self.numColors):
            string += " %42d" % self.colorNum[i]
            string += " %-11s" % self.colorLet[i]
            string += " %3d %3d %3d" % tuple(self.colorRGB[i])
            if i < (self.numColors - 1): string += "\n"
        return string
    
    def matlab(self):
        "Return a string that will create a matlab struct array"
        string  = "cs(i).colorName = '%s';\n" % self.colorName
        string += "cs(i).numColors = %d;\n" % self.numColors
        string += "cs(i).colorType = '%s';\n" % self.colorType
        string += "cs(i).critVal = %3.1f;\n" % self.critVal
        string += "cs(i).colorNum = %s;\n" % str(self.colorNum)
        string += "cs(i).colorLet = %s;\n" % str(self.colorLet)
        string += "cs(i).colorRGB = [%3d %3d %3d" % tuple(self.colorRGB[0])
        for i in range(1, self.numColors):
            string += ";\n%17s %3d %3d %3d" % (('',) + tuple(self.colorRGB[i]))
        string += "] / 255;\n"
        return string


# open and read the color brewer data file
colorBrewFileName = 'colorBrewScales.txt'
colorBrewFileName = os.path.join(os.getcwd(), colorBrewFileName)
colorBrewFile = file(colorBrewFileName,"r")
colorBrewLines = colorBrewFile.readlines()
colorBrewFile.close()

lines = [s.rstrip() for s in colorBrewLines]

# find all the comment lines in the data file
commentRE = re.compile('^#')           
commentLinesN = []
for lineN in range(len(lines)):
    if commentRE.findall(lines[lineN]):
        commentLinesN.append(lineN)

# find all the line numbers with a new color scale
colorNameLineRE = re.compile('^\w')     # any alphanum at start
scaleLinesN = []
for lineN in range(len(lines)):
    if commentRE.findall(lines[lineN]):
        pass
    elif colorNameLineRE.findall(lines[lineN]):
        scaleLinesN.append(lineN)


if verify:
    for lineN in commentLinesN:
        print lines[lineN]

for lineN in scaleLinesN:
    line = lines[lineN]
    lineSplit = line.split()
    
    # extract data from the first line of the new color scale
    cs = colorScale()
    cs.colorName = lineSplit[0]
    cs.numColors = int(lineSplit[1])
    cs.colorType = lineSplit[2]
    if len(lineSplit) > 8: cs.critVal = float(lineSplit[3])
    cs.colorNum.append(int(lineSplit[-5]))
    cs.colorLet.append(lineSplit[-4])
    cs.colorRGB.append([int(d) for d in lineSplit[-3:]])
    
    # read the rest of the data for this color scale
    for i in range(1,cs.numColors):
        line = lines[lineN+i]
        lineSplit = line.split()
        cs.colorNum.append(int(lineSplit[-5]))
        cs.colorLet.append(lineSplit[-4])
        cs.colorRGB.append([int(d) for d in lineSplit[-3:]])
    
    if verify:
        print cs
    
    if printMatlab:
        print "i = %d;" % printMatlab
        print cs.matlab()
        printMatlab += 1



suptitle.m
Original suptitle file.

suptitle_withpatch.m
Updated suptitle file. Bug occurs when suptitle is called, if the current axes has a legend, it is hidden.  
This is because SUPTITLE saves the current axes (as |haold|) on entry and calls |axes(haold)| on exit, hiding 
the legend. Patch which fixes this is in this m-file. These changes were made by Mark Histed (histed@mit.edu).


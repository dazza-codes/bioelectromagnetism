Talairach Daemon Client v 2.0
http://ric.uthscsa.edu/resources
Programmers: Peter Kochunov, Angela Uecker
Copyright The Research Imaging Center - UTHSCSA

The Talairach Daemon Client will read tab or space delimited records from text files containing lists of Talairach coordinates arranged in x-y-z order. Or, using the Single Point Processing dialog, one can input in a single coordinate to label. It will then look up the coordinate in the Talairach Daemon database for the Talairach label. There are options to search for the single point, search range or nearest gray matter. The output is written to a file which can be viewed in the program, a third-party text editor or imported into a third-party spreadsheet.

Installation
Windows:
IMPORTANT. The TD Java search client requires the Microsoft Java virtual machine. The download instructions show you how to get the most recent version. The Java TD search client will not work with others such as Sun's Java virtual machine.

Download the software at http://ric.uthscsa.edu/resources.

Install by double-clicking on the installer file, tdc.pkg.
NOTE: The downloaded files will be installed in the directory "C:tdc" by default. The executable file in this directory is "TDC.exe", and you can make a shortcut and place it on the desktop for more convenient usage.

Macintosh:
Download the software at http://ric.uthscsa.edu/resources.

Double-click on the tdc.pkg file and following the installation instructions that follow.
NOTE: The default installation directory is the system Applications directory. You can also choose your installation directory during the installation process.

Solaris:
Download the software at http://ric.uthscsa.edu/resources.

Unzip the tdc.zip file.  To run the program, type "tdc" at the command-line.


Batch File Format
Before you begin searching you must create a text file with a list of the x-y-z Talairach coordinates that you want to label. The file "test.txt" is an example of the required format. The current version of the TD Java Search client will read files containing up to 40,000 records (x-y-z coordinates).

General Instructions
Batch Processing:

Once file(s) are created for batch processing, select your file(s) from the File menu. The number of coordinates will be printed in the information window.
Press the "Process" button to begin. Several messages will be printed in the information window to inform you of the file loading progress. A progress bar is updated every 10th coordinate to give you an idea of how long it will take to retrieve all of the labels. As each label is retrieved it is appended to the output file along with the Talairach coordinates. The output file will be created dependent on the input file name (i.e. If your input file is test.txt, your output file will be test.td. Note the different extensions.)
When all labels have been retrieved you are given the option to view the result file in the program's text area. Also, the result file can be viewed by any text editor or import it into a spreadsheet.

Single Point Processing:

From the Option menu, choose "Single Point Processing". A dialog will appear.
Enter in your coordinates into the fields for x, y, and z.
Choose your options for processing (i.e. Single Point, Range Search or Nearest Gray Matter)
Choose your options for saving to a file. If you want to save to a file, you must check the checkbox and pick a file. The file will stay open for writing as long as the Single Point Processing dialog is open. Once it closes, the file is closed. If the dialog is opened again for processing, you would need to pick another file for the next set of processing or the program will write over the current file.
Click Process to process the point.

Suggested Uses
Brodmann Area (BA) Labels for Cortical Activations:

When seeking Brodmann Area labels for cortical activation sites it is possible to extend the search diameter to find the nearest BA label. For an experiment designed to activate the M1 mouth motor region only 38% of the sites were found to fall within Brodmann Areas, but as the search diameter was increased to 3 mm, 5 mm, and 7 mm, BA labels were obtained for 62%, 92%, and 100% of the sites. This example shows the utility of obtaining BA labels for cortical activation sites. While the appropriate search range may vary from site to site within the brain, this is easy to test. The current release of the Java TD client only supports the 5 mm search diameter, but ranges from 3-11 mm will be made available soon. As the search range increases the number of labels found increases. This presents two problemsto the user. First, the label retrieval process becomes much slower. For a 5x5x5 mm search range 125 voxels are searched and many labels are found. The TD client reduces the number of labels by only responding with unique labels. The unique labels are organized by incidence within the search range, with the highest label incidence being the first label returned. Along with each of the unique labels is the number of voxel within the search range with that label. Second, with large search ranges, the user has to deal with much larger files to sort through. The recommended strategy is to search with a small search region initially, remove coordinates with BA labels, and proceed using smaller input files for larger search ranges.

Anatomical Organization of Coordinate Data:

For application that provide Talairach coordinates (SPM, MEDx, etc.), it is helpful to use anatomical labels from the TD to organize findings anatomically. The data saved by the Java Talairach Daemon client can be used to create a labeled file. A common use of this labeled file is to open it in Excel and rearrange the original slice-ordered data into 3-D anatomical groupings by sorting by lobe and gyrus. This provides users with a good sense of which activation sites to group together anatomically.

Feedback
We welcome your comments regarding the layout and usage of the dialogs and windows. Also, should the program not work correctly send a bug report. 

Thank you,
The Talairach Daemon Team

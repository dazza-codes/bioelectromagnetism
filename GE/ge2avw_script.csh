
#!/bin/csh

#PATH PARAMETERS
set airPath = /data/fmri/bin/AIR3.08

#DEFINE THE PARAMETERS FOR THE ANALYSIS
# MR header size is 7904

echo "exam number ?" 
@ studyNum = $<
echo "series number ?"
@ seriesNum = $<
echo "number of slices"
@ numSlices = $<
echo "matrix size (e.g. 128) -- [assumes square matrix]"
@ matrix = $<

#CALCULATE SLICE SIZE (IN BYTES) - MINUS THE GE HEADER
@ imgSize = (2 * $matrix * $matrix)

echo Uncompressing files...
uncompress *.Z
gunzip *.gz

ls E${studyNum}S${seriesNum}I*.MR | awk -F"." '{print substr($1,10,3)}' > tempList
sort -n tempList > sortedTempList
@ first = `head -1 sortedTempList`
/bin/rm tempList
/bin/rm sortedTempList

#EXTRACT VOXEL DIMENSIONS FROM GE HEADER FILE
set sliceThickness = `thickness E${studyNum}S${seriesNum}I${first}.MR`
set pixX = `pixelsize E${studyNum}S${seriesNum}I${first}.MR x`
set pixY = `pixelsize E${studyNum}S${seriesNum}I${first}.MR y`

#.MR FILE CONVERSION TO 3D .IMG
echo Converting series number $seriesNum
set x = $first
while ($x <= $numSlices)
         tail -${imgSize}c E${studyNum}S${seriesNum}I${x}.MR > ${x}temp.img
         $airPath/makeaheader ${x}temp.img 1 $matrix $matrix 1 $pixX $pixY $sliceThickness
         echo ${x}temp.img >> imgVol
         @ x++
end
set volContents = `cat imgVol`  
/bin/rm imgVol
$airPath/reunite E${studyNum}S${seriesNum}_AIR3D y $volContents
/bin/rm *temp.*
$airPath/reorient E${studyNum}S${seriesNum}_AIR3D E${studyNum}S${seriesNum}_AIR3Dspm o xx yy x
/bin/rm E${studyNum}S${seriesNum}_AIR3D.*

exit

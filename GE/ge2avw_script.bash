#!/bin/bash

echo ""

#PATH PARAMETERS
air=`which reslice | sed s#reslice##`
if [ $air = '' ]; then
    echo 'Cannot locate AIR functions - install it or modify this script'
    exit
else
    airPath=$air
    echo "Found AIR functions in $airPath"
fi


#DEFINE THE PARAMETERS FOR THE ANALYSIS
# MR header size is 7904

#echo usage, if no arguments given
if [ $# == 0 ]; then
  echo "No Input Parameters..."
  echo ""
  echo ""
  echo "ge2avw_script -e examN -s seriesN -nslices Nslices -m matrix_size"
  echo " "
  echo "examN, exam number"
  echo "seriesN, series number"
  echo "Nslices, number of slices"
  echo "matrix_size (e.g. 128) -- [assumes square matrix]"
  echo ""
  echo "This script employs AIR and shell utilities to extract"
  echo "GE data into analyze image files."
  echo ""
  exit
fi

#init
examNum=1
seriesNum=1
numSlices=1
matrix=""

#collect arguments
#ge2avw_script -e examN -s seriesN -nslices Nslices -m matrix_size
while [ $# != 0 ]; do
    case "$1" in
	-e)
	    shift; examNum=$1; shift
	    ;;
	
	-s)
	    shift; seriesNum=$1; shift
	    ;;
	
	-nslices)
	    shift; numSlices=$1; shift
	    ;;
	
	-m)
	    shift; matrix=$1; shift
	    ;;
	
	*)
	    echo "Unknown option '$1'"
	    exit 1
	    ;;
	
    esac
done

echo "processing exam $examNum, series $seriesNum, $numSlices slices, matrix $matrix"

#CALCULATE SLICE SIZE (IN BYTES) - MINUS THE GE HEADER
imgSize=`echo "2 * $matrix * $matrix" | bc`
echo "Image size is $imgSize bytes"

if   [ -e *.Z ]; then
    echo Uncompressing files...
    uncompress *.Z
elif [ -e *.gz ]; then
    echo Uncompressing files...
    gunzip *.gz
fi

path=`pwd`

files=`ls E${examNum}S${seriesNum}I*.MR`
case "$files" in
    "")
	echo "No files found, check examN & seriesN"
	exit
	;;
esac

echo "$files" | awk -F"." '{print substr($1,10,3)}' > tempList

sort -n tempList > sortedTempList
first=`head -1 sortedTempList`
/bin/rm -f tempList sortedTempList




#EXTRACT VOXEL DIMENSIONS FROM GE HEADER FILE
sliceThickness=`thickness E${examNum}S${seriesNum}I${first}.MR`
pixX=`pixelsize E${examNum}S${seriesNum}I${first}.MR x`
pixY=`pixelsize E${examNum}S${seriesNum}I${first}.MR y`

if [ "$pixX" = "" ]; then
    echo "pixX is not defined"
    exit
fi



#.MR FILE CONVERSION TO 3D .IMG
echo Converting series number $seriesNum
set x = $first
while ($x <= $numSlices)
         tail -${imgSize}c E${examNum}S${seriesNum}I${x}.MR > ${x}temp.img
         $airPath/makeaheader ${x}temp.img 1 $matrix $matrix 1 $pixX $pixY $sliceThickness
         echo ${x}temp.img >> imgVol
         @ x++
end
/bin/rm -f *temp.*

volContents=`cat imgVol`; /bin/rm -f imgVol

$airPath/reunite E${examNum}S${seriesNum}_AIR3D y $volContents
$airPath/reorient E${examNum}S${seriesNum}_AIR3D E${examNum}S${seriesNum}_AIR3Dspm o xx yy x

/bin/rm -f E${examNum}S${seriesNum}_AIR3D.*

exit

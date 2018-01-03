#

#   bet.tcl - GUI for BET - Brain Extraction Tool
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 1999-2001 University of Oxford
#
#   TCLCOPYRIGHT

#{{{ setups

source [ file dirname [ info script ] ]/fslstart.tcl

set VARS(history) {}

#}}}
#{{{ bet

proc bet { w } {

    #{{{ vars and setup

global bet FSLDIR INMEDX argc argv PWD

toplevel $w
wm title $w "BET - Brain Extraction Tool - v1.1"
wm iconname $w "BET"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

frame $w.f

#}}}
    #{{{ input image (if not in medx)

if { $INMEDX } {
    set bet($w,input) "dummy"
    set bet($w,output) "dummy"
} else {
    if { $argc > 0 && [ string length [ lindex $argv 0 ] ] > 0 } {
        set inputname [ file rootname [ lindex $argv 0 ] ].hdr
        if { [ string first / $inputname ] == 0 || [ string first ~ $inputname ] == 0 } {
            set bet($w,input) $inputname
        } else {
            set bet($w,input) ${PWD}/$inputname
        }
        set bet($w,output) [ file rootname $bet($w,input) ]_brain.hdr
    }

    FSLFileEntry $w.f.input \
            -variable bet($w,input) \
            -pattern "*.hdr" \
            -directory $PWD \
            -label "Input image   " \
            -title "Select the input image" \
            -width 50 \
            -filterhist VARS(history) \
	    -command "bet:select $w"

    FSLFileEntry $w.f.output \
            -variable bet($w,output) \
            -pattern "*.hdr" \
            -directory $PWD \
            -label "Output image" \
            -title "Select the output image" \
            -width 50 \
            -filterhist VARS(history)

    pack $w.f.input $w.f.output -in $w.f -side top -padx 5 -pady 5 -anchor w
}

#}}}
    #{{{ generate segmented image

frame $w.f.segment

label $w.f.segment.label -text "Generate image with non-brain matter removed"

set bet($w,segment_yn) 1
checkbutton $w.f.segment.yn -variable bet($w,segment_yn)

pack $w.f.segment.label $w.f.segment.yn -in $w.f.segment -side left -padx 5

#}}}
    #{{{ generate overlay image

frame $w.f.overlay

label $w.f.overlay.label -text "Generate image with estimated brain surface overlaid on original"

set bet($w,overlay_yn) 0
checkbutton $w.f.overlay.yn -variable bet($w,overlay_yn)

pack $w.f.overlay.label $w.f.overlay.yn -in $w.f.overlay -side left -padx 5

#}}}
    #{{{ Advanced Options

collapsible frame $w.f.opts -title "Advanced Options"
set bet($w,xtopol_yn) 0
set bet($w,cost_yn) 0

#{{{ generate binary brain image

frame $w.f.opts.mask

label $w.f.opts.mask.label -text "Generate binary brain mask image"

set bet($w,mask_yn) 0
checkbutton $w.f.opts.mask.yn -variable bet($w,mask_yn)

pack $w.f.opts.mask.label $w.f.opts.mask.yn -in $w.f.opts.mask -side left -padx 5

#}}}
#{{{ apply thresholding

frame $w.f.opts.threshold

label $w.f.opts.threshold.label -text "Apply thresholding to segmented brain image (and mask if required)"

set bet($w,threshold_yn) 0
checkbutton $w.f.opts.threshold.yn -variable bet($w,threshold_yn)

pack $w.f.opts.threshold.label $w.f.opts.threshold.yn -in $w.f.opts.threshold -side left -padx 5

#}}}
#{{{ generate skull image

frame $w.f.opts.skull

label $w.f.opts.skull.label -text "Generate exterior skull surface image"
set bet($w,skull_yn) 0
checkbutton $w.f.opts.skull.yn -variable bet($w,skull_yn) -command "bet:updateskullcombine $w"

frame $w.f.opts.combine_skull
label $w.f.opts.combine_skull.label -text "Combine skull image with original"
set bet($w,combine_skull_yn) 0
checkbutton $w.f.opts.combine_skull.yn -variable bet($w,combine_skull_yn)
pack $w.f.opts.combine_skull.label $w.f.opts.combine_skull.yn -in $w.f.opts.combine_skull -side left -padx 5

pack $w.f.opts.skull.label $w.f.opts.skull.yn -in $w.f.opts.skull -side left -padx 5

#}}}
#{{{ COMMENT generate raw cost image

#frame $w.f.opts.cost
#
#label $w.f.opts.cost.label -text "Generate raw cost function image"
#
#set bet($w,cost_yn) 0
#checkbutton $w.f.opts.cost.yn -variable bet($w,cost_yn)
#
#pack $w.f.opts.cost.label $w.f.opts.cost.yn -in $w.f.opts.cost -side left -padx 5
#

#}}}
#{{{ COMMENT generate XTopol output

#frame $w.f.opts.xtopol
#
#label $w.f.opts.xtopol.label -text "Generate output for viewing in XTopol surface viewer"
#
#set bet($w,xtopol_yn) 0
#checkbutton $w.f.opts.xtopol.yn -variable bet($w,xtopol_yn)
#
#pack $w.f.opts.xtopol.label $w.f.opts.xtopol.yn -in $w.f.opts.xtopol -side left -padx 5
#

#}}}
#{{{ fractional brain threshold

set bet($w,fraction) 0.5
tixControl $w.f.opts.fraction -label "Fractional intensity threshold; smaller values give larger brain outline estimates" \
        -variable bet($w,fraction) -step 0.05 -min 0.0001 -max 0.9999 -selectmode immediate

#}}}
#{{{ gradient brain threshold

set bet($w,gradient) 0
tixControl $w.f.opts.gradient -label "Threshold gradient; positive values give larger brain outline at bottom, smaller at top" \
        -variable bet($w,gradient) -step 0.05 -min -1 -max 1 -selectmode immediate

#}}}

#pack $w.f.opts.mask $w.f.opts.threshold $w.f.opts.skull $w.f.opts.cost $w.f.opts.xtopol $w.f.opts.fraction $w.f.opts.gradient -in $w.f.opts.b -side top -anchor w -pady 0
pack $w.f.opts.mask $w.f.opts.threshold $w.f.opts.skull $w.f.opts.fraction $w.f.opts.gradient -in $w.f.opts.b -side top -anchor w -pady 0

#}}}

    pack $w.f.segment $w.f.overlay $w.f.opts -in $w.f -side top -padx 5 -pady 0 -anchor w

    #{{{ Button Frame

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
 
    button $w.ok \
        -text "OK" -width 5 \
        -command "bet:apply $w destroy"
    bind $w.ok <Return> {
        [winfo toplevel %W].ok invoke}
 
    button $w.apply     -command "bet:apply $w keep" \
        -text "Apply" -width 5
    bind $w.apply <Return> {
        [winfo toplevel %W].apply invoke}
 
    button $w.cancel    -command "bet:destroy $w" \
        -text "Cancel" -width 5
    bind $w.cancel <Return> {
        [winfo toplevel %W].cancel invoke}
 
    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/bet/index.html" \
            -text "Help" -width 5
    bind $w.help <Return> {
        [winfo toplevel %W].help invoke}

    pack $w.btns.b -side bottom -fill x
    pack $w.ok $w.apply $w.cancel $w.help -in $w.btns.b \
        -side left -expand yes -padx 3 -pady 10 -fill y
 
    pack $w.f $w.btns -expand yes -fill both

#}}}
}

#}}}
#{{{ bet:select

proc bet:select { w dummy } {

    global bet

    set bet($w,output) [ file rootname $bet($w,input) ]_brain.hdr

#    if { [ string length $bet($w,output) ] == 0 } {
#	set bet($w,output) [ file rootname $bet($w,input) ]_brain.hdr
#    }
}

#}}}
#{{{ bet:apply

proc bet:apply { w dialog } {
    global bet

    bet_proc $bet($w,input) $bet($w,output) $bet($w,segment_yn) $bet($w,overlay_yn) $bet($w,mask_yn) $bet($w,threshold_yn) $bet($w,xtopol_yn) $bet($w,cost_yn) $bet($w,skull_yn) $bet($w,combine_skull_yn) $bet($w,fraction) $bet($w,gradient)

    update idletasks

    if {$dialog == "destroy"} {
        bet:destroy $w
    }
}

#}}}
#{{{ bet:destroy

# Summary:      Destroys bet dialog box
proc bet:destroy { w } {
    destroy $w
}

#}}}
#{{{ bet:updateskullcombine

proc bet:updateskullcombine { w } {
    global bet INMEDX

    if { $bet($w,skull_yn) == 0 || $INMEDX == 0 } {
	set bet($w,combine_skull_yn) 0
	pack forget $w.f.opts.combine_skull
    } else {
	pack $w.f.opts.combine_skull -in $w.f.opts.skull -side left -padx 5
    }
}

#}}}

wm withdraw .
bet .rename
tkwait window .rename


DEFINE !ERPOPEN (	 Fin  = !TOKENS(1)
			/Fout = !TOKENS(2)	).

GET DATA  /TYPE = TXT
 /FILE = !Fin
 /DELCASE = LINE
 /DELIMITERS = "\t"
 /ARRANGEMENT = DELIMITED
 /FIRSTCASE = 2
 /IMPORTCASE = ALL
 /VARIABLES =
 Subjects   A15
 Group      F8.2
 L_IPFdif   F10.2
 L_IFdif    F10.2
 L_SPFdif   F10.2
 L_SFdif    F10.2
 L_SCdif    F10.2
 L_ICdif    F10.2
 L_SPdif    F10.2
 L_IPdif    F10.2
 L_ATdif    F10.2
 L_PTdif    F10.2
 L_OTdif    F10.2
 R_IPFdif   F10.2
 R_IFdif    F10.2
 R_SPFdif   F10.2
 R_SFdif    F10.2
 R_SCdif    F10.2
 R_ICdif    F10.2
 R_SPdif    F10.2
 R_IPdif    F10.2
 R_ATdif    F10.2
 R_PTdif    F10.2
 R_OTdif    F10.2
 X          1X
 .
CACHE.
EXECUTE.

VARIABLE LEVEL group
 l_ipfdif l_ifdif  l_spfdif l_sfdif  l_scdif  l_icdif
 l_spdif  l_ipdif  l_atdif  l_ptdif  l_otdif  r_ipfdif
 r_ifdif  r_spfdif r_sfdif  r_scdif  r_icdif  r_spdif
 r_ipdif  r_atdif  r_ptdif  r_otdif (SCALE) .

VALUE LABELS group 1 'CONT' 2 'PTSD' .

SAVE OUTFILE = !Fout
  /COMPRESSED.

!ENDDEFINE.

!ERPOPEN Fin ='D:\data_emse\ptsdpet\link14hz\regions\eeg_regions_dif_n0300_amp.txt'
         Fout='D:\data_emse\ptsdpet\link14hz\regions\eeg_regions_dif_n0300_amp.sav' .
!ERPOPEN Fin ='D:\data_emse\ptsdpet\link14hz\regions\eeg_regions_dif_n0300_lat.txt'
         Fout='D:\data_emse\ptsdpet\link14hz\regions\eeg_regions_dif_n0300_lat.sav' .

!ERPOPEN Fin ='D:\data_emse\ptsdpet\link14hz\regions\eeg_regions_dif_p0550_amp.txt'
         Fout='D:\data_emse\ptsdpet\link14hz\regions\eeg_regions_dif_p0550_amp.sav' .
!ERPOPEN Fin ='D:\data_emse\ptsdpet\link14hz\regions\eeg_regions_dif_p0550_lat.txt'
         Fout='D:\data_emse\ptsdpet\link14hz\regions\eeg_regions_dif_p0550_lat.sav' .


script file="d:\mydocuments\programming\spss_scripts\SaveClose.sbs"
 ("D:\data_emse\ptsdpet\link14hz\regions\regions_dif_load.spo") .

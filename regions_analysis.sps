

DEFINE !ERPDES ().

SORT CASES BY group .
SPLIT FILE
  SEPARATE BY group .

DESCRIPTIVES
  VARIABLES=l_ipfexp l_ifexp l_spfexp l_sfexp l_scexp l_icexp l_spexp l_ipexp
  l_atexp l_ptexp l_otexp r_ipfexp r_ifexp r_spfexp r_sfexp r_scexp r_icexp
  r_spexp r_ipexp r_atexp r_ptexp r_otexp l_ipfcon l_ifcon l_spfcon l_sfcon
  l_sccon l_iccon l_spcon l_ipcon l_atcon l_ptcon l_otcon r_ipfcon r_ifcon
  r_spfcon r_sfcon r_sccon r_iccon r_spcon r_ipcon r_atcon r_ptcon r_otcon
  /STATISTICS=MEAN STDDEV VARIANCE MIN MAX SEMEAN SKEWNESS
  /SORT=MEAN (D) .

SPLIT FILE
  OFF.

!ENDDEFINE.



DEFINE !ERPGLM (	 V1 = !TOKENS(1)
			/V2 = !TOKENS(1)
			/V3 = !TOKENS(1)
			/V4 = !TOKENS(1)	).

GLM
  !V1 !V2 !V3 !V4 BY group
  /WSFACTOR = wm 2 Polynomial hemi 2 Polynomial
  /CONTRAST (group)=Simple(1)
  /METHOD = SSTYPE(3)
  /EMMEANS = TABLES(group)
  /EMMEANS = TABLES(wm)
  /EMMEANS = TABLES(hemi)
  /EMMEANS = TABLES(group*wm)
  /EMMEANS = TABLES(group*hemi)
  /EMMEANS = TABLES(wm*hemi)
  /EMMEANS = TABLES(group*wm*hemi)
  /PLOT = PROFILE( wm*group )
  /PRINT = DESCRIPTIVE ETASQ OPOWER PARAMETER HOMOGENEITY
  /PLOT = RESIDUALS
  /CRITERIA = ALPHA(.05)
  /WSDESIGN = wm hemi wm*hemi
  /DESIGN = group .

!ENDDEFINE.


DEFINE !ERPPLOT ().

GRAPH
  /ERRORBAR( CI 95 )=l_ipfexp l_ifexp l_spfexp l_sfexp l_scexp l_icexp
  l_spexp l_ipexp l_atexp l_ptexp l_otexp r_ipfexp r_ifexp r_spfexp r_sfexp
  r_scexp r_icexp r_spexp r_ipexp r_atexp r_ptexp r_otexp l_ipfcon l_ifcon
  l_spfcon l_sfcon l_sccon l_iccon l_spcon l_ipcon l_atcon l_ptcon l_otcon
  r_ipfcon r_ifcon r_spfcon r_sfcon r_sccon r_iccon r_spcon r_ipcon r_atcon
  r_ptcon r_otcon
  /MISSING=VARIABLEWISE .

!ENDDEFINE.


DEFINE !ERPRUN ( Fin  = !TOKENS(1)
		    /Fout = !TOKENS(1)  ).

GET
  FILE = !Fin .
!ERPDES .
!ERPPLOT .
!ERPGLM	V1=l_ipfexp	V2=r_ipfexp	V3=l_ipfcon	V4=r_ipfcon	.
!ERPGLM	V1=l_ifexp	V2=r_ifexp	V3=l_ifcon	V4=r_ifcon	.
!ERPGLM	V1=l_spfexp	V2=r_spfexp	V3=l_spfcon	V4=r_spfcon	.
!ERPGLM	V1=l_sfexp	V2=r_sfexp	V3=l_sfcon	V4=r_sfcon	.
!ERPGLM	V1=l_scexp	V2=r_scexp	V3=l_sccon	V4=r_sccon	.
!ERPGLM	V1=l_icexp	V2=r_icexp	V3=l_iccon	V4=r_iccon	.
!ERPGLM	V1=l_spexp	V2=r_spexp	V3=l_spcon	V4=r_spcon	.
!ERPGLM	V1=l_ipexp	V2=r_ipexp	V3=l_ipcon	V4=r_ipcon	.
!ERPGLM	V1=l_atexp	V2=r_atexp	V3=l_atcon	V4=r_atcon	.
!ERPGLM	V1=l_ptexp	V2=r_ptexp	V3=l_ptcon	V4=r_ptcon	.
!ERPGLM	V1=l_otexp	V2=r_otexp	V3=l_otcon	V4=r_otcon	.

script file="d:\mydocuments\programming\spss_scripts\SaveClose.sbs"
 (!Fout) .

!ENDDEFINE.




!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0040_amp.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0040_amp.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0040_lat.sav' .
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0040_lat.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0040_amp.sav' .
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0040_amp.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0040_lat.sav' .
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0040_lat.spo' .

!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0100_amp.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0100_amp.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0100_lat.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0100_lat.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0100_amp.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0100_amp.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0100_lat.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0100_lat.spo' .

!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0220_amp.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0220_amp.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0220_lat.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0220_lat.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0220_amp.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0220_amp.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0220_lat.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0220_lat.spo' .

!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0320_amp.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0320_amp.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0320_lat.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0320_lat.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0320_amp.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0320_amp.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0320_lat.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0320_lat.spo' .

!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0450_amp.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0450_amp.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0450_lat.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_n0450_lat.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0450_amp.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0450_amp.spo' .
!ERPRUN Fin ='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0450_lat.sav'
        Fout='D:\data_emse\ptsdpet\link14hz\regions\wm_regions_p0450_lat.spo' .

function rdb_hdr = GE_readRawHeaderRdbRec(fid)
%
% rdb_header = GE_readRawHeaderRdbRec(fid)
%
% Loads the rdb_header located in a file with file id fid
% and returns it as a structure. 
%
% Souheil J. Inati
% Dartmouth College
% August 2000
% souheil.inati@dartmouth.edu

% define the structure and read in the data
% to overcome the line length limit
% break up the assignment into pieces using the setfield function
rdb_hdr = struct('rdb_hdr_rdbm_rev', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_run_int', fread(fid,1,'int32')); % Rdy pkt Run Number %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_scan_seq', fread(fid,1,'int16')); % Rdy pkt Sequence Number %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_run_char', fread(fid,6,'char')); % Rdy pkt Run no in fread(fid,1,'char') %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_scan_date', fread(fid,10,'char')); %%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_scan_time', fread(fid,8,'char')); %%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_logo', fread(fid,10,'char')); % rdbmused to verify file %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_file_contents', fread(fid,1,'int16')); % Data type 0=emp 1=nrec 2=rw 	0, 1, 2 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_lock_mode', fread(fid,1,'int16')); % unused %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_dacq_ctrl', fread(fid,1,'int16')); % rhdacqctrl bit mask		15 bits %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_recon_ctrl', fread(fid,1,'int16')); % rhrcctrl bit mask 		15 bits %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_exec_ctrl', fread(fid,1,'int16')); % rhexecctrl bit mask 		15 bits %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_scan_type', fread(fid,1,'int16')); % bit mask 			15 bits %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_data_collect_type', fread(fid,1,'int16')); % rhtypebit mask		15 bits %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_data_format', fread(fid,1,'int16')); % rhformatbit mask 		15 bits %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_recon', fread(fid,1,'int16')); % rhrecon proc-a-son recon	0 - 100 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_datacq', fread(fid,1,'int16')); % rhdatacq proc-a-son dacq %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_npasses', fread(fid,1,'int16')); % rhnpasses passes for a scan0 - 256 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_npomp', fread(fid,1,'int16')); % rhnpomp pomp group slices	1,2 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_nslices', fread(fid,1,'int16')); % rhnslices slices in a pass	0 - 256 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_nechoes', fread(fid,1,'int16')); % rhnecho echoes of a slice	1 - 32 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_navs', fread(fid,1,'int16')); % rhnavs num of excitiations	1 - 32727 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_nframes', fread(fid,1,'int16')); % rhnframes yres 0 - 1024 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_baseline_views', fread(fid,1,'int16')); % rhbline baselines 0 - 1028 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_hnover', fread(fid,1,'int16')); % rhhnover overscans 0 - 1024 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_frame_size', fread(fid,1,'int16')); % rhfrsize xres 0 - 1024 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_point_size', fread(fid,1,'int16')); % rhptsize 2 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vquant', fread(fid,1,'int16')); % rhvquant 3d volumes	1 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_cheart', fread(fid,1,'int16')); % RX Cine heart phases 1 - 32 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ctr', fread(fid,1,'float32')); % RX Cine TR in sec	0 - 3.40282e38%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ctrr', fread(fid,1,'float32')); % RX Cine RR in sec	0 - 30.0 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_initpass', fread(fid,1,'int16')); % rhinitpass allocate passes0 - 32767 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_incrpass', fread(fid,1,'int16')); % rhincrpass tps autopauses	0 - 32767 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_method_ctrl', fread(fid,1,'int16')); % rhmethod0=recon, 1=psd	0, 1 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_da_xres', fread(fid,1,'int16')); % rhdaxres 0 - 1024 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_da_yres', fread(fid,1,'int16')); % rhdayres 0 - 2049 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_rc_xres', fread(fid,1,'int16')); % rhrcxres 0 - 1024 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_rc_yres', fread(fid,1,'int16')); % rhrcyres 0 - 1024 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_im_size', fread(fid,1,'int16')); % rhimsize 0 - 512 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_rc_zres', fread(fid,1,'int32')); % power of 2 > rhnslices	0 - 128 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_raw_pass_size', fread(fid,1,'int32')); % rhrawsize 0 - 2147483647%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_sspsave', fread(fid,1,'int32')); % rhsspsave 0 - 2147483647%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_udasave', fread(fid,1,'int32')); % rhudasave 0 - 2147483647%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_fermi_radius', fread(fid,1,'float32')); % rhfermr fermi radius		0 - 3.40282e38%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_fermi_width', fread(fid,1,'float32')); % rhfermw fermi width		0 - 3.40282e38%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_fermi_ecc', fread(fid,1,'float32')); % rhferme fermi excentiricty	0 - 3.40282e38%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_clip_min', fread(fid,1,'float32')); % rhclipmin 4x IP limit		+-16383 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_clip_max', fread(fid,1,'float32')); % rhclipmax 4x IP limit		+-16383 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_default_offset', fread(fid,1,'float32')); % rhdoffset default offset = 0	+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_xoff', fread(fid,1,'float32')); % rhxoff scroll img in x 	+-256 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_yoff', fread(fid,1,'float32')); % rhyoff scroll img in y	+-256 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_nwin', fread(fid,1,'float32')); % rhnwin hecho window width	0 - 256 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ntran', fread(fid,1,'float32')); % rhntran hecho trans width	0 - 256 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_scalei', fread(fid,1,'float32')); % PS rhscalei			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_scaleq', fread(fid,1,'float32')); % PS rhscaleqdef = 0		+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_rotation', fread(fid,1,'int16')); % RX 0 90 180 270 deg		0 - 3 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_transpose', fread(fid,1,'int16')); % RX 0, 1 n / y transpose 	0 - 1%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_kissoff_views', fread(fid,1,'int16')); % rhblank zero image views	0 - 512 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_slblank', fread(fid,1,'int16')); % rhslblankslice blank 3d	0 - 128 % 
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_gradcoil', fread(fid,1,'int16')); % RX 0=off 1=Schnk 2=Rmr	0 - 2 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ddaover', fread(fid,1,'int16')); % rhddaover unused %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_sarr', fread(fid,1,'int16')); % SARR bit mask 		15 bits %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_fd_tr', fread(fid,1,'int16')); % SARR feeder timing info %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_fd_te', fread(fid,1,'int16')); % SARR feeder timing info %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_fd_ctrl', fread(fid,1,'int16')); % SARR control of feeder %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_algor_num', fread(fid,1,'int16')); % SARR df decimation ratio %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_fd_df_dec', fread(fid,1,'int16')); % SARR which feeder algor %
% NOTE: this is different from the rdbm.h style %
% rdb_hdr_dab is a 4 element array of structures of type RDB_MULTI_RCV_TYP
% RDB_MULTI_RCV_TYP is a struct of two elements
% typedef struct
% {
% 	short start_rcv;
% 	short stop_rcv;
% } RDB_MULTI_RCV_TYPE;
% instead of this I'm using 8 shorts named as expected
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_dab1_start_rcv', fread(fid,1,'int16')); % rhdab0s rhdab0e st, stp rcv 	0 - 15 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_dab1_stop_rcv', fread(fid,1,'int16')); % rhdab0s rhdab0e st, stp rcv 	0 - 15 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_dab2_start_rcv', fread(fid,1,'int16')); % rhdab0s rhdab0e st, stp rcv 	0 - 15 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_dab2_stop_rcv', fread(fid,1,'int16')); % rhdab0s rhdab0e st, stp rcv 	0 - 15 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_dab3_start_rcv', fread(fid,1,'int16')); % rhdab0s rhdab0e st, stp rcv 	0 - 15 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_dab3_stop_rcv', fread(fid,1,'int16')); % rhdab0s rhdab0e st, stp rcv 	0 - 15 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_dab4_start_rcv', fread(fid,1,'int16')); % rhdab0s rhdab0e st, stp rcv 	0 - 15 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_dab4_stop_rcv', fread(fid,1,'int16')); % rhdab0s rhdab0e st, stp rcv 	0 - 15 %
% end change
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user0', fread(fid,1,'float32')); % rhuser0 			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user1', fread(fid,1,'float32')); % rhuser1 			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user2', fread(fid,1,'float32')); % rhuser2 			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user3', fread(fid,1,'float32')); % rhuser3 			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user4', fread(fid,1,'float32')); % rhuser4 			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user5', fread(fid,1,'float32')); % rhuser5 			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user6', fread(fid,1,'float32')); % rhuser6 			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user7', fread(fid,1,'float32')); % rhuser7 			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user8', fread(fid,1,'float32')); % rhuser8 			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user9', fread(fid,1,'float32')); % rhuser9 			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user10', fread(fid,1,'float32')); % rhuser10			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user11', fread(fid,1,'float32')); % rhuser11			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user12', fread(fid,1,'float32')); % rhuser12			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user13', fread(fid,1,'float32')); % rhuser13			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user14', fread(fid,1,'float32')); % rhuser14			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user15', fread(fid,1,'float32')); % rhuser15			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user16', fread(fid,1,'float32')); % rhuser16			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user17', fread(fid,1,'float32')); % rhuser17			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user18', fread(fid,1,'float32')); % rhuser18			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user19', fread(fid,1,'float32')); % rhuser19			+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_type', fread(fid,1,'int32'));	% rhvtypebit mask		31 bits %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefxa', fread(fid,1,'float32'));	 % RX x flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefxb', fread(fid,1,'float32'));	 % RX x flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefxc', fread(fid,1,'float32'));	 % RX x flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefxd', fread(fid,1,'float32'));	 % RX x flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefya', fread(fid,1,'float32'));	 % RX y flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefyb', fread(fid,1,'float32'));	 % RX y flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefyc', fread(fid,1,'float32'));	 % RX y flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefyd', fread(fid,1,'float32'));	 % RX y flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefza', fread(fid,1,'float32'));	 % RX z flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefzb', fread(fid,1,'float32'));	 % RX z flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefzc', fread(fid,1,'float32'));	 % RX z flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_coefzd', fread(fid,1,'float32'));	 % RX z flow direction control	0 - 4 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vm_coef1', fread(fid,1,'float32'));	 % RX weight for mag image 1	0 - 1 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vm_coef2', fread(fid,1,'float32'));	 % RX weight for mag image 2	0 - 1 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vm_coef3', fread(fid,1,'float32'));	 % RX weight for mag image 3	0 - 1 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vm_coef4', fread(fid,1,'float32'));	 % RX weight for mag image 4	0 - 1 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_v_venc', fread(fid,1,'float32')); % RX vel encodeing cm / sec	0.001 - 5000 %
rdb_hdr = setfield(rdb_hdr, 'spectral_width', fread(fid,1,'float32')); % specwidthfilter width kHz	500 - 3355432 %
rdb_hdr = setfield(rdb_hdr, 'csi_dims', fread(fid,1,'int16')); % spectro %
rdb_hdr = setfield(rdb_hdr, 'xcsi', fread(fid,1,'int16')); % rhspecrescsix		2 - 64 %
rdb_hdr = setfield(rdb_hdr, 'ycsi', fread(fid,1,'int16')); % rhspecrescsiy		2 - 64 %
rdb_hdr = setfield(rdb_hdr, 'zcsi', fread(fid,1,'int16')); % spectro %
rdb_hdr = setfield(rdb_hdr, 'roilenx', fread(fid,1,'float32')); % RX x csi volume dimension %
rdb_hdr = setfield(rdb_hdr, 'roileny', fread(fid,1,'float32')); % RX y csi volume dimension %
rdb_hdr = setfield(rdb_hdr, 'roilenz', fread(fid,1,'float32')); % RX z csi volume dimension %
rdb_hdr = setfield(rdb_hdr, 'roilocx', fread(fid,1,'float32')); % RX x csi volume center %
rdb_hdr = setfield(rdb_hdr, 'roilocy', fread(fid,1,'float32')); % RX y csi volume center %
rdb_hdr = setfield(rdb_hdr, 'roilocz', fread(fid,1,'float32')); % RX z csi volume center %
rdb_hdr = setfield(rdb_hdr, 'numdwell', fread(fid,1,'float32')); % specdwells			0 - 3.40282e38%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_command', fread(fid,1,'int32'));	 % PS internal use only	%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_mps_r1', fread(fid,1,'int32'));	 % PS MPS R1 setting		1 - 7 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_mps_r2', fread(fid,1,'int32'));	 % PS MPS R2 setting		1 - 30 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_mps_tg', fread(fid,1,'int32'));	 % PS MPS Transmit gain setting	0 - 200%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_mps_freq', fread(fid,1,'int32')); % PS MPS Center frequency hz	+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_aps_r1', fread(fid,1,'int32'));	 % PS APS R1 setting		1 - 7 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_aps_r2', fread(fid,1,'int32'));	 % PS APS R2 setting		1 - 30 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_aps_tg', fread(fid,1,'int32'));	 % PS APS Transmit gain setting	0 - 200%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_aps_freq', fread(fid,1,'int32')); % PS APS Center frequency hz	+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_scalei', fread(fid,1,'float32')); % PS rational scaling 		+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_scaleq', fread(fid,1,'float32')); % PS unused %			
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_snr_warning', fread(fid,1,'int32')); % PS noise test 0=16 1=32 bits	0, 1 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_aps_or_mps', fread(fid,1,'int32')); % PS prescan order logic	0 - 5 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_mps_bitmap', fread(fid,1,'int32')); % PS bit mask			4 bits%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_powerspec', fread(fid,256,'char')); % PS %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_filler1', fread(fid,1,'int32'));	 % PS filler %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ps_filler2', fread(fid,1,'int32'));	 % PS filler %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_rec_noise_mean', fread(fid,16,'float32')); % PS mean noise each receiver +-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_rec_noise_std', fread(fid,16,'float32')); % PS noise calc for muti rec	+-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'halfecho', fread(fid,1,'int16')); % spectro full, half echo 0, 1 %
% 858 bytes %
% New fields 02-19-92 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_im_size_y', fread(fid,1,'int16')); % rh???? 			0 - 512 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_data_collect_type1', fread(fid,1,'int32')); % rh???? bit mask		31 bits %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_freq_scale', fread(fid,1,'float32')); % rh???? freq k-space step +-3.40282e38 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_phase_scale', fread(fid,1,'float32')); % rh???? freq k-space step +-3.40282e38 %
% 14 bytes %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ovl', fread(fid,1,'int16')); % rhovl - overlaps for MOTSA % 
% Phase Correction Control Param. %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pclin', fread(fid,1,'int16')); % Linear Corr. 0:off, 1:linear, 2:polynomial %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pclinnpts', fread(fid,1,'int16')); % fit number of points %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pclinorder', fread(fid,1,'int16')); % fit order %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pclinavg', fread(fid,1,'int16')); % linear phase corr avg 0:off, 1:on %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pccon', fread(fid,1,'int16')); % Const Corr. 0:off, 1:Ky spec., 2:polyfit(2/ilv), 3:polyfit(1/ilv) %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcconnpts', fread(fid,1,'int16')); % fit number of points %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcconorder', fread(fid,1,'int16')); % fit order %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcextcorr', fread(fid,1,'int16')); % external correction file 0:don't use, 1: use %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcgraph', fread(fid,1,'int16')); % Phase Correction coef. image 0:off, 1:linear & constant %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcileave', fread(fid,1,'int16')); % Interleaves to use for correction: 0=all, 1=only first %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_hdbestky', fread(fid,1,'int16')); % bestky view for fractional Ky scan %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcctrl', fread(fid,1,'int16')); % phase correction research control %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcthrespts', fread(fid,1,'int16')); % 2..512 adjacent points %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcdiscbeg', fread(fid,1,'int16')); % 0..512 beginning point to discard %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcdiscmid', fread(fid,1,'int16')); % 0..512 middle point to discard %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcdiscend', fread(fid,1,'int16')); % 0..512 ending point to discard %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcthrespct', fread(fid,1,'int16')); % Threshold percentage %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcspacial', fread(fid,1,'int16')); % Spacial best ref scan index 0..512 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pctemporal', fread(fid,1,'int16')); % Temporal best ref scan index 0..512 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcspare', fread(fid,1,'int16')); % spare for phase correction %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ileaves', fread(fid,1,'int16')); % Nunmber of interleaves %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_kydir', fread(fid,1,'int16')); % Ky traversal dircetion 0: top-down, 1:center out %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_alt', fread(fid,1,'int16')); % Alt read sign 0=no, 1=odd/even, 2=pairs %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_reps', fread(fid,1,'int16')); % Number of scan repetitions %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ref', fread(fid,1,'int16')); % Ref Scan 0: off 1: on %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcconnorm', fread(fid,1,'float32')); % Constant S term normalization factor %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcconfitwt', fread(fid,1,'float32')); % Constant polyfit weighting factor %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pclinnorm', fread(fid,1,'float32')); % Linear S term normalization factor %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pclinfitwt', fread(fid,1,'float32')); % Linear polyfit weighting factor %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcbestky', fread(fid,1,'float32')); % Best Ky location %
% VRG Filter param %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vrgf', fread(fid,1,'int32')); % control word for VRG filter %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vrgfxres', fread(fid,1,'int32')); % control word for VRGF final x resolution %
% Bandpass AsymmetryCorrection Param. %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_bp_corr', fread(fid,1,'int32')); % control word for bandpass asymmetry %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_recv_freq_s', fread(fid,1,'float32')); % starting frequency (+62.5) %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_recv_freq_e', fread(fid,1,'float32')); % ending frequency (-62.5) %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_hniter', fread(fid,1,'int32')); % Selects the number of
% iterations used in homodyne processing %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_fast_rec', fread(fid,1,'int32')); % Added for homodyne II, tells if
% the fast receiver is being used
% and the lpf setting of the fast
				% receiver, 0: fast receiver off,
				% 1 - 5: lpf settings %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_refframes', fread(fid,1,'int32')); % total # of frames for ref scan %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_refframep', fread(fid,1,'int32')); % # of frames per pass for a ref scan %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_scnframe', fread(fid,1,'int32')); % total # of frames for a entire scan %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pasframe', fread(fid,1,'int32')); % # of frames per pass %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user_usage_tag', fread(fid,1,'uint32')); % for spectro %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user_fill_mapMSW', fread(fid,1,'uint32')); % for spectro %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user_fill_mapLSW', fread(fid,1,'uint32')); % for Spectro %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user20', fread(fid,1,'float32')); % all following usercv are for spectro %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user21', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user22', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user23', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user24', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user25', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user26', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user27', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user28', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user29', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user30', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user31', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user32', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user33', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user34', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user35', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user36', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user37', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user38', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user39', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user40', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user41', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user42', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user43', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user44', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user45', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user46', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user47', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_user48', fread(fid,1,'float32'));
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcfitorig', fread(fid,1,'int16')); % Adjust view indexes if set so bestky view = 0 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcshotfirst', fread(fid,1,'int16')); % First view within an echo group used for fit%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcshotlast', fread(fid,1,'int16')); % Last view within an echo group used for fit %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcmultegrp', fread(fid,1,'int16')); % If = 1, force pts from other egrps to be used %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pclinfix', fread(fid,1,'int16')); % If = 2, force slope to be set to pclinslope %
 % If = 1, neg readout slope = pos readout slope %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcconfix', fread(fid,1,'int16')); % If = 2, force slope to be set to pcconslope %
 % If = 1, neg readout slope = pos readout slope %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pclinslope', fread(fid,1,'float32')); % Value to set lin slope to if forced %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pcconslope', fread(fid,1,'float32')); % Value to set con slope to if forced %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_pccoil', fread(fid,1,'int16')); % If 1,2,3,4, use that coil's results for all %
% Variable View Sharing %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vvsmode', fread(fid,1,'int16')); % Variable view sharing mode %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vvsaimgs', fread(fid,1,'int16')); % number of original images %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vvstr', fread(fid,1,'int16')); % TR in microseconds %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vvsgender', fread(fid,1,'int16')); % gender: male or female %
% 3D Slice ZIP %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_zip_factor', fread(fid,1,'int16')); % Slice ZIP factor: 0=OFF, 2, or 4 %
% Maxwell Term Correction Coefficients %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef1a', fread(fid,1,'float32')); % Coefficient A for flow image 1 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef1b', fread(fid,1,'float32')); % Coefficient B for flow image 1 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef1c', fread(fid,1,'float32')); % Coefficient C for flow image 1 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef1d', fread(fid,1,'float32')); % Coefficient D for flow image 1 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef2a', fread(fid,1,'float32')); % Coefficient A for flow image 2 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef2b', fread(fid,1,'float32')); % Coefficient B for flow image 2 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef2c', fread(fid,1,'float32')); % Coefficient C for flow image 2 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef2d', fread(fid,1,'float32')); % Coefficient D for flow image 2 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef3a', fread(fid,1,'float32')); % Coefficient A for flow image 3 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef3b', fread(fid,1,'float32')); % Coefficient B for flow image 3 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef3c', fread(fid,1,'float32')); % Coefficient C for flow image 3 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_maxcoef3d', fread(fid,1,'float32')); % Coefficient D for flow image 3 %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_ut_ctrl', fread(fid,1,'int32')); % System utility control variable %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_dp_type', fread(fid,1,'int16')); % EPI II diffusion control cv %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_arw', fread(fid,1,'int16')); % Arrhythmia rejection window(percentage:1-100)%
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_vps', fread(fid,1,'int16'));	% View Per Segment for FastCine % 
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_mcReconEnable', fread(fid,1,'int16')); % N-Coil recon map %
rdb_hdr = setfield(rdb_hdr, 'rdb_hdr_excess', fread(fid,420,'int16')); % free space for later expansion %

return;

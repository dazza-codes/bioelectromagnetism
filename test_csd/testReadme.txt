
Hi,

 If you create a model file with the dipole simulator and save the data, 
you can plot the data drom an average or a csd/laplacian reference (among 
others). I used 'eeg_interp_sph_spline_test_150.m' to calculate the csd
at each time point, I then compared the waveform to the dipole simulator results.

You will have to erase the top line from the diple simulator generated .avr file.
I have included an electrodes file (dip_sim_elecs.elp) with spherical coordinates. This file 
was modified from the one output by the dipole simulator to include the additional information
the eeg_toolbox code looks for. NOTE: theta and phi are inverted from your sph2cart function... 
I have also included the model I used to generate the data, but I added some noise independently.

Please let me know if you have any questions,

Signe


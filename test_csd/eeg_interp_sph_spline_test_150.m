% eeg_interp_sph_sline_test_150 - 
%	script to test eeg_interp_sph_spline_scd at 150 time points
%

clear

p = eeg_toolbox_defaults;

testPath = fullfile(eeg_toolbox_path,'test_csd','');
cd(testPath)

p.elec.path = testPath;
p.elec.file = 'dip_sim_elecs.elp';
p.elec.type = 'Spherical1';	%SIGNE: modified from 'scan3ddasc'
p.elec.plot = 0;

p.volt.path = testPath;
p.volt.file = 's1_cue1_noisy_mod_new.avr'; % doesn't exist?
p.volt.file = 's1_cue1_noisy_mod.avr';
p.volt.type = 'ascii';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load electrode co-ordinates

p = elec_open(p);
p.elec.plot = 1;

x = p.elec.data.x;
y = p.elec.data.y;
z = p.elec.data.z;

p.volt.sampleMsec = 2; %1;
p.volt.sampleTime = 1 * p.volt.sampleMsec;
p.volt.samplePoint = 1;
p.volt.epochStart = 0;
p.volt.epochEnd = 300;
p.volt.timeArray = [0:2:300];

p = eeg_open(p);

for i = 1:150,
  
  p.volt.samplePoint = i;
  
  V = p.volt.data(i,:)';
  
  % scd interpolation
  
  FV = eeg_interp_sph_spline_scd_sig(V,[x y z]);
  
  data(i,:) = FV.Cdata;
end

% data(i, j) contains the value at the ith time point of the jth electrode

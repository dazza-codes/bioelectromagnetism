function [F,time_labels,Rcoils,Ocoils,stim_num,coilWeights,nc_montage,bad_channels_list,baseline_variance,EEGwf, ...
   EEGpickupLoc,EEGreferenceLoc,EEGpickupRadius,EEGreferenceRadius,success] = load_AvgNetCDFv5( ...
   ncfnames,stim_num)

% load_AvgNetCDFv5 :  reads a set of netMEG files, and returns the good data for the specified stimulus
%       for Matlab 5.0 and larger
% 
% useage [F,time_labels,Rcoils,Ocoils,stim_num,coilWeights,nc_montage,bad_channels_list,baseline_variance,EEGwf, ...
%  EEGpickupLoc,EEGreferenceLoc,EEGpickupRadius,EEGreferenceRadius,success] = load_AvgNetCDF( ...
%  ncnames,stim_num)
%
% Input:
%    ncfnames ... array of names of all netCDF files in the set
%    stim_num ...  stim to read from each file in the set (scalar)
%
% Output:
% F       : Data matrix of MEG measurements - probes x latencies
% time_labels : Row vector of latencies 
% Rcoils  : MEG sensor Location matrix - probes X (3 [x y z ])*ncoils 
%           contains location info for each coil
% Ocoils  : MEG sensor Orientation matrix - probes X (3 [x' y' z'])*ncoils 
%           orientation info for each coil
% stim_num: Stimulus number for the experiment
% coilWeights : vector of the weight to be applied to the values calculated at each MEG sensor coil. (ncoils)
%  (assumes all sensors have the same weighting. This could be changed -- the netMEG file saves the weights
%   separately for each sensor)
% nc_montage : string, name of sensor system
% bad_channels_list : string, list of channels whose data was not used
% baseline_variance : column vector of noise for the selected stim (nsensors,1)
% EEGwf : Data matrix of EEG measurements - sensors x latencies
% EEGpickupLoc : location of EEG pickup sensor, 3 X EEGsensors 
% EEGreferenceLoc : location of EEG reference sensor,  3 X EEGsensors 
% EEGpickupRadius : Radius of EEG pickup sensor,  EEGsensors
% EEGreferenceRadius : Radius of EEG reference sensor,  EEGsensors
% success : if the file couldn't be read for any reason, success = 0
%
% Elaine Best, Los Alamos National Laboratory, 10/95
% Altered 1/6/98 to read version 1.4 files, and warn of old ones
% Altered 7/98 to also work with new version of netCDF library for MATLAB 5.0
% 

success = 1;

% loop to get the data from each file
[nfiles,d] = size(ncfnames);
Rcoils = [];
Ocoils = [];

  % loop over netMEG files
  for ifile = 1:nfiles
  
    % open the file in read only mode
    filename = ncfnames(ifile,:);
    nc = netcdf(filename,'nowrite');

    % get a list of the variable names
    variables = var(nc);
    var_names = ncnames(variables);
  
    % find out which sensor array (montage) and netCDF file version we are using
      nc_montage = nc.MontageName(:);
      bad_channels_list = nc.BadChannelsDeleted(:);
      nc_version = nc.netCDFfileVersion(:);
  
       samp_int = nc{'SamplingInterval'}(:);
       prestim = nc{'LengthOfPrestim'}(stim_num);
       npts = nc{'numSamples'}(stim_num);
  
       time_labels = (0:npts-1)*samp_int - prestim;
  
     if (strcmp(nc_version(7:9),'1.1') | strcmp(nc_version(7:9),'1.2') )
       disp(' ')
       disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ')
       disp('ERROR: ')
       disp('Your netMEG file is an obsolete version.')
       disp('You can read it into MEGAN to update it.')
       disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ')
       disp(' ')
       success = 0;
       return
     end
     if (strcmp(nc_version(7:9),'1.3') )
       disp('Your netMEG file is a bad version.')
       disp('Only a few of these files exist.  There is no')
       disp('way to update them. Please discard the file.')
       success = 0;
       return
     end
  
    % get the channel type array 
    ind = find(strcmp(var_names,'ChannelTypes'));
    if length(ind) ~= 0 
      ChannelTypes = nc{'ChannelTypes'}(:);
    else 
      disp('Returning: netMEG Error reading  reading ChannelTypes')
      success = 0;
      return;
    end
    [nchans,dummy] = size(ChannelTypes);
  
    % get the channel status array if it exists, otherwise make one that indicates
    % that all of the data is good
    ind = find(strcmp(var_names,'ChannelStatus'));
    if length(ind) ~= 0 
      ChannelStatus = nc{'ChannelStatus'}(:);
    else 
      ChannelStatus = ones(nchans(1),1) 
    end
  
    % See if we have any EEG sensor location data, read it if we do.
    % (It is possible to have EEG channel data without EEG sensor locations.
    %  However, we can't localize without the sensor locations.)

    % EEG pickup  and reference locations and radii
    ind = find(strcmp(var_names,'EEGpickupLocation'));
    if length(ind) ~= 0 
      EEGpickupLoctmp = nc{'EEGpickupLocation'}(:);
      [neeg_sensors,dummy] = size(EEGpickupLoctmp);
    else 
      neeg_sensors = 0;
    end

    if neeg_sensors > 0 
        ind = find(strcmp(var_names,'EEGreferenceLocation'));
        if length(ind) ~= 0 
          EEGreferenceLoctmp = nc{'EEGreferenceLocation'}(:);
        else 
          disp('Returning: netMEG Error reading  EEGreferenceLocation')
          success = 0;
          return;
        end

        ind = find(strcmp(var_names,'EEGpickupRadius'));
        if length(ind) ~= 0 
          EEGpickupRadiustmp = nc{'EEGpickupRadius'}(:);
        else 
          disp('Returning: netMEG Error reading  EEGpickupRadius')
          success = 0;
          return;
        end

        ind = find(strcmp(var_names,'EEGreferenceRadius'));
        if length(ind) ~= 0 
          EEGreferenceRadiustmp = nc{'EEGreferenceRadius'}(:);
        else 
          disp('Returning: netMEG Error reading  EEGreferenceRadius')
          success = 0;
          return;
        end
    end

    % figure out which channels are MEG and which (if any) are EEG, and
    % select data only for channels with a status of 1.
      
    good_MEG_indices = find(ChannelTypes(:,1)=='M' & ChannelTypes(:,2) == 'E' & ...
             ChannelTypes(:,3) == 'G' & ChannelStatus == 1);
    ngood_MEG = size(good_MEG_indices,1);

    if (neeg_sensors > 0)
        good_EEG_indices = find(ChannelTypes(:,1)=='E' & ChannelTypes(:,2) == 'E' & ...
             ChannelTypes(:,3) == 'G' & ChannelStatus == 1);
        EEG_indices = find(ChannelTypes(:,1)=='E' & ChannelTypes(:,2) == 'E' & ...
             ChannelTypes(:,3) == 'G' );
        EEG_sensor_status = ChannelStatus(EEG_indices);
        good_EEG_sensors = find(EEG_sensor_status == 1);
         
        ngood_EEG = size(good_EEG_indices,1);
        if ifile == 1
           EEGpickupLoc = EEGpickupLoctmp(good_EEG_sensors,:)';
           EEGreferenceLoc = EEGreferenceLoctmp(good_EEG_sensors,:)';
           EEGpickupRadius = EEGpickupRadiustmp(good_EEG_sensors);
           EEGreferenceRadius = EEGreferenceRadiustmp(good_EEG_sensors);
        else
           EEGpickupLoc = [EEGpickupLoc,EEGpickupLoctmp(good_EEG_sensors,:)'];
           EEGreferenceLoc = [EEGreferenceLoc,EEGreferenceLoctmp(good_EEG_sensors,:)'];
           EEGpickupRadius = [EEGpickupRadius;EEGpickupRadiustmp(good_EEG_sensors)];
           EEGreferenceRadius = [EEGreferenceRadius;EEGreferenceRadiustmp(good_EEG_sensors)];
        end
    else
        ngood_EEG = 0;
        EEGpickupLoc = 0;
        EEGreferenceLoc = 0;
        EEGpickupRadius = 0;
        EEGreferenceRadius = 0;
    end
  
    if (ngood_MEG > 0 )
      % read the coil location and orientation data
      % and concatenate it with the previously read data, if any
    
      ind = find(strcmp(var_names,'SensorElementsOrient'));
      if length(ind) ~= 0 
        SensorElementsOrient = nc{'SensorElementsOrient'}(:);
      else 
        disp('Returning: netMEG Error reading  SensorElementsOrient')
        success = 0;
        return;
      end
    
      ind = find(strcmp(var_names,'SensorElementsLoc'));
      if length(ind) ~= 0 
        SensorElementsLoc = nc{'SensorElementsLoc'}(:);
      else 
        disp('Returning: netMEG Error reading  SensorElementsLoc')
        success = 0;
        return;
      end
  
      [dum1,max_sens_elem,dum2] = size(SensorElementsLoc);
    
      % loop over the coils, and load the data for each coil into
      % the Rcoils and Ocoils arrays
      % (Matlab 4 could not handle arrays larger than 2 dimensions, so we store
      % the data in a 2 D array.)
      for igood = 1:ngood_MEG
        isn = good_MEG_indices(igood);
        coil_loc = reshape(SensorElementsLoc(isn,:,:),max_sens_elem,3);
        sens_orient = reshape(SensorElementsOrient(isn,:,:),max_sens_elem,3);
      
        % load the coil locations and orientations into the Rcoils and Ocoils arrays
        if ifile == 1 & igood == 1
          Rcoils = [coil_loc(1,:)];
          Ocoils = [sens_orient(1,:)];
          for j = 2: max_sens_elem
            Rcoils = [Rcoils,coil_loc(j,:)];
            Ocoils = [Ocoils,sens_orient(j,:)];
          end
        else
          Rcoils_row = [coil_loc(1,:)];
          Ocoils_row = [sens_orient(1,:)];
          for j = 2: max_sens_elem
            Rcoils_row = [Rcoils_row,coil_loc(j,:)];
            Ocoils_row = [Ocoils_row,sens_orient(j,:)];
          end
          Rcoils = [Rcoils;Rcoils_row];
          Ocoils = [Ocoils;Ocoils_row];
        end
      end %for

      % get the coil weights
      % we store the weights for each sensor, but for now are assuming that they
      % are the same for each sensor, and are just using a vector.
       ind = find(strcmp(var_names,'CoilWeight'));
       if length(ind) ~= 0 
         coilWeight = nc{'CoilWeight'}(:);
       else 
         disp('Returning: netMEG Error reading  CoilWeight')
         success = 0;
         return;
       end
       coilWeights = coilWeight(1,:)';
  
      % read the MEG noise
       ind = find(strcmp(var_names,'BaselineVariance'));
       if length(ind) ~= 0 
         baseline_variance_tmp = nc{'BaselineVariance'}(:);
       else 
         disp('Returning: netMEG Error reading  BaselineVariance')
         success = 0;
         return;
       end

    end % ngood_MEG
  
    % read the amplitudes
    % get the full array of data for this stim
    % from the file, and then extract the data we want
    ind = find(strcmp(var_names,'Waveforms'));
    if length(ind) ~= 0 
      Ftmp = nc{'Waveforms'}(stim_num,1:npts,:);
      [np,nch,dum] = size(Ftmp);
      Ftmp = Ftmp';
    else 
      disp('Returning: netMEG Error reading  Waveforms')
      success = 0;
      return;
    end
  
    % MEG waveform data
      if (ngood_MEG > 0)
        if ifile == 1
           F = Ftmp(good_MEG_indices,:);
           baseline_variance = baseline_variance_tmp(stim_num,good_MEG_indices)';
        else
           F = [F ; Ftmp(good_MEG_indices,:)];
           baseline_variance = [baseline_variance;baseline_variance_tmp(stim_num,good_MEG_indices)'];
        end
      end

    % EEG waveform data
      if (ngood_EEG > 0)
        % EEG Waveforms
        if ifile == 1
           EEGwf = Ftmp(good_EEG_indices,:);
        else
           EEGwf = [EEGwf ; Ftmp(good_EEG_indices,:)];
        end
      else
        EEGwf = 0;
      end
  
    % close the file
    nc = close(nc);
  
  end % for loop over files
  time_labels = time_labels';

end % code for version 5

return

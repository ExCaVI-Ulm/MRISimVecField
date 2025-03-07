function [settings] = calcNyquist2D(settings)
%calcNyquist2D depending on the sequence type chosen, Nyquists theorem is
%applied and used to calculate timings, durations and sampling frequencies
%as for conventional MRI (downsampled!)
%   Input:  -settings struct
%  
%   Output:
%           -settings struct


settings.trajectory.dK = 1/settings.reco.FOV;
settings.trajectory.pixel_width = settings.reco.FOV / settings.reco.matrixsize_reco;
settings.trajectory.k_FOV = 1 / settings.trajectory.pixel_width;

%Nyquist sampling
settings.trajectory.Tread = settings.trajectory.k_FOV *2*pi / ( settings.general.gamma * settings.trajectory.GRstrength_est );   %readout length

settings.trajectory.TOversamp = 1;  %Oversampling factor
settings.trajectory.BW = settings.trajectory.TOversamp*settings.general.gamma * settings.trajectory.GRstrength_est * settings.reco.FOV / (2*pi);% sampling Bandwidth

switch settings.trajectory.type
    case 'ImportFiles'
        settings.CST.path_txtFiles = [settings.general.curr_path, '/FieldData/'];     %path of folder where .txt files are saved
        
        %determine number of txt-files in folder to get number of rotations/"phase-encoding"-steps
        tmp_folderinfo = dir([settings.CST.path_txtFiles '/*.txt']);        
        settings.trajectory.N_PhaseEnc = length(tmp_folderinfo);
        
        settings.trajectory.BW =settings.trajectory.BW; %either take values calculated from Nyquist Sampling / standard encoding or set new values
        settings.trajectory.Tread = 0.5*settings.trajectory.Tread;  %above: "cartesian" calculation but it's rather UTE-like -> half readout length
        settings.trajectory.Nsamples = ceil(settings.trajectory.Tread * settings.trajectory.BW); %number of samples
        
        clearvars tmp_folderinfo
    case 'Rect1D'
        settings.trajectory.N_PhaseEnc = 1;
        settings.trajectory.BW =settings.trajectory.BW; %either take values calculated from Nyquist Sampling / standard encoding or set new values
        settings.trajectory.Tread = settings.trajectory.Tread;  %above: "cartesian" calculation but it's rather UTE-like -> half readout length
        settings.trajectory.Nsamples = (settings.trajectory.Tread * settings.trajectory.BW); %number of samples
    case 'Cart2D'
        settings.trajectory.N_PhaseEnc = round(settings.trajectory.k_FOV / settings.trajectory.dK); %in case of cartesian sampling
        settings.trajectory.BW =settings.trajectory.BW;
        settings.trajectory.Tread = settings.trajectory.Tread;
        settings.trajectory.Nsamples = round(settings.trajectory.Tread * settings.trajectory.BW);
    otherwise
        warning('Unrecognized sampling pattern')
end

end
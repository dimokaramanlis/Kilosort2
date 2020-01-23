function [ops] = getKs2OptionsMEA(metadata)
%GETKSOPTIONSMEA Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
% options for reading/saving files
ops.root        = metadata.root;
ops.fbinary     = metadata.binpath;
ops.fproc       = metadata.whpath;
%=========================================================================
ops.trange = [0 Inf]; % time range to sort
%=========================================================================
% array-specific options
switch metadata.meatype
    case '252MEA10030'
        meaChannelMap([16 16], 100,  fullfile(ops.root, 'ks_sorted'), 1); 
        ops.nfilt_factor        = 4;
        ops.NchanTOT            = 252;
    case '252MEA20030'    
        meaChannelMap([16 16], 200,  fullfile(ops.root, 'ks_sorted'), 1); 
        ops.nfilt_factor        = 4;
        ops.NchanTOT            = 252;
    case '60MEA10030'
        meaChannelMap([8 8], 100,  fullfile(ops.root, 'ks_sorted'), 1);
        ops.nfilt_factor        = 4;
        ops.NchanTOT            = 60;
    case '60MEA20030'
        meaChannelMap([8 8], 200,  fullfile(ops.root, 'ks_sorted'), 1);
        ops.nfilt_factor        = 4;
        ops.NchanTOT            = 60;
    case '60MEA10010'    
        meaChannelMap([8 8], 100,  fullfile(ops.root, 'ks_sorted'), 1); 
        ops.nfilt_factor        = 4;        
        ops.NchanTOT            = 60;
end
ops.chanMap             = fullfile(ops.root, 'ks_sorted','chanMap.mat'); % make this file using createChannelMapFile.m		
%==========================================================================
%extra ops
ops.min_NchanNear = 32;
ops.min_Nnearest  = 32;
%==========================================================================
% sample rate
ops.fs                  = metadata.bininfo.fs; %sampling frequency		

% frequency for high pass filtering (150)
ops.fshigh = 150;   

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0.1; 

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [8 4];  

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.lam = 60;  

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.9; 

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 1/50; 

% number of samples to average over (annealed from first to second value) 
ops.momentum = [40 800]; 

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30; 

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8; 

% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
ops.reorder         = 1;       % whether to reorder batches for drift correction. 
ops.nskip           = 20;  % how many batches to skip for determining spike PCs

ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 64*round(1024*ops.fs/1e4) + ops.ntbuff;% this is the batch size (try decreasing if out of memory) 		

ops.whiteningRange      = 32; % number of channels to use for whitening each channel
ops.nSkipCov            = 10; % compute whitening matrix from every N-th batch
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.nPCs                = 3; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available
%==========================================================================	
%hidden options
ops.nt0             = floor(round(21*ops.fs/1e4)/2)*2+1; %spike template time bins
ops.nt0min          = floor(round(7*ops.fs/1e4)/2)*2; %spike template time bins
ops.loc_range       = [round(2*ops.fs/1e4) 1];  % ranges to detect peaks; plus/minus in time and channel ([5 4])		
ops.long_range      = [round(12*ops.fs/1e4) 1]; % ranges to detect isolated peaks ([30 6])		
%==========================================================================
if not(isempty(metadata.exptypes))
    ops = updateOpsMEA(metadata, ops);
end
end


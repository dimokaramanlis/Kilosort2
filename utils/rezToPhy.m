
function rezToPhy(rez, savePath)
% pull out results from kilosort's rez to either return to workspace or to
% save in the appropriate format for the phy GUI to run on. If you provide
% a savePath it should be a folder, and you will need to have npy-matlab
% available (https://github.com/kwikteam/npy-matlab)


% spikeTimes will be in samples, not seconds
rez.W = gather(single(rez.Wphy));
rez.U = gather(single(rez.U));
rez.mu = gather(single(rez.mu));

if size(rez.st3,2)>4
    rez.st3 = rez.st3(:,1:4);
end

% not needed, happened already
% [~, isort]   = sort(rez.st3(:,1), 'ascend'); 
% rez.st3      = rez.st3(isort, :);
% rez.cProj    = rez.cProj(isort, :);
% rez.cProjPC  = rez.cProjPC(isort, :, :);


fs = dir(fullfile(savePath, '*.npy'));
for i = 1:length(fs)
   delete(fullfile(savePath, fs(i).name));
end
if exist(fullfile(savePath, '.phy'), 'dir')
    rmdir(fullfile(savePath, '.phy'), 's');
end

spikeTimes = uint64(rez.st3(:,1));
% [spikeTimes, ii] = sort(spikeTimes);
spikeTemplates = uint32(rez.st3(:,2));
if size(rez.st3,2)>4
    spikeClusters = uint32(1+rez.st3(:,5));
end
amplitudes = rez.st3(:,3);

Nchan = rez.ops.Nchan;

xcoords     = rez.xcoords(:);
ycoords     = rez.ycoords(:);
chanMap     = rez.ops.chanMap(:);
chanMap0ind = chanMap - 1;

nt0 = size(rez.W,1);
U = rez.U;
W = rez.W;

Nfilt = size(W,2);

templates = zeros(Nchan, nt0, Nfilt, 'single');
for iNN = 1:size(templates,3)
   templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))';
end
templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
templatesInds = repmat([0:size(templates,3)-1], size(templates,1), 1); % we include all channels so this is trivial

templateFeatureInds = uint32(rez.iNeigh);
pcFeatureInds = uint32(rez.iNeighPC);

whiteningMatrix = rez.Wrot/rez.ops.scaleproc;
whiteningMatrixInv = whiteningMatrix^-1;

% here we compute the amplitude of every template...

% unwhiten all the templates
tempsUnW = zeros(size(templates));
for t = 1:size(templates,1)
    tempsUnW(t,:,:) = squeeze(templates(t,:,:))*whiteningMatrixInv;
end

% The amplitude on each channel is the positive peak minus the negative
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

% The template amplitude is the amplitude of its largest channel
tempAmpsUnscaled = max(tempChanAmps,[],2);

% assign all spikes the amplitude of their template multiplied by their
% scaling amplitudes
spikeAmps = tempAmpsUnscaled(spikeTemplates).*amplitudes;

% take the average of all spike amps to get actual template amps (since
% tempScalingAmps are equal mean for all templates)
ta = clusterAverage(spikeTemplates, spikeAmps); clear spikeAmps;
tids = unique(spikeTemplates);
tempAmps(tids) = ta; % because ta only has entries for templates that had at least one spike
gain = getOr(rez.ops, 'gain', 1);
tempAmps = gain*tempAmps'; % for consistency, make first dimension template number

if ~isempty(savePath)

    writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy')); clear spikeTimes;
    writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_templates.npy')); % -1 for zero indexing
    if size(rez.st3,2)>4
        writeNPY(uint32(spikeClusters-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    else
        writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    end
    writeNPY(amplitudes, fullfile(savePath, 'amplitudes.npy')); clear amplitudes;
    writeNPY(templates, fullfile(savePath, 'templates.npy')); clear templates;
    writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));

    chanMap0ind = int32(chanMap0ind);

    writeNPY(chanMap0ind, fullfile(savePath, 'channel_map.npy'));
    writeNPY([xcoords ycoords], fullfile(savePath, 'channel_positions.npy'));

    writeNPY(templateFeatureInds'-1, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
    writeNPY(pcFeatureInds'-1, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
    
    writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
    writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));

   

    if ~isempty(rez.cProj)
        writeNPY(rez.cProj, fullfile(savePath, 'template_features.npy'));
        writeNPY(rez.cProjPC, fullfile(savePath, 'pc_features.npy'));
    else
        NchanNear = rez.ops.min_Nnearest; 
        mfileF = memmapfile(rez.cProjpath, 'Format',{'single',[NchanNear, numel(rez.idiscard)],'x'});
        writeNPY(mfileF.Data.x(:,~rez.idiscard)', fullfile(savePath, 'template_features.npy'));
        clear mFileF;
        
        mfilePC = memmapfile(rez.cProjPCpath, 'Format',{'single',[NchanNear, 3, numel(rez.idiscard)],'x'});
        writeNPY(permute(mfilePC.Data.x(:,:,~rez.idiscard), [3 2 1]), fullfile(savePath, 'pc_features.npy'));
        clear mfilePC;
    end
    
    
     if isfield(rez, 'simScore')
        similarTemplates = rez.simScore;
        writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
    end


    % save a list of "good" clusters for Phy
%     fileID = fopen(fullfile(savePath, 'channel_names.tsv'), 'w');
%     fprintf(fileID, 'cluster_id%sKSLabel', char(9));
%     for j = 1:Nchan
%         fprintf(fileID, '%d%s%d', j-1,char(9),chanMap0ind(j));
%         fprintf(fileID, char([13 10]));
%     end
%     fclose(fileID);

    fileID = fopen(fullfile(savePath, 'cluster_KSLabel.tsv'),'w');
    fprintf(fileID, 'cluster_id%sKSLabel', char(9));
    fprintf(fileID, char([13 10]));

    fileIDCP = fopen(fullfile(savePath, 'cluster_ContamPct.tsv'),'w');
    fprintf(fileIDCP, 'cluster_id%sContamPct', char(9));
    fprintf(fileIDCP, char([13 10]));

    fileIDA = fopen(fullfile(savePath, 'cluster_Amplitude.tsv'),'w');
    fprintf(fileIDA, 'cluster_id%sAmplitude', char(9));
    fprintf(fileIDA, char([13 10]));

    rez.est_contam_rate(isnan(rez.est_contam_rate)) = 1;
    for j = 1:length(rez.good)
        if rez.good(j)
            fprintf(fileID, '%d%sgood', j-1, char(9));
        else
            fprintf(fileID, '%d%smua', j-1, char(9));
        end
        fprintf(fileID, char([13 10]));

        fprintf(fileIDCP, '%d%s%.1f', j-1, char(9), rez.est_contam_rate(j)*100);
        fprintf(fileIDCP, char([13 10]));

        fprintf(fileIDA, '%d%s%.1f', j-1, char(9), tempAmps(j));
        fprintf(fileIDA, char([13 10]));

    end
    fclose(fileID);
    fclose(fileIDCP);
    fclose(fileIDA);

     %make params file
    if ~exist(fullfile(savePath,'params.py'),'file')
        fid = fopen(fullfile(savePath,'params.py'), 'w');

        [~, fname, ext] = fileparts(rez.ops.fbinary);

        fprintf(fid,['dat_path = ''',fname ext '''\n']);
        fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
        fprintf(fid,'dtype = ''int16''\n');
        fprintf(fid,'offset = 0\n');
        if mod(rez.ops.fs,1)
            fprintf(fid,'sample_rate = %i\n',rez.ops.fs);
        else
            fprintf(fid,'sample_rate = %i.\n',rez.ops.fs);
        end
        fprintf(fid,'hp_filtered = False\n');
        fprintf(fid,'template_scaling = 20.0\n');
        fclose(fid);
    end
end

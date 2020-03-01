function savePrincipalComponents(rez, savepath_pc, savepath_feat)
%SAVEPRINCIPALCOMPONENTS Summary of this function goes here
%   Detailed explanation goes here

%==========================================================================
min_NchanNear = getOr(rez.ops, 'min_NchanNear', 32); 
NchanNear     = min(rez.ops.Nchan, min_NchanNear);
Nrank         = 3;
[mpath, ~]    = fileparts(rez.ops.fproc);

Nspikes       = size(rez.st3, 1);
Nunits        = max(rez.st3(:,2));
spkinds       = accumarray(rez.st3(:,2), 1:size(rez.st3,1),[],@(x) {x});
%==========================================================================
shape_pc    = [Nspikes Nrank NchanNear];
header_pc   = constructNPYheader('single', shape_pc);
shape_feat  = [Nspikes NchanNear];
header_feat = constructNPYheader('single', shape_feat);
Nbatch = 10; batchsize = ceil(Nspikes/Nbatch);
%==========================================================================
% make NPY files and fill with zeros
fid_pc = fopen(savepath_pc, 'W');
fwrite(fid_pc, header_pc, 'uint8');
for ibatch = 1:Nbatch
    swrite = min(batchsize, Nspikes - (ibatch-1)*batchsize);
    fwrite(fid_pc, zeros(NchanNear*Nrank*swrite, 1, 'single'), 'single');
end
fclose(fid_pc);

fid_feat = fopen(savepath_feat, 'W');
fwrite(fid_feat, header_feat, 'uint8');
for ibatch = 1:Nbatch
    swrite = min(batchsize, Nspikes - (ibatch-1)*batchsize);
    fwrite(fid_feat, zeros(NchanNear*swrite, 1, 'single'), 'single');
end
fclose(fid_feat);
%==========================================================================
% save PCs
mfilePC = memmapfile(savepath_pc, 'Format',{'single',[Nspikes, Nrank, NchanNear],'x'},...
    'Writable',true, 'Offset', numel(header_pc));

for iunit = 1:Nunits
    pcpath = fullfile(mpath, sprintf(rez.cProjPCpath, iunit));
	
	if ~exist(pcpath, 'file'), continue; end 
	
    fid_pc = fopen(pcpath, 'r');
    pcmat  = fread(fid_pc, '*single');
    fclose(fid_pc);
    pcmat = reshape(pcmat, [], Nrank, NchanNear);
    
    isp = sort(spkinds{iunit},'ascend');
    
    mfilePC.Data.x(isp,:,:) = pcmat;
    delete(pcpath);
end

clear mfilePC;
%==========================================================================
% save features

mfileFeat = memmapfile(savepath_feat, 'Format',{'single',[Nspikes, NchanNear],'x'},...
    'Writable',true, 'Offset', numel(header_feat));

for iunit = 1:Nunits
    featpath = fullfile(mpath, sprintf(rez.cProjpath, iunit));
	
	if ~exist(featpath, 'file'), continue; end 
	
    fid_feat   = fopen(featpath, 'r');
    featmat    = fread(fid_feat, '*single');
    fclose(fid_feat);
    
    featmat  = reshape(featmat, [], NchanNear);
    
    isp = sort(spkinds{iunit},'ascend');
    
    mfileFeat.Data.x(isp,:) = featmat;
    delete(featpath);
end

clear mfileFeat;


%==========================================================================

end


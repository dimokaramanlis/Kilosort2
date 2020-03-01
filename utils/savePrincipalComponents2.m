function savePrincipalComponents2(rez, savepath_pc, savepath_feat)
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
spkinds       = cellfun(@(x) sort(x,'ascend'), spkinds, 'un', 0);
nspks         = cellfun(@numel, spkinds);
%==========================================================================
shape_pc    = [Nspikes Nrank NchanNear];
header_pc   = constructNPYheader('single', shape_pc);
shape_feat  = [Nspikes NchanNear];
header_feat = constructNPYheader('single', shape_feat);
%==========================================================================
% save PC file

fid_pc = fopen(savepath_pc, 'W');
fwrite(fid_pc, header_pc, 'uint8'); % write header first
Ndims = Nrank * NchanNear;
Nbatch = 8; batchsize = ceil(Ndims/Nbatch);

for ibatch = 1:Nbatch

    dimswrite = min(batchsize, Ndims - (ibatch-1)*batchsize);
	istart    = (ibatch - 1)*batchsize + 1;
	iend      = istart + dimswrite - 1;
	
	pcwrite = zeros(Nspikes, dimswrite, 'single');
	
	for iunit = 1:Nunits
	    pcpath = fullfile(mpath, sprintf(rez.cProjPCpath, iunit));
		if ~exist(pcpath, 'file'), continue; end 
	
		fid_pc_unit = fopen(pcpath, 'r');
	    pcmat       = fread(fid_pc_unit, [nspks(iunit) Ndims],'*single');
		fclose(fid_pc_unit);
		
		pcwrite(spkinds{iunit}, :) = pcmat(:, istart:iend);
	end
	
    fwrite(fid_pc, pcwrite, 'single');	
	
end

fclose(fid_pc); clear pcwrite;
%==========================================================================
% save features file

fid_feat = fopen(savepath_feat, 'W');
fwrite(fid_feat, header_feat, 'uint8');
Ndims  = NchanNear;
Nbatch = 4; batchsize = ceil(Ndims/Nbatch);

for ibatch = 1:Nbatch

    dimswrite = min(batchsize, Ndims - (ibatch-1)*batchsize);
	istart    = (ibatch - 1)*batchsize + 1;
	iend      = istart + dimswrite - 1;
	
	featwrite = zeros(Nspikes, dimswrite, 'single');
	
	for iunit = 1:Nunits
	    featpath = fullfile(mpath, sprintf(rez.cProjpath, iunit));
		if ~exist(featpath, 'file'), continue; end 
	
		fid_feat_unit = fopen(featpath, 'r');
	    featmat       = fread(fid_feat_unit, [nspks(iunit) Ndims],'*single');
		fclose(fid_feat_unit);
		
		featwrite(spkinds{iunit}, :) = featmat(:, istart:iend);
	end
	
    fwrite(fid_feat, featwrite, 'single');	
	
end

fclose(fid_feat);
%==========================================================================
% delete single unit files
for iunit = 1:Nunits
	pcpath   = fullfile(mpath, sprintf(rez.cProjPCpath, iunit));
	featpath = fullfile(mpath, sprintf(rez.cProjpath, iunit));
    if exist(pcpath,   'file'), delete(pcpath);   end 
    if exist(featpath, 'file'), delete(featpath); end 
end
%==========================================================================

end


function deletePrincipalComponents(rez, iunit, idelete)

%--------------------------------------------------------------------------
min_NchanNear = getOr(rez.ops, 'min_NchanNear', 32); 
NchanNear     = min(rez.ops.Nchan, min_NchanNear);
Nrank  = 3;
[mpath, ~]     = fileparts(rez.ops.fproc);
pcpath   = fullfile(mpath, sprintf(rez.cProjPCpath, iunit));
featpath = fullfile(mpath, sprintf(rez.cProjpath,   iunit));
%--------------------------------------------------------------------------
% load pcs
fid_pc = fopen(pcpath, 'r');
pcmat  = fread(fid_pc, '*single');
fclose(fid_pc);

pcmat = reshape(pcmat, [], NchanNear*Nrank);
assert(size(pcmat, 1) == numel(idelete))
%--------------------------------------------------------------------------
% load feats
fid_feat = fopen(featpath, 'r');
featmat  = fread(fid_feat, '*single');
fclose(fid_feat);

featmat = reshape(featmat, [], NchanNear);
assert(size(featmat, 1) == numel(idelete))
%--------------------------------------------------------------------------
% delete PCs, feats and paths
pcmat  (idelete, :) = []; delete(pcpath); 
featmat(idelete, :) = []; delete(featpath);
%--------------------------------------------------------------------------
% write again
fid_pc = fopen(pcpath, 'W');
fwrite(fid_pc, pcmat, 'single');
fclose(fid_pc);

fid_feat = fopen(featpath, 'W');
fwrite(fid_feat, featmat, 'single');
fclose(fid_feat);
%--------------------------------------------------------------------------


end
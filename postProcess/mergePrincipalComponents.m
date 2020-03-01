function mergePrincipalComponents(rez, ilow, ihigh)
%MERGEPRINCIPALCOMPONENTS Merges files corresponding to single units
%--------------------------------------------------------------------------
% setup options
min_NchanNear = getOr(rez.ops, 'min_NchanNear', 32); 
NchanNear     = min(rez.ops.Nchan, min_NchanNear);
Nrank  = 3;
[mpath, ~]     = fileparts(rez.ops.fproc);
%--------------------------------------------------------------------------
% load spikes
stlow   = rez.st3(rez.st3(:,2) == ilow,  1);
sthigh  = rez.st3(rez.st3(:,2) == ihigh, 1);
%--------------------------------------------------------------------------
% load pcs
pathlow_pc = fullfile(mpath, sprintf(rez.cProjPCpath, ilow));
fidlow_pc  = fopen(pathlow_pc, 'r');
pclow   = fread(fidlow_pc,  [numel(stlow) NchanNear*Nrank ], '*single');
fclose(fidlow_pc);

pathhigh_pc = fullfile(mpath, sprintf(rez.cProjPCpath, ihigh));
fidhigh_pc = fopen(pathhigh_pc, 'r');
pchigh  = fread(fidhigh_pc, [ numel(sthigh) NchanNear*Nrank], '*single');
fclose(fidhigh_pc);
%--------------------------------------------------------------------------
% load features
pathlow_feat  = fullfile(mpath, sprintf(rez.cProjpath, ilow));
fidlow_feat   = fopen(pathlow_feat, 'r');
featlow       = fread(fidlow_feat,  [numel(stlow) NchanNear], '*single');
fclose(fidlow_feat);

pathhigh_feat = fullfile(mpath, sprintf(rez.cProjpath, ihigh));
fidhigh_feat  = fopen(pathhigh_feat, 'r');
feathigh      = fread(fidhigh_feat, [numel(sthigh) NchanNear], '*single');
fclose(fidhigh_feat);
%--------------------------------------------------------------------------
% merge spikes, pcs and features to get the correct order
stall     = [stlow; sthigh];
pcmerge   = [pclow; pchigh];
featmerge = [featlow; feathigh];

% sort pcs to keep the correct order
[~, isort]  = sort(stall,'ascend');
pcmerge     = pcmerge(isort, :);
featmerge   = featmerge(isort, :);
%--------------------------------------------------------------------------
% write merged pcs
fidhigh_pc = fopen(pathhigh_pc, 'W');
fwrite(fidhigh_pc, pcmerge, 'single');
fclose(fidhigh_pc);

% write merged features
fidhigh_feat = fopen(pathhigh_feat, 'W');
fwrite(fidhigh_feat, featmerge, 'single');
fclose(fidhigh_feat);

% delete pcs of ilow
delete(pathlow_pc); delete(pathlow_feat);
%--------------------------------------------------------------------------
end
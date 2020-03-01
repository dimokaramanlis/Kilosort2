function splitPrincipalComponents(rez, ilow, ihigh)
%--------------------------------------------------------------------------
% get options
min_NchanNear = getOr(rez.ops, 'min_NchanNear', 32); 
NchanNear     = min(rez.ops.Nchan, min_NchanNear);
Nrank         = 3;
[mpath, ~]    = fileparts(rez.ops.fproc);
%--------------------------------------------------------------------------
% find spk locs
spkinds = find(ismembc(rez.st3(:,2), [ihigh, ilow]));
assert(issorted(spkinds)); % important to check

indslow  = rez.st3(spkinds, 2)== ilow;
indshigh = rez.st3(spkinds, 2)== ihigh;
%--------------------------------------------------------------------------
% load pcs
pathhigh_pc = fullfile(mpath, sprintf(rez.cProjPCpath, ihigh));
fidhigh_pc = fopen(pathhigh_pc, 'r');
pchigh  = fread(fidhigh_pc, [numel(spkinds) NchanNear*Nrank], '*single');
fclose(fidhigh_pc);
%--------------------------------------------------------------------------
% load features
pathhigh_feat = fullfile(mpath, sprintf(rez.cProjpath, ihigh));
fidhigh_feat  = fopen(pathhigh_pc, 'r');
feathigh      = fread(fidhigh_feat, [numel(spkinds) NchanNear], '*single');
fclose(fidhigh_feat);
%--------------------------------------------------------------------------
% split in memory

pclow    = pchigh  (indslow,  :);
pchigh   = pchigh  (indshigh, :);
featlow  = feathigh(indslow,  :);
feathigh = feathigh(indshigh, :);

% delete high
delete(pathhigh_pc); delete(pathhigh_feat);
%--------------------------------------------------------------------------
% write split pcs
fidhigh_pc = fopen(pathhigh_pc, 'W');
fwrite(fidhigh_pc, pchigh, 'single');
fclose(fidhigh_pc);

pathlow_pc = fullfile(mpath, sprintf(rez.cProjPCpath, ilow));
fidlow_pc = fopen(pathlow_pc, 'W');
fwrite(fidlow_pc, pclow, 'single');
fclose(fidlow_pc);
%--------------------------------------------------------------------------
% write split features

fidhigh_feat = fopen(pathhigh_feat, 'W');
fwrite(fidhigh_feat, feathigh, 'single');
fclose(fidhigh_feat);

pathlow_feat = fullfile(mpath, sprintf(rez.cProjpath, ilow));
fidlow_feat  = fopen(pathlow_feat, 'W');
fwrite(fidlow_feat, featlow, 'single');
fclose(fidlow_feat);
%--------------------------------------------------------------------------
end
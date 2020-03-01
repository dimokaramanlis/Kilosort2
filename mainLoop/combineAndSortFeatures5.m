function rez = combineAndSortFeatures5(rez, fW_0, fWpc_0, fW_1, fWpc_1, nspk0, nspk1, st3)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 
% %==========================================================================
Nnearest = getOr(rez.ops, 'min_Nnearest', 32);
Nrank   = 3;

[mpath, ~]     = fileparts(fWpc_0);
Nunits = max(st3(:,2));
spkinds = accumarray(st3(:,2), 1:size(st3,1),[],@(x) {x});
nspks   = cellfun(@numel, spkinds);
Nbatch  = 6;

fprintf('Writing principal components and features to hard drive...\n');

%--------------------------------------------------------------------------
% go through first part in batches
fprintf('Going through first half...\n');

batchspikes = ceil(nspk0/Nbatch);

for ibatch = 1: Nbatch
    
    spksread = min(batchspikes, nspk0 - (ibatch-1)*batchspikes);
    istart = batchspikes * (ibatch-1) + 1;
    iend   = istart + spksread - 1;
    cspikes = st3(istart:iend, 2);
    
    fid_pc0 = fopen(fWpc_0, 'r');
    offset_pc = 4 * batchspikes * (ibatch-1) * Nnearest * Nrank;
    fseek(fid_pc0, offset_pc, 'bof');
    dat_pc = fread(fid_pc0, [Nnearest * Nrank spksread], '*single');
    fclose(fid_pc0);
    
    fid_feat0 = fopen(fW_0, 'r');
    offset_feat = 4 * batchspikes * (ibatch-1) * Nnearest;
    fseek(fid_feat0, offset_feat, 'bof');
    dat_feat = fread(fid_feat0, [Nnearest spksread], '*single');
    fclose(fid_feat0);
    
    for iunit = 1:Nunits
        unitpath_pc   = fullfile(mpath,   sprintf('pcs_%04d.dat', iunit));
        unitpath_feat = fullfile(mpath, sprintf('feats_%04d.dat', iunit));

        % change parts
        unitpcs   = dat_pc  (:, cspikes == iunit);
        unitfeats = dat_feat(:, cspikes == iunit);

        % write pcs
        fidunit_pc = fopen(unitpath_pc, 'A');
        fwrite(fidunit_pc, unitpcs, 'single');
        fclose(fidunit_pc);
        
        % write features
        fidunit_feat = fopen(unitpath_feat, 'A');
        fwrite(fidunit_feat, unitfeats, 'single');
        fclose(fidunit_feat);
    end

end
delete(fWpc_0); delete(fW_0); toc;
%--------------------------------------------------------------------------
% go through second part in batches
fprintf('Going through second half...\n');

batchspikes = ceil(nspk1/Nbatch);
for ibatch = 1: Nbatch
    
    spksread = min(batchspikes, nspk1 - (ibatch-1)*batchspikes);
    istart = nspk0 + batchspikes * (ibatch-1) + 1;
    iend   = istart + spksread - 1;
    cspikes = st3(istart:iend, 2);
    
    % read pcs
    fid_pc1 = fopen(fWpc_1, 'r');
    offset = 4 * batchspikes * (ibatch-1) * Nnearest * Nrank;
    fseek(fid_pc1, offset, 'bof');
    dat_pc = fread(fid_pc1, [Nnearest * Nrank spksread], '*single');
    fclose(fid_pc1);
    
    % read features
    fid_feat1 = fopen(fW_1, 'r');
    offset_feat = 4 * batchspikes * (ibatch-1) * Nnearest;
    fseek(fid_feat1, offset_feat, 'bof');
    dat_feat = fread(fid_feat1, [Nnearest spksread], '*single');
    fclose(fid_feat1);

    for iunit = 1:Nunits
        unitpath_pc   = fullfile(mpath, sprintf('pcs_%04d.dat', iunit));
        unitpath_feat = fullfile(mpath, sprintf('feats_%04d.dat', iunit));

        % change parts
        unitpcs   = dat_pc(:, cspikes == iunit);
        unitfeats = dat_feat(:, cspikes == iunit);
        
        % write pcs
        fidunit_pc = fopen(unitpath_pc, 'A');
        fwrite(fidunit_pc, unitpcs, 'single');
        fclose(fidunit_pc);
        
         % write features
        fidunit_feat = fopen(unitpath_feat, 'A');
        fwrite(fidunit_feat, unitfeats, 'single');
        fclose(fidunit_feat);
    end

end
delete(fWpc_1); delete(fW_1); toc;
%--------------------------------------------------------------------------
% sort and permute for each cell
fprintf('Sorting components for each cell...\n');
for iunit = 1:Nunits
    
    currinds   = sort(spkinds{iunit}, 'ascend');
    [~, csort] = sort(st3(currinds,1), 'ascend');
    
    %----------------------------------------------------------------------
    % PC PART
    
    unitpath_pc   = fullfile(mpath, sprintf('pcs_%04d.dat', iunit));

    % read
    fidunit_pc = fopen(unitpath_pc, 'r');
    unitpcs = fread(fidunit_pc,'*single');
    fclose(fidunit_pc);
    
    unitpcs = reshape(unitpcs, Nnearest, Nrank, []);
    assert( size(unitpcs, 3) == nspks(iunit));

    % change parts
    unitpcs    = permute(unitpcs, [3 2 1]);
    unitpcs    = unitpcs(csort,:,:);
    
    % write again
    fidunit_pc = fopen(unitpath_pc, 'W');
    fwrite(fidunit_pc, unitpcs, 'single');
    fclose(fidunit_pc);
    %----------------------------------------------------------------------
    % FEATURE PART
    
    unitpath_feat = fullfile(mpath, sprintf('feats_%04d.dat', iunit));
    
     % read
    fidunit_feat = fopen(unitpath_feat, 'r');
    unitfeats    = fread(fidunit_feat,'*single');
    fclose(fidunit_feat);
    
    unitfeats = reshape(unitfeats, Nnearest, []);
    assert( size(unitfeats, 2) == nspks(iunit));

    % change parts
    unitfeats  = unitfeats';
    unitfeats  = unitfeats(csort,:);
    
    % write again
    fidunit_feat = fopen(unitpath_feat, 'W');
    fwrite(fidunit_feat, unitfeats, 'single');
    fclose(fidunit_feat);
    %----------------------------------------------------------------------
end

toc;
%==========================================================================
rez.cProjPC = []; rez.cProj   = [];
rez.cProjPCpath = 'pcs_%04d.dat';
rez.cProjpath   = 'feats_%04d.dat';
%==========================================================================
end


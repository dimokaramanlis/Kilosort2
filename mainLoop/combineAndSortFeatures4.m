function rez = combineAndSortFeatures4(rez, isort, fW_0, fWpc_0, fW_1, fWpc_1, nspk0, nspk1, st3)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 
% %==========================================================================
Nnearest = getOr(rez.ops, 'min_Nnearest', 32);
Nrank   = 3;
%==========================================================================
fprintf('Writing template features to hard drive...\n');

[mpath, ~]     = fileparts(fWpc_0);
ntot = nspk0 + nspk1;
Nb = 10;
batchsize = ceil(ntot/Nb);

% deal with Features

mfileF0 = memmapfile(fW_0, 'Format',{'single',[Nnearest, nspk0],'x'});
mfileF1 = memmapfile(fW_1, 'Format',{'single',[Nnearest, nspk1],'x'});

fWfinalpath  = fullfile(mpath, 'fW.dat');
fid_fW_final = fopen(fWfinalpath, 'W');

for ii = 1:Nb
    swrite = min(batchsize, ntot - (ii-1)*batchsize);
    istart = (ii-1)*batchsize + 1;
    iend   = istart + swrite - 1;
    currinds = isort(istart:iend);
    
    fwtowrite = zeros(Nnearest, swrite, 'single');
    
    currinds0 = currinds(currinds<=nspk0);
    f0 = mfileF0.Data.x(:, currinds0);
    fwtowrite(:, currinds<=nspk0) = f0;
    
    currinds1 = currinds(currinds>nspk0) - nspk0;
    f1 = mfileF1.Data.x(:, currinds1);
    fwtowrite(:,currinds>nspk0) = f1;
    
    fwrite(fid_fW_final, fwtowrite, 'single');

end

fclose(fid_fW_final);

clear mfileF0; clear mfileF1; delete(fW_0); delete(fW_1); 
toc;
%==========================================================================
% Deal with PCs
fprintf('Writing principal components to hard drive...\n');

Nunits = max(st3(:,2));
spkinds = accumarray(st3(:,2), 1:size(st3,1),[],@(x) {x});
nspks   = cellfun(@numel, spkinds);
Nbatch  = 5;

%--------------------------------------------------------------------------
% go through first part in batches
fprintf('Going through first half...\n');

batchspikes = ceil(nspk0/Nbatch);

for ibatch = 1: Nbatch
    
    spksread = min(batchspikes, nspk0 - (ibatch-1)*batchspikes);
    
    fid1 = fopen(fWpc_0, 'r');
    offset = 4 * batchspikes * (ibatch-1) * Nnearest * Nrank;
    fseek(fid1, offset, 'bof');
    dat = fread(fid1, [Nnearest * Nrank spksread], '*single');
    fclose(fid1);

    istart = batchspikes * (ibatch-1) + 1;
    iend   = istart + spksread - 1;
    cspikes = st3(istart:iend, 2);
    
    for iunit = 1:Nunits
        unitpath = fullfile(mpath, sprintf('pcs_%04d.dat', iunit));
        
        % change parts
        unitpcs = dat(:, cspikes == iunit);

        % write again
        fidunit = fopen(unitpath, 'A');
        fwrite(fidunit, unitpcs, 'single');
        fclose(fidunit);
    end

end
delete(fWpc_0); toc;
%--------------------------------------------------------------------------
% go through second part in batches
fprintf('Going through second half...\n');

batchspikes = ceil(nspk1/Nbatch);
for ibatch = 1: Nbatch
    
    spksread = min(batchspikes, nspk1 - (ibatch-1)*batchspikes);
    
    fid1 = fopen(fWpc_1, 'r');
    offset = 4 * batchspikes * (ibatch-1) * Nnearest * Nrank;
    fseek(fid1, offset, 'bof');
    dat = fread(fid1, [Nnearest * Nrank spksread], '*single');
    fclose(fid1);

    istart = nspk0 + batchspikes * (ibatch-1) + 1;
    iend   = istart + spksread - 1;
    cspikes = st3(istart:iend, 2);
    for iunit = 1:Nunits
        unitpath = fullfile(mpath, sprintf('pcs_%04d.dat', iunit));
        
        % change parts
        unitpcs = dat(:, cspikes == iunit);

        % write again
        fidunit = fopen(unitpath, 'A');
        fwrite(fidunit, unitpcs, 'single');
        fclose(fidunit);
    end

end
delete(fWpc_1); toc;
%--------------------------------------------------------------------------
% sort and permute for each cell
fprintf('Sorting components for each cell...\n');
for iunit = 1:Nunits
    
    unitpath = fullfile(mpath, sprintf('pcs_%04d.dat', iunit));

    % read
    fidunit = fopen(unitpath, 'r');
    unitpcs = fread(fidunit,'*single');
    fclose(fidunit);
    
    unitpcs = reshape(unitpcs, Nnearest, Nrank, []);
    assert( size(unitpcs, 3) == nspks(iunit));

    % change parts

    unitpcs    = permute(unitpcs, [3 2 1]);
    currinds   = sort(spkinds{iunit}, 'ascend');
    [~, csort] = sort(st3(currinds,1), 'ascend');
    unitpcs    = unitpcs(csort,:,:);
    
    % write again
    fidunit = fopen(unitpath, 'W');
    fwrite(fidunit, unitpcs, 'single');
    fclose(fidunit);
end

toc;
%==========================================================================
rez.cProjPC = []; rez.cProj   = [];
rez.cProjPCpath = 'pcs_%04d.dat';
rez.cProjpath   = fWfinalpath;
%==========================================================================
end


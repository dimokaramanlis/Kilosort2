function rez = combineAndSortFeatures3(rez, isort, fW_0, fWpc_0, fW_1, fWpc_1, nspk0, nspk1, st3)
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

mfilePC0 = memmapfile(fWpc_0, 'Format',{'single',[Nnearest * Nrank, nspk0],'x'});
mfilePC1 = memmapfile(fWpc_1, 'Format',{'single',[Nnearest * Nrank, nspk1],'x'});


for iunit = 1:Nunits

	currinds = spkinds{iunit};
	fwpctowrite = zeros(Nnearest*Nrank, numel(currinds), 'single');
    
    currinds0 = currinds(currinds<=nspk0);
    pc0 = mfilePC0.Data.x(:, currinds0);
    fwpctowrite(:, currinds<=nspk0) = pc0;
    
    currinds1 = currinds(currinds>nspk0) -nspk0;
    pc1 = mfilePC1.Data.x(:, currinds1);
    fwpctowrite(:,currinds>nspk0) = pc1;
	
    [~, csort] = sort(st3(currinds,1), 'ascend');
    
	unitpath = fullfile(mpath, sprintf('pcs_%04d.dat', iunit));
    
	fwpctowrite = fwpctowrite(:, csort);
    fwpctowrite = permute(reshape(fwpctowrite, Nnearest, Nrank,  numel(currinds)), [3 2 1]);
    
	fidunit = fopen(unitpath, 'W');
    fwrite(fidunit, fwpctowrite, 'single');
	fclose(fidunit);
	
end
clear mfilePC0; clear mfilePC1; delete(fWpc_0); delete(fWpc_1); 
toc;
%==========================================================================
rez.cProjPC = []; rez.cProj   = [];
rez.cProjPCpath = 'pcs_%04d.dat';
rez.cProjpath   = fWfinalpath;
%==========================================================================
end


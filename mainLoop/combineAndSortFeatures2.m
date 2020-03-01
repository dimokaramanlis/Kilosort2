function rez = combineAndSortFeatures2(rez, isort, fW_0, fWpc_0, fW_1, fWpc_1, nspk0, nspk1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 
% %==========================================================================
Nnearest = getOr(rez.ops, 'min_Nnearest', 32);
Nrank   = 3;
%==========================================================================
% deal with PCs

[mpath, ~]     = fileparts(fWpc_0);
fWpcfinalpath  = fullfile(mpath, 'fWpc.dat');
ntot = nspk0 + nspk1;
Nb = 20;
batchsize = ceil(ntot/Nb);

mfilePC0 = memmapfile(fWpc_0, 'Format',{'single',[Nnearest * Nrank, nspk0],'x'});
mfilePC1 = memmapfile(fWpc_1, 'Format',{'single',[Nnearest * Nrank, nspk1],'x'});

fid_fWpc_final = fopen(fWpcfinalpath, 'W');

for ii = 1:Nb
    swrite = min(batchsize, ntot - (ii-1)*batchsize);
    istart = (ii-1)*batchsize + 1;
    iend   = istart + swrite - 1;
    currinds = isort(istart:iend);
    
    fwpctowrite = zeros(Nnearest*Nrank, swrite, 'single');
    
    currinds0 = currinds(currinds<=nspk0);
    pc0 = mfilePC0.Data.x(:, currinds0);
    fwpctowrite(:, currinds<=nspk0) = pc0;
    
    currinds1 = currinds(currinds>nspk0) -nspk0;
    pc1 = mfilePC1.Data.x(:, currinds1);
    fwpctowrite(:,currinds>nspk0) = pc1;

    
    fwrite(fid_fWpc_final, fwpctowrite, 'single');

end

fclose(fid_fWpc_final);
clear mfilePC0; clear mfilePC1; delete(fWpc_0); delete(fWpc_1); 
%==========================================================================

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
    
    currinds1 = currinds(currinds>nspk0) -nspk0;
    f1 = mfileF1.Data.x(:, currinds1);
    fwtowrite(:,currinds>nspk0) = f1;

    
    fwrite(fid_fW_final, fwtowrite, 'single');

end

fclose(fid_fW_final);
clear mfileF0; clear mfileF1; delete(fW_0); delete(fW_1); 
%==========================================================================
rez.cProjPC = []; rez.cProj   = [];
rez.cProjPCpath = fWpcfinalpath;
rez.cProjpath = fWfinalpath;
%==========================================================================
end


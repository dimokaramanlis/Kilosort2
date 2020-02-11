function rez = combineAndSortFeatures(rez, isort, fW_0, fWpc_0, fW_1, fWpc_1, nspk0, nspk1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nnearest = rez.ops.min_Nnearest;
Nrank   = 3;
%==========================================================================
% deal with PCs

[mpath, ~]     = fileparts(fWpc_0);
fWpcfinalpath  = fullfile(mpath, 'fWpc.dat');

% create final file
fid_fWpc_final = fopen(fWpcfinalpath, 'W');
ntot = nspk0 + nspk1;
nb = 20;
batchsize = ceil(ntot/nb);
for ii = 1:nb
    swrite = min(batchsize, ntot - (ii-1)*batchsize);
    fwrite(fid_fWpc_final, zeros(Nnearest*Nrank, swrite, 'single'), 'single');
end
fclose(fid_fWpc_final);

%memory map final file
mfilePC = memmapfile(fWpcfinalpath, 'Format',{'single',[Nnearest * Nrank, ntot],'x'}, 'Writable', true);

%write first part

fidpc0 = fopen(fWpc_0, 'r');
f0 = fread(fidpc0, [Nnearest*Nrank nspk0], '*single');
fclose(fidpc0);
f0 = f0(:, isort(isort<=nspk0));


mfilePC.Data.x(:, isort<=nspk0) = f0;
clear f0;

%write second part

fidpc1 = fopen(fWpc_1, 'r');
f1 = fread(fidpc1, [Nnearest*Nrank nspk1], '*single');
fclose(fidpc1);
f1 =  f1(:, isort(isort>nspk0)-nspk0);

mfilePC.Data.x(:, isort>nspk0) = f1;
clear f1;

delete(fWpc_0); delete(fWpc_1); clear mfilePC;
%==========================================================================
% deal with template features
[mpath, ~]     = fileparts(fW_0);
fWfinalpath  = fullfile(mpath, 'fW.dat');

% create final file
fid_fW_final = fopen(fWfinalpath, 'W');
ntot = nspk0 + nspk1;
nb = 10;
batchsize = ceil(ntot/nb);
for ii = 1:nb
    swrite = min(batchsize, ntot - (ii-1)*batchsize);
    fwrite(fid_fW_final, zeros(Nnearest, swrite, 'single'), 'single');
end
fclose(fid_fW_final);

%memory map final file
mfileF = memmapfile(fWfinalpath, 'Format',{'single',[Nnearest, ntot],'x'}, 'Writable', true);

%write first part

fid0 = fopen(fW_0, 'r');
f0 = fread(fid0, [Nnearest nspk0], '*single');
fclose(fid0);
f0 = f0(:, isort(isort<=nspk0));

mfileF.Data.x(:, isort<=nspk0) = f0;
clear f0;

%write second part

fid1 = fopen(fW_1, 'r');
f1 = fread(fid1, [Nnearest nspk1], '*single');
fclose(fid1);
f1 =  f1(:, isort(isort>nspk0)-nspk0);

mfileF.Data.x(:, isort>nspk0) = f1;
clear f1;

delete(fW_0); delete(fW_1); clear mfileF;

%==========================================================================


rez.cProjPC = [];
rez.cProj   = [];
rez.cProjPCpath = fWpcfinalpath;
rez.cProjpath = fWfinalpath;
%==========================================================================


end


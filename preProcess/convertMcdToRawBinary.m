function [bininfo] = convertMcdToRawBinary(ops, targetpath)
%CONVERTMCDTORAWBINARY Save mcd files as binary without analog channels
%Neuroshare functions and dlls have to be in MATLAB path
%We can potentially add what's needed in Kilosort's folder
%Also saves the time samples correponding to the beginning of each stimulus
%--------------------------------------------------------------------------
tic;
%get mcd filenames
mcdfilenames = dir([ops.root,filesep,'*.mcd']);
[~, reindex]=sort(str2double(regexp(({mcdfilenames(:).name}),'\d+','match','once')));
mcdfilenames={mcdfilenames(reindex).name}'; Nfiles=numel(mcdfilenames);
%--------------------------------------------------------------------------
%load the dll file
[dllpath,libtoload] = getMCSdllPath();
nsresult=mexprog(18, [dllpath, filesep, libtoload]);  %set dll library
%--------------------------------------------------------------------------
%get information about the recording time
filesamples = zeros(numel(mcdfilenames),1);
stimids = zeros(numel(mcdfilenames),1);

for imcd=1:numel(mcdfilenames)
    mcdpathname = [ops.root,filesep,mcdfilenames{imcd}]; %get mcd path
    [nsresult, hfile] = mexprog(1, mcdpathname); %open file
    [nsresult, mcdfileInfo] = mexprog(3, hfile); %get file info
    filesamples(imcd)=mcdfileInfo.TimeSpan/mcdfileInfo.TimeStampResolution;
    nsresult = mexprog(14, hfile);%close file
    stimids (imcd) = cell2mat(textscan(mcdfilenames{imcd}, '%d_'));
end
NchanTOT = mcdfileInfo.EntityCount;
filesamples = floor(filesamples);
%--------------------------------------------------------------------------
% make bininfo file for splitting later
bininfo.stimsamples = accumarray(stimids, filesamples,[],@sum);
fs = round(1/mcdfileInfo.TimeStampResolution); % sampling frequency
bininfo.fs=fs;
bininfo.NchanTOT = NchanTOT-4;
fprintf('Total length of recording is %2.2f min...\n', sum(filesamples)/fs/60);
%--------------------------------------------------------------------------
%get information about the array arrangement and the signal
[nsresult, hfile] = mexprog(1, [ops.root,filesep,mcdfilenames{1}]); %open file
[nsresult, chinfos] = mexprog(4, hfile,0:(NchanTOT-1)); %get channel info
[nsresult, volinfos] = mexprog(7, hfile,0:(NchanTOT-1)); % get general info
nsresult=mexprog(14, hfile);%close data file.
labellist = {chinfos.EntityLabel}; clear chinfos; %extract labels of the entities
maxVoltage=volinfos(1).MaxVal; minVoltage=volinfos(1).MinVal;
resVoltage=volinfos(1).Resolution; clear volinfos;
newRange=2^15*[-1 1]; multFact=range(newRange)/(maxVoltage-minVoltage);
%--------------------------------------------------------------------------
% get the channel names based on the map of the array
chanMap = getChannelMapForRawBinary(labellist,'dataformat','mcd','channelnumber',NchanTOT,'meatype',ops.meatype);
%--------------------------------------------------------------------------
fprintf('Saving .mcd data as .dat...\n');

maxSamples = 48e5;

fidOut= fopen(targetpath, 'W'); %using W (capital), makes writing ~4x faster
msg=[];
if fs == 1e4,   upsampfac = 3;  else,    upsampfac = 1;      end % upsample 10K to 30K     
for iFile=1:Nfiles
    mcdpathname = [ops.root,filesep,mcdfilenames{iFile}];
    nsamples = filesamples(iFile);
    Nchunk=ceil(nsamples/maxSamples);
    
    [nsresult, hfile] = mexprog(1, mcdpathname);  %open file
    for iChunk=1:Nchunk
        offset = max(0, (maxSamples * (iChunk-1)));
        sampstoload=min(nsamples-offset,maxSamples);
        
        [~,~,dat]=mexprog(8,hfile, chanMap, offset, sampstoload);%read data
        dat=int16(dat*multFact)';        
        if upsampfac>1, dat = kron(dat,ones(1,upsampfac,'int16')); end     % upsample 10K to 30K        

        
        fwrite(fidOut, dat, 'int16');
    end
    nsresult = mexprog(14, hfile); %close file
    
    %report status
    fprintf(repmat('\b', 1, numel(msg)));
    msg=sprintf('Time %3.0f min. Mcd files processed %d/%d \n',...
        toc/60, iFile,Nfiles);
    fprintf(msg);
    
end
fclose(fidOut); clear mexprog; %unload DLL
if upsampfac>1, bininfo.fs = fs*upsampfac; end % upsample 10K to 30K  
%--------------------------------------------------------------------------
end


function [dllpath,libtoload] = getMCSdllPath()
%GETMCSDLLPATH Summary of this function goes here

dlllocation = which('load_multichannel_systems_mcd');
dllpath = fileparts(dlllocation);

switch computer()
    case 'PCWIN'; libtoload = 'nsMCDLibraryWin32.dll';
    case 'GLNX86'; libtoload = 'nsMCDLibraryLinux32.so';
    case 'PCWIN64'; libtoload = 'nsMCDLibraryWin64.dll';
    case 'GLNXA64'; libtoload = 'nsMCDLibraryLinux64.so';
    case 'MACI64'; libtoload = 'nsMCDLibraryMacIntel.dylib';
    otherwise
        disp('Your architecture is not supported'); return;
end
end
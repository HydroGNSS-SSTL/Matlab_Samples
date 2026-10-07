%% ************************************************************************
% Read data from a HydroGNSS L1 DDM netCDF file
%% ************************************************************************
% MODULE NAME:      readL1DDM.m
% SOFTWARE NAME:    readL1DDM
% SOFTWARE VERSION: 1
%% ************************************************************************
% FUNCTIONS:
% Reads the scaled L1 DDM file
%% ************************************************************************
% USE:
% [l1_DDM] = readL1DDM('C:\PDGS_NAS_folder\HydroGNSS-1\DataRelease\L1A_L1B\2026-07\', 10, 10)
% folders_nc, Month folder name containing .nc files.
% Optional start_day and end_day required from the day/s of the month selected
%% ************************************************************************
% EXAMPLE:
% Read folders of DDMS from a start_day to end_day
% [l1_BB] = readL1Global(folder_nc,'blackbodyNadirL1E1LHCP.nc',start_day, end_day);
%% ************************************************************************
function [l1DDM] = readL1DDM(folders_nc, start_day, end_day)
disp('Reading L1 DDM file...')

%% input and settings
l1DDM = [];    % structure
l1DDM.count = 0;
l1DDM.Channel0 = [];
l1DDM.Channel1 = [];
l1DDM.Channel2 = [];
l1DDM.Channel3 = [];
l1DDM.stats.NumbOfTrack = 0; %numb of all the track in l1bdata
l1DDM.stats.NumbOfDDMS = 0; %total numb of DDMs
l1DDM.stats.numbOfTrackPerDataset = []; %vector containing the numb of tracks for each 6hr dataset present in the big dataset l1bdata
l1DDM.stats.datasetNumb = 0;  %numb of 6hr dataset contained in l1bdata

fileNameDDM = 'DDMs.nc';

% Code only when folders_nc is a folder and not a .nc file
if isfile(folders_nc)   %check if it exists as a file
    if contains(folders_nc,'nc')  %check if it is the  .nc FILE of interest
        fprintf('reading %s \n', folders_nc);
        [l1DDM] = readnc(folders_nc,l1DDM);
    else
        fprintf(' File %s not in the correct format', files_nc)
    end
else
    %for each month
    % monthlist = dir(folders_nc);
    % if isempty(monthlist); error('Month folders not found'); end
    % monthNames = {monthlist([monthlist.isdir]).name};
    % monthNames = monthNames(~ismember(monthNames ,{'.','..'}));
    % for m = 1:length(monthNames)

    % for each day of the month
    dayslist = dir(fullfile(folders_nc));
    if isempty(dayslist)
        error('Day folders not found');
    end
    daysNames = {dayslist([dayslist.isdir]).name};
    daysNames = daysNames(~ismember(daysNames ,{'.','..'}));
    daysNumbers = str2double(daysNames);

    % Wanted day checks
    if(isnan(daysNumbers))  % no data found
        error('File not in the right format, problem with daysNames reading');
    end
    if (nargin==1)
        start = 1;
        stop = length(daysNames);
    else % only selected days
        start = find(daysNumbers == start_day);
        stop = find(daysNumbers == end_day);

        % handling possible problems with days selection
        if(isempty(start))
            start=find(daysNumbers>start_day,1);
            if(isempty(start))
                error('No available startday');
            else
                disp('Selected start day not found, it is taken the first following day');
            end
        end
        if(isempty(stop))
            stop=find(daysNumbers>end_day,1);
            if(isempty(stop))
                error('No available enday');
            else
                disp('Selected end day not found,it is taken the first following day');
            end
        end
        if(stop<start)
            error('error with days selection');
        end
    end
    % Process days
    for day = start:stop
        hourlist = dir(fullfile(folders_nc,daysNames{day}));
        hourNames = {hourlist([hourlist.isdir]).name};
        hourNames = hourNames(~ismember(hourNames ,{'.','..'}));
        for h = 1:1:length(hourNames)   %for all the hours in a day
            tempfolder = fullfile(folders_nc, daysNames{day},hourNames{h});
            if (exist(fullfile(tempfolder,fileNameDDM),'file') == 2) % check there are corresponding files
                %% Read DDMs -----------------------------------------------
                fprintf('Reading DDMS in day %s hour %s\n',daysNames{day}, hourNames{h});
                dmfile = fullfile(folders_nc, daysNames{day},hourNames{h},fileNameDDM);   %DDMS.nc file opening
                [l1DDM] = readnc(dmfile,l1DDM);

            else
                fprintf('no %s in day %d hour %d, no ddms will be read \n',fileRead, day, h);
            end    % check if there is the correspondent mergedmetada

        end %for day hours

    end  %end  for day
end

fprintf('Finished Reading L1 DDMs.\n')
end

%% Function to process single files and add to l1data set
function [l1DDM] = readnc(files_nc,l1DDM)

ids = []; % reset for each hour of the day  % IDs, open the file and extract each parameter dataset by index
% DATASIZE, List all parameters and attributes, use first track as a guide
ids.ncid = netcdf.open(files_nc, 'NC_NOWRITE'); % Return ID of named group.
ids.globalNcids = netcdf.inqVarIDs(ids.ncid);
ni = ncinfo(files_nc); %info about the file
ids.trackNcids = netcdf.inqGrps(ids.ncid);  % tracks for each GPS and gal satellite
l1DDM.stats.datasetNumb = l1DDM.stats.datasetNumb + 1;      %number of 6hrs dataset taken into consideration
l1DDM.stats.numbOfTrackPerDataset(l1DDM.stats.datasetNumb) = length(ids.trackNcids);  %numb of tracks in each 6hr dataset

for track = 1:length(ids.trackNcids)  % for each track track 000000 for example
    ids.channelNcids{track} = netcdf.inqGrps(ids.trackNcids(track));  %select the channels, can be 4 or less ( ch0: Nadir_L1E1_LHCP, 1: Nadir_L1E1_RHCP,  2: Nadir_L5E5_LHCP,  3: Nadir_L5E5_RHCP)
    for chan = 1:length(ids.channelNcids{track})    %for each channel
        l = size(netcdf.inqGrps(ids.channelNcids{track}(chan)),2);   %l=2 can represent the fact that there are coherent and incoherent measurements
        ids.coinNcids{track}(chan,1:l) = netcdf.inqGrps(ids.channelNcids{track}(chan));
        if l == 1        % to solve problem of size compatibility
            ids.coinNcids{track}(chan,2)=NaN;
        end
    end
end

for inco = 1:length(ids.coinNcids) % for all the reflection tracks
    l1DDM.stats.NumbOfTrack = l1DDM.stats.NumbOfTrack+1;    %total number of tracks for all the days, hours and months
    if inco==round(length(ids.coinNcids)/2)
        disp('    processing...');
    end   %at half of the processing
    l1DDM.count = l1DDM.count + 1;
    for channel = 1:1: size(ids.channelNcids{inco},2)   % for all the 4 or less channels in each track  (they are not ordered)

        % DDM interrogation, initial bounding box analysis
        varTemp = netcdf.inqVarID(ids.coinNcids{inco}(channel,1), 'DDM'); % DDMS data per track and chan
        rawDDMuint  = (netcdf.getVar(ids.coinNcids{inco}(channel,1), varTemp));
        rawDDM = double(netcdf.getVar(ids.coinNcids{inco}(channel,1), varTemp));  % N raw ddms for a single channel in a single track
        ChanName=ni.Groups(inco).Groups(channel).Name;

        FilteredrawDDMuint = rawDDMuint(:,:,:);   
        FilteredrawDDM = rawDDM(:,:,:);
        l1DDM.stats.NumbOfDDMS = l1DDM.stats.NumbOfDDMS+size(FilteredrawDDMuint,3);  %total numb of DDMS for chan of interest associated to or water/land/ice/all content for all hours,days

        %find peaks or filtered ddms properties
        maxValue = 65535;  % max value for a 16bits ADC
        maskMaxVal = (FilteredrawDDMuint == maxValue);
        c = find(maskMaxVal);  % return the index by unrolling the matrix into a column vector
        [row,col,ddmIdx] = ind2sub(size(FilteredrawDDMuint),c); %returns the index associated to c in an index of row and column and z for a 3d matrix (3 by 3 matrix), there can be equal ddms indeces if a ddm contains more than 1 peak

        % find indeces of DDM with multiple peaks, consider only 1 index for 1 max in a ddm, removes multiple indeces associated to the same ddm
        [~,IA,~] = intersect(ddmIdx,1:1:length(FilteredrawDDMuint(1,1,:)),'stable');

        %allocation of vector for info about position of max peak in a ddm
        DDMpeakpos=zeros(size(FilteredrawDDMuint,3),2);
        if length(IA)~=size(FilteredrawDDMuint,3)
            fprintf('problem numb of ddms per track not = to n. of peaks= %d',maxValue)
            missIdx=setdiff(IA,1:1:length(FilteredrawDDMuint(1,1,:)));  %find missing ddms

            DDMpeakpos(IA,:) =  [row(IA,:),col(IA,:)];
            DDMpeakpos(missIdx,:) = 0;   %fill missing ddms with a zero value

        else
            DDMpeakpos =  [row(IA,:),col(IA,:)];
        end

        if(~ isempty(FilteredrawDDM))   %if there are still ddms in this track after the filters
            %% save raw DDM
            l1DDM.(ChanName).DDMs{l1DDM.count,1} = FilteredrawDDMuint; %cat(DIM,A,B) concatenates the arrays A and B along the dimension DIM.
            %l1bDDM.(ChanName).SVN(l1bDDM.count,1) = SVN{NumbOfTrack}(1);  %SVN
            % peak position
            l1DDM.(ChanName).DDMpeakpos{l1DDM.count,1} = DDMpeakpos;   %index of rows and columns where there is a peak

        else
            %fprintf('no filtered ddms for track %d \n',NumbOfTrack )
        end
    end   %end for channels

end  %end for tracks (inco)

netcdf.close(ids.ncid)  %close netcdf

end

%% ************************************************************************
% Read all data from a netCDF file
%% ************************************************************************
% MODULE NAME:      readL1Global.m
% SOFTWARE NAME:    readL1Global
% SOFTWARE VERSION: 1
%% ************************************************************************
% FUNCTIONS
% Reads the global parameters from a L1 NetCDF file
%% ************************************************************************
% USE
% [l1Data] = readL1Global(files_nc, file_type, startDay, endDay)
% files_nc, Month folder name containing .nc files.
% File type allows for selection of direct or blackbody files
% Optional start_day and end_day are the day/s of the month selected
% Associated DDMs and mergedL1 file is required in the same directory.
%% ************************************************************************
% EXAMPLE:
% Read L1 direct data in day 10 of July
% [l1_Direct] = readL1Global('C:\PDGS_NAS_folder\HydroGNSS-1\DataRelease\L1A_L1B\2026-07\10\H12','directSignalPowerL1E1.nc',10, 10); %
% Read the stated file in the current working directory
% [l1_blackbody] = readL1Global('blackbodyNadirL1E1LHCP.nc');
%% ************************************************************************
function [l1Data] = readL1Global(folders_nc,file_type,start_day, end_day)

% list of fields to ignore if causing problems, size changes for example
ignoreFields = {''};
% checks
if ~exist("file_type",'var')
    file_type = 'directSignalPowerL1E1.nc'; %default to L1 direct
end
fprintf('Reading L1 %s...',file_type)

l1Data = [];  % final structure
l1Data.Global = [];
if contains(folders_nc,'.nc')
    [l1Data] = readnc(folders_nc,ignoreFields, 1, l1Data);
else
    %Iterate over a directory of files
    %     monthlist = dir(files_nc);
    %     if isempty(monthlist); error('Month folders not found'); end
    %     monthNames = {monthlist([monthlist.isdir]).name};
    %     monthNames = monthNames(~ismember(monthNames ,{'.','..'}));
    %     %[indx,tf] = listdlg('PromptString',{'Select a directory to process:',...
    %     %    'Month of dataset',''},...
    %     %    'ListString',monthNames);
    %     %if tf == 0
    %     %    error('No file selected');
    %     %end
    %     %monthNames = monthNames(indx);
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
    for d = start:stop
        hourlist = dir(fullfile(folders_nc,daysNames{d}));
        hourNames = {hourlist([hourlist.isdir]).name};
        hourNames = hourNames(~ismember(hourNames ,{'.','..'}));
        for h = 1:length(hourNames)
            filenc = fullfile(folders_nc,daysNames{d},hourNames{h},file_type);
            [l1Data] = readnc(filenc,ignoreFields, h, l1Data);
        end
    end
end
disp('Finished.')
end

%% Function to process single files and add to l1data set
function [l1bData] = readnc(files_nc,~, ~, l1bData)

ncfile = fullfile(files_nc);
catdim = 1; % dimension to concatenate along

%% DATASIZE, List all parameters and attributes, use first track as a guide
ni = ncinfo(ncfile);
if ~isempty(ni.Variables) % merged data
    data.globalVariables = string({ni.Variables.Name}'); % global
end
l1bData.globalAttributes = {ni.Attributes.Name;ni.Attributes.Value}';

%% IDs, open the file and extract each parameter dataset by index
ids.ncid = netcdf.open(ncfile, 'NC_NOWRITE');
% obtain 'NCIDs' for each track, channel and data type
ids.globalNcids = netcdf.inqVarIDs(ids.ncid);

%% EXTRACT, then address each track or Co/Incoherent variable depending on the track and channel required:
if ~isempty(ni.Variables) % NON-merged data
    % fill global
    for globalIdx = 1:length(data.globalVariables) % global
        nmTemp = data.globalVariables(globalIdx); % name
        if(~(isfield(l1bData.Global, (nmTemp))))  % allocation
            l1bData.Global.(nmTemp) =[];
        end
        varTemp = netcdf.inqVarID(ids.ncid, nmTemp); % get data from name
        datasetsVal = netcdf.getVar(ids.ncid, varTemp);
        if ndims(datasetsVal) == 3 
            catdim = 3; % dimension to concatenate along
        end
        l1bData.Global.(nmTemp) = cat(catdim,l1bData.Global.(nmTemp),datasetsVal);
    end
end

netcdf.close(ids.ncid)
end
%% ************************************************************************
% Read all data from a netCDF file
%% ************************************************************************
% MODULE NAME:      readL1merged.m
% SOFTWARE NAME:    readL1merged
% SOFTWARE VERSION: 1
%% ************************************************************************
% AUTHORS:      T. Norris
% @copyright:   SSTL
%% ************************************************************************
% FUNCTIONS
% Accumulates the data in one huge vector per field for easy trending 
% and filtering
%% ************************************************************************
% USE
% Either provide a folder name to extract all data, or a single file name.
% Run by command line: l1bData = readL1merged('metadata_L1_merged.nc');
%% ************************************************************************
function [l1bData] = readL1merged(files_nc)
disp('Reading L1merged...')

% list of fields to ignore if causing problems, size changes for example
ignoreFields = {'IntegratedDelay','IntegratedDoppler','Temp Axis'};

l1bData = [];
if contains(files_nc,'.nc')
    [l1bData] = readnc(files_nc,ignoreFields, 1, l1bData);
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
   % for m = 1:length(monthNames)
%         dayslist = dir(fullfile(files_nc));
%         if isempty(dayslist); error('Day folders not found'); end
%         daysNames = {dayslist([dayslist.isdir]).name};
%         daysNames = daysNames(~ismember(daysNames ,{'.','..'}));
%         for d = 1:length(daysNames)
%             %if ~isnan(str2double(dayslist(d).name))
%             hourlist = dir(fullfile(files_nc,daysNames{d}));
%             hourNames = {hourlist([hourlist.isdir]).name};
%             hourNames = hourNames(~ismember(hourNames ,{'.','..'}));
%             for h = 1:length(hourNames)
%                 filenc = fullfile(files_nc,daysNames{d},hourNames{h},'metadata_L1_merged.nc');
                if isfile(filenc)
                    [l1bData] = readnc(filenc,ignoreFields, (d*4) + (h-4), l1bData);
                else 
                    %break
                end
%             end
%         end
% 
%     %end
end
disp('Finished Reading L1merged.')
end

%% Funciton to process single files and add to l1data set
function [l1bData] = readnc(files_nc,ignoreFields, count, l1bData)

ncfile = fullfile(files_nc);

%% DATASIZE, List all parameters and attributes, use first track as a guide
ni = ncinfo(ncfile);
if ~isempty(ni.Variables) % merged data
    data.globalVariables = string({ni.Variables.Name}'); % global
end
l1bData.globalAttributes = {ni.Attributes.Name;ni.Attributes.Value}';
data.trackNumber = length(ni.Groups); % tracks
data.trackVariables = string({ni.Groups(1).Variables.Name}');
data.chanNumber = length(ni.Groups(1).Groups); % channels
data.chanVariables = string({ni.Groups(1).Groups(1).Variables.Name}');
data.inCoVariables = string({ni.Groups(1).Groups(1).Groups(1).Variables.Name}'); % incoherent

%% IDs, open the file and extract each parameter dataset by index
ids.ncid = netcdf.open(ncfile, 'NC_NOWRITE');
% obtain 'NCIDs' for each track, channel and data type
ids.globalNcids = netcdf.inqVarIDs(ids.ncid);
ids.trackNcids = netcdf.inqGrps(ids.ncid);
for track = 1:length(ids.trackNcids)
    ids.channelNcids(track,:) = netcdf.inqGrps(ids.trackNcids(track));
    for chan = 1:length(ids.channelNcids(track,:))

        ids.coinNcids{track}(chan,:) = netcdf.inqGrps(ids.channelNcids(track,chan));
    end
end

%% EXTRACT, then address each track or Co/Incoherent variable depending on the track and channel required:

% preallocate
if count==1
    if ~isempty(ni.Variables) % NON-merged data
    % fill global
        for globalIdx = 1:length(data.globalVariables) % global
            nmTemp = strrep(data.globalVariables(globalIdx),' ','_'); % name
            l1bData.Global.(nmTemp) = [];
        end
    end
    for trackIdx = 1:length(data.trackVariables) % track
        % attributes
        trackattributes = {ni.Groups(trackIdx).Attributes.Name}';
        for att = 1:length(trackattributes)
            l1bData.TrackAttributes.(trackattributes{att}) = [];
        end
        % variables
        nmTemp = strrep(data.trackVariables(trackIdx),' ','_'); % variable
        l1bData.Track.(nmTemp) = [];
        for channel = 1:size(ids.channelNcids,2)
            % Chan attributes
            chanattributes = {ni.Groups(1).Groups(1).Attributes.Name}';
            for att = 1:length(chanattributes)
                l1bData.ChanAttributes.(chanattributes{att}) = [];
            end
            % Inco attributes
            if ~isempty(ni.Variables) % NON-merged data
                incoattributes = {ni.Groups(trackIdx).Groups(channel).Groups(1).Attributes.Name}';
                for att = 1:length(incoattributes)
                    l1bData.InCoAttributes.(incoattributes{att}) = [];
                end
            end
            % variables
            for inco = 1:length(data.inCoVariables)
                nmTemp = strrep(data.inCoVariables(inco),' ','_'); % name
                l1bData.(sprintf('IncoherentChan%d',channel)).(nmTemp) = [];
            end
        end
    end
end

if ~isempty(ni.Variables) % NON-merged data
    % fill global
    for globalIdx = 1:length(data.globalVariables) % global
        nmTemp = strrep(data.globalVariables(globalIdx),' ','_'); % name
        varTemp = netcdf.inqVarID(ids.ncid, data.globalVariables(globalIdx)); % get data from name
        l1bData.Global.(nmTemp) = [l1bData.Global.(nmTemp); netcdf.getVar(ids.ncid, varTemp)];
    end
end

% fill all
for trackIdx = 1:data.trackNumber % track 000000 for example
    %  track variables
    for trackVar = 1:length(data.trackVariables) % each track variable
        % read and store all track parameters
        nmTemp = strrep(data.trackVariables(trackVar),' ','_'); % name
        varTemp = netcdf.inqVarID(ids.trackNcids(trackIdx), data.trackVariables(trackVar)); % data
        dataarray = netcdf.getVar(ids.trackNcids(trackIdx), varTemp);
        l1bData.Track.(nmTemp) = [l1bData.Track.(nmTemp);dataarray];
    end
    % track attributes
    trackattributes = {ni.Groups(trackIdx).Attributes.Name}';
    for att = 1:length(trackattributes)
        l1bData.TrackAttributes.(trackattributes{att}) = ...
            [l1bData.TrackAttributes.(trackattributes{att}); repmat(ni.Groups(trackIdx).Attributes(att).Value,length(dataarray),1)];
    end
    for channel = 1:size(ids.channelNcids,2)
        % channel variables
        %nmTemp = data.trackVariables(trackVar); % name
        %varTemp = netcdf.inqVarID(ids.trackNcids(trackIdx), nmTemp); % data
        %l1bData.Track.(nmTemp) = [l1bData.Track.(nmTemp);netcdf.getVar(ids.trackNcids(trackIdx), varTemp)];
        if ~isempty(ni.Variables) % NON-merged data
            % chan attributes
            chanattributes = {ni.Groups(trackIdx).Groups(channel).Attributes.Name}';
            for att = 1:length(chanattributes)
                l1bData.ChanAttributes.(chanattributes{att}) = [l1bData.ChanAttributes.(chanattributes{att});ni.Groups(trackIdx).Groups(channel).Attributes(att).Value];
            end
            % inco attributes
            if ~isempty(ni.Variables) % NON-merged data
                incoattributes = {ni.Groups(trackIdx).Groups(channel).Groups(1).Attributes.Name}';
                for att = 1:length(incoattributes)
                    l1bData.InCoAttributes.(incoattributes{att}) = [l1bData.InCoAttributes.(incoattributes{att});ni.Groups(trackIdx).Groups(channel).Groups(1).Attributes(att).Value];
                end
            end
        end
        for inco = 1:length(data.inCoVariables) % read and store all incoherent parameters
            % read and store all channel parameters
            nmTemp = strrep(data.inCoVariables(inco),' ','_'); % name
            if any(strcmp(data.inCoVariables(inco),ignoreFields))
                % field ignored
            else
                varTemp = netcdf.inqVarID(ids.coinNcids{trackIdx}(channel,1), data.inCoVariables(inco)); % data
                l1bData.(sprintf('IncoherentChan%d',channel)).(nmTemp) = ...
                    [l1bData.(sprintf('IncoherentChan%d',channel)).(nmTemp);netcdf.getVar(ids.coinNcids{trackIdx}(channel,1), varTemp)];
            end
        end
    end
end

netcdf.close(ids.ncid)
end
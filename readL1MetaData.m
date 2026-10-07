%% ************************************************************************
% Read data from a HydroGNSS L1 merged metadata with H5library (low level)
%% ************************************************************************
% MODULE NAME:      readL1MetaData.m
% SOFTWARE NAME:    readL1MetaData
% SOFTWARE VERSION: 1
%% ************************************************************************
% FUNCTIONS
% Reads the L1 merged metadata file
%% ************************************************************************
% USE
% [l1bData] = readL1MetaData(files_nc, startDay, endDay, outputFields, preservetracks, combineconstellations)
% files_nc, Month folder name containing .nc files.
% Optional startDay and endDay are the day/s of the month selected
% outputFields, only extract parameters defined in this cell list, empty if all parameters are to be extracted
% preservetracks, bool
% combineconstellations, bool
% Associated DDMs and mergedL1 file is required in the same directory.
%% ************************************************************************
% EXAMPLE:
% Read a specfied file, maintain track and channel structure in output:
% [l1b_Data_file] = readL1MetaData(fullfile(folder_nc,'metadata_L1_merged.nc'), 1, 1, [], 1, 1);
%
% Read day 6 of a month (Concatenate multiple hourly files), with output groups tagged with channel names
% [l1b_Data_month] = readL1MetaData(folder_nc, 6, 6, [], 0, 0)
%% ************************************************************************
function [l1bData] = readL1MetaData(files_nc, ...
    start_Day, end_Day, output_fields, preserve_tracks, combine_constellations)

% init and checks
count = 0;       % numb of 6hr datasets
l1bData = [];  % final structure
l1bData.Global = [];
endpath = '.nc';

% list of fields to ignore that have caused problems in the past, if present
ignoreFields = {'IntegratedDelay','IntegratedDoppler','Temp Axis', 'DelayAxis', 'DopplerAxis'};

if ~exist('outputFields', 'var')
    output_fields = [];
end
if ~exist('preservetracks', 'var')
    preserve_tracks = 0; % default to concatenated list
end
if ~exist('combineconstellations', 'var')
    combine_constellations = 0; % default to output channels rather than constellation specific structures
end
tagFields = {'PRN','SVN','GnssBlock','MasterNavChannelNumber'}; % list of fields from track parameters to add to the channel output

if isfile(files_nc)   %check if it exists as a file
    if contains(files_nc,endpath)  %check if it is the  .nc FILE of interest
        fprintf('reading %s \n', files_nc);
        [l1bData] = readnc(files_nc,ignoreFields, tagFields, output_fields, l1bData, preserve_tracks,combine_constellations);
        fprintf('\nReading %s file...\n',endpath);
    else
        fprintf(' File %s not in the correct format', files_nc)
    end

else  %check if it is a folder
    filetype = 'metadata_L1_merged';
    if (~isfolder(files_nc))    % if folder does not exist it returns a NaN structure 
        disp(' Input directory does not exist.' );
        [l1bData] = NaN;

    else   % %folder exists: start reading files in the folder

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
        dayslist = dir(fullfile(files_nc));
        if isempty(dayslist); error('Day folders not found'); end
        daysNames = {dayslist([dayslist.isdir]).name};
        daysNames = daysNames(~ismember(daysNames ,{'.','..'}));

        if (nargin == 1 )   % if no inputs are given, all the days are read
            start=1;
            stop=length(daysNames);

        else
            daysNumbers=str2double(daysNames);
            %Days check
            if(isnan(daysNumbers))  
                error('Days not found or not in the right format'); 
            end
            start=find(daysNumbers==start_Day);
            stop=find(daysNumbers==end_Day);

            if(isempty(start))  
                start=find(daysNumbers>start_Day,1);
                if(isempty(start)) 
                    error('No available startday');
                else 
                    disp('Selected start day not found, using the first successive day');
                end
            end

            if(isempty(stop))  
                stop=find(daysNumbers>end_Day,1);
                if(isempty(stop)) 
                    error('No available endDay');
                else 
                    disp('Selected end day not found, using the next successive day');
                end
            end
            if(stop<start)
                error('error');
            end

        end
        for d=start:stop  % for the days selected
            hourlist = dir(fullfile(files_nc,daysNames{d}));
            hourNames = {hourlist([hourlist.isdir]).name};
            hourNames = hourNames(~ismember(hourNames ,{'.','..'}));
            for h = 1:length(hourNames)
                filenc = fullfile(files_nc,daysNames{d},hourNames{h},[filetype,endpath]);
                if (isfile(filenc) && contains(filenc,endpath))
                    count=count+1;  % numb of 6hr dataset read
                    str1 = sprintf('Reading day %d track ', daysNumbers(d));
                    str2=hourNames(h);
                    disp(append(str1, str2));
                    [l1bData] = readnc(filenc,ignoreFields, tagFields, output_fields, ...
                        l1bData, preserve_tracks, combine_constellations);  
                else
                    disp('File is not valid.')
                    % break
                end
            end
        end

    end
end
fprintf('\n')
disp('Finished Reading L1 dataset.')
end

%% Function to read a single file contained in files_nc  and to concatenate all the files read in the l1bdata5 structure in input,
% reading implemented by using the h5 library

function [l1bDataH5] = readnc(files_nc,ignoreFields, tagFields, outputFields, l1bDataH5, preservetracks,combineconstellations)
fid = H5F.open(files_nc, 'H5F_ACC_RDONLY', 'H5P_DEFAULT'); % Return ID

%info about the variables datasets attributes
%freqPolStr_vect={};   %vect of names of all channels available per track
stats=[];  %structure with all info about names and num of fields
stats.groups.num =0 ;
stats.datasets.num = 0;
stats.attributes.num = 0;
stats.datasetsPerTrack.num =0;
%nameGr={};  %name for groups
chanName={}; %name for channels
%VariablesPerchan=[];
% start from root group
path='/';
Gidfile = H5G.open(fid, path);
info = H5G.get_info(Gidfile);   %global group info
n = info.nlinks;   %numb of groups /00000 (tracks) + numb of datasets variables

% track loop
for k = 0:n-1 %numb of groups /00000 (tracks) + numb of datasets variables
    stopflag=0;  %flag to stop loop wether there are errors
    stats.datasetsPerTrack.num =0;
    stats.datasetsPerTrack.names={};
    stats.chan.num{k+1}=0;  %numb of channel per track
    %find name by idx
    name = H5L.get_name_by_idx(fid, path, ...
        'H5_INDEX_NAME','H5_ITER_INC', ...
        k,'H5P_DEFAULT');

    num = str2double(name);  % Convert  name string in number
    % if it is not a number -> global datasets variable name
    if isnan(num)
        %datasets Variables
        %       stats.datasets.num =stats.datasets.num +1;       %numb of variables
        %       stats.datasets.names{stats.datasets.num}=name;   % names
        nmTemp = strrep(name,' ','_');

        if(~(isfield(l1bDataH5.Global, (nmTemp))))  % allocation
            l1bDataH5.Global.(nmTemp) =[];
        end
        datasetsVal=H5D.read(H5D.open(fid,name)); %reading data values
        l1bDataH5.Global.(nmTemp)=[l1bDataH5.Global.(nmTemp);datasetsVal ];

        %if num is a number -> track
    else
		if k == 0
			fprintf('Track %4d',k)
		elseif k == n-1
			fprintf('\b\b\b\b%4d \n',k)
		else
			fprintf('\b\b\b\b%4d',k)
		end
        stats.groups.num = stats.groups.num + 1;
        fullpath = [path name];  %eg fullpath= /00000
        gid=H5G.open(fid,fullpath);  %group info per track
        infoGr=H5G.get_info(gid);
        nameGr={};
        %numb of groups + numb of datasets per track
        for g=0:infoGr.nlinks-1
            nameGr{g+1}=H5L.get_name_by_idx(gid,'.','H5_INDEX_NAME','H5_ITER_NATIVE',g,'H5P_DEFAULT');  % names of groups +channels
        end
        chnIdx=find(contains(nameGr, 'Channel'));        %idx of channel
        chanName(k+1, 1:size(chnIdx,2))=nameGr(chnIdx);  %chan names per track
        stats.chan.num{k+1}=size(chnIdx,2);  %numb of chan per track

        %% for each channel
        for ich=1:size(chnIdx,2)  %numb of channels per track
            Nvar=0;  %numb of dataValues per each channel
            chanstr = chanName{k+1, ich};
            fullpathChan=[fullpath, path,chanstr];
            
            gidCh=H5G.open(fid,fullpathChan);
            %attributes reading from the name
            if combineconstellations
                freqPolStr = chanstr; % use file channel names
                chanID = H5A.read(H5A.open(gidCh, 'ChannelID'));
            else
                freq = H5A.read(H5A.open(gidCh, 'SignalFrequency'));
                pol = H5A.read(H5A.open(gidCh, 'SignalPolarisation'));
                if iscell(freq); freq = freq{1}; end
                if iscell(pol); pol = pol{1}; end
                freqPolStr=sprintf('%s_%s',freq,pol);  %freqPolStr_vect{k+1, ich};
            end
            if ~(isfield(l1bDataH5, (freqPolStr)))   %allocation
                l1bDataH5.(freqPolStr)=[];
            end
            fullpathIncChan=[fullpathChan, path,'Incoherent'];
            gidperChan=H5G.open(fid,fullpathIncChan);
            infoGr=H5G.get_info(gidperChan);
            stats.VarPerchan.num=infoGr.nlinks;   %numb of datasets variables per channel
            nameVarperChan=[];
            %datasets variables per track
            datasetIdx=find(~contains(nameGr, 'Channel') & ~contains(nameGr, 'IntegrationMidPointTime') ); %idx
            stats.datasetsPerTrack.names=nameGr(datasetIdx);  %names of variables
            stats.datasetsPerTrack.num =size(datasetIdx,2);   %numb of variables
            
            %% for each dataset in a channel
            if ~isempty(outputFields)
                [LIA,~] = ismember(stats.datasetsPerTrack.names,outputFields); % check for wanted parameters
                dataSetTrack = stats.datasetsPerTrack.names(LIA);
            else
                dataSetTrack = stats.datasetsPerTrack.names;
            end
            if combineconstellations
                if ~(isfield(l1bDataH5.(freqPolStr),'ChannelID')) 
                    l1bDataH5.(freqPolStr).ChannelID = []; % hardware channel 1 to 16 in L1
                    l1bDataH5.(freqPolStr).ChannelID{k+1,1} = double(chanID);
                else
                        % cell array of tracks
                        l1bDataH5.(freqPolStr).ChannelID{k+1,1} = double(chanID);
                end
            end
            for idat=1:length(dataSetTrack)             
                Namedat= dataSetTrack{idat};
                nmTempDat = strrep(Namedat,' ','_');
                pathGrdataset=[fullpath, path, Namedat];
                if ~(isfield(l1bDataH5.(freqPolStr), (nmTempDat))) && ~preservetracks %allocation
                    l1bDataH5.(freqPolStr).(nmTempDat)=[];
                end               
                %check possible errors in data variables values 
                try
                    datasetArray=H5D.read(H5D.open(fid,pathGrdataset));
                    dim=size(datasetArray,1);
                catch ME
                    fprintf('error in datasetsvar track %d \n', k)
                    stopflag=1;
                    break;
                end
                % empty array
                if dim==0
                    if Nvar==0
                        stopflag=1;
                        break;
                    else
                        datasetArray=nan(Nvar,1);
                    end
                end
                % right numb of values to allocate per channel in track k
                while(Nvar==0)
                    Nvar=dim;
                end
                if preservetracks
                    % cell array of tracks
                    l1bDataH5.CommonVariables.(nmTempDat){k+1,1} = datasetArray;
                else
                    % concatenate tracks
                    l1bDataH5.(freqPolStr).(nmTempDat) ...
                        = [l1bDataH5.(freqPolStr).(nmTempDat) ; datasetArray];
                end
            end
            % error in values allocation -> skip to next track to have consistence
            if stopflag
                fprintf('track %d not considered \n',k)
                break;
            end

            %% for each group and freq and pol find the attributes belonging to interestingFields
            for ifield=1:size(tagFields,2)
                % allocate interesting fields (attributes)
                if ~(isfield(l1bDataH5.(freqPolStr), (tagFields{ifield}))) && ~preservetracks
                    l1bDataH5.(freqPolStr).(tagFields{ifield}) =[];
                end
                % check problems
                try
                    valAtt=H5A.read(H5A.open(gid,tagFields{ifield}));
                    AttributesSat= repmat(valAtt,Nvar,1);
                catch ME
                    fprintf('error in %s -> filled with Nan values.\n',tagFields{ifield});
                    AttributesSat= nan(Nvar,1);
                end
                if preservetracks
                    % cell array of tracks
                    l1bDataH5.CommonVariables.(tagFields{ifield}){k+1,1} = AttributesSat;
                else
                    % concatenate tracks
                    l1bDataH5.(freqPolStr).(tagFields{ifield}) ...
                        = [l1bDataH5.(freqPolStr).(tagFields{ifield}); AttributesSat ];
                end
            end

            %% for each variable varying per each channel
            for varCh=0:stats.VarPerchan.num-1
                nameVarperChan{varCh+1} = H5L.get_name_by_idx(fid, fullpathIncChan, ...
                    'H5_INDEX_NAME','H5_ITER_INC', ...
                    varCh,'H5P_DEFAULT');
                if ~isempty(outputFields) % check if the parameter is required in the output
                    matchFields = any(contains(outputFields, nameVarperChan{varCh+1}));
                else
                    matchFields = true;
                end

                if (matchFields && ~any(contains(ignoreFields,nameVarperChan{varCh+1})))

                    fullpathIncChanVar=[fullpathIncChan, path, nameVarperChan{varCh+1}];
                    if(~(isfield(l1bDataH5.(freqPolStr), (nameVarperChan{varCh+1}))))   %allocation
                        l1bDataH5.(freqPolStr).(nameVarperChan{varCh+1})=[];
                    end
                    VariablesPerchan=H5D.read(H5D.open(fid,fullpathIncChanVar));
                    % check error in size correspondence for each channel
                    if size(VariablesPerchan,1)~= Nvar
                        VariablesPerchan( size(VariablesPerchan,1)+1, Nvar)=NaN;
                        fprintf('prob Variablesperchan size -> filled with Nan values.\n')
                    end
                    if preservetracks
                        % cell array of tracks
                        l1bDataH5.(freqPolStr).(nameVarperChan{varCh+1}){k+1,1} ...
                            = VariablesPerchan;
                    else
                        % concatenate tracks
                        l1bDataH5.(freqPolStr).(nameVarperChan{varCh+1}) ...
                            = [l1bDataH5.(freqPolStr).(nameVarperChan{varCh+1}); VariablesPerchan];
                    end
                else
                    %ignorefields
                end
            end

        end

    end

end
%close H5
H5F.close(fid);
end

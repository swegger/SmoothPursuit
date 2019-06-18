function d = SetupSmoothPursuitProject(sname,projectName,varargin)
%% SetupSmoothPursuitProject
%
%   d = SetupSmoothPursuitProject(sname)
%
%   Sets of a data structure for smooth pursuit. Returns data structure
%   with the following fields:
%       sname ? subject name
%       projectDir - project directory
%       meta - meta data about each data collection run
%
%%

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'sname')
addRequired(Parser,'projectName')

parse(Parser,sname,projectName,varargin{:})

sname = Parser.Results.sname;
projectName = Parser.Results.projectName;

%% Project info
directory = ['~/Projects/' projectName];

%% Check if project file for sname exists; if not create 
cd(directory);
test = exist(sname,'file');
if ~test
    mkdir(sname);
end

%% Check if project data file exists for sname; if not create
cd([directory '/' sname])
datafile = [directory '/' sname '/' sname '_' projectName '.mat'];
test = exist(datafile,'file');
if ~test
    d.sname = sname;
    d.projectName = projectName;
    d.projectFnc = eval(['@' projectName]);
    d.projectDir = directory;
    save(datafile,'-struct','d')
    
else
    d = load(datafile);
end

%% Determine if any data files exist that have not been processed
files = dir([d.projectDir '/' d.sname '/']);
datafileInds = find(vertcat(files.isdir));

if length(datafileInds) > 2
    for i = 3:length(datafileInds)
        irun = i-2;
        fileName = files(datafileInds(i)).name;
        
        if isfield(d,'meta')
            if length(d.meta) < irun
                dataDir = [d.projectDir '/' sname '/' fileName];
                d.meta(irun).datafile = dataDir;
                
                d = d.projectFnc(irun,d,fileName);
            elseif ~strcmp(d.meta(irun).datafile(end-7:end),fileName)
                dataDir = [d.projectDir '/' sname '/' fileName];
                d.meta(irun).datafile = dataDir;
                
                d = d.projectFnc(irun,d,fileName);
            end
        else
            dataDir = [d.projectDir '/' sname '/' fileName];
            d.meta(irun).datafile = dataDir;
            d = d.projectFnc(irun,d,fileName);
        end
    end
else
    disp(['No data files for subject ' sname ', project ' projectName])
end

%% Save output
save(datafile,'-struct','d','-v7.3')



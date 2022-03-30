%% trainingObject
%
%   Defines properties and methods of object for analysis of training dta.
%
%%

classdef trainingObject
    % trainingObject class
    properties (SetAccess = private)
        sname;
        datapath;
        trials = 0;
        trialNumbers;
        trialDataFiles;
        trialtype;
        directions;
        speeds;
        onOff;
        eye;
        calib;
    end
    
    methods
        
        %% Core methods
        function obj = trainingObject(sname,datapath)
        %Constructor
            obj.sname = sname;
            obj.datapath = datapath;
            
            % For affine transformation of eye signal data
            obj.calib.t = 1:150;
            obj.calib.posGain = 0.025;
            obj.calib.speedGain = 0.09189;
            obj.calib.accThres = 1.1;
            obj.calib.windowSize = 40;
            
        end
        
        function out = returnDatapath(obj)
            out = obj.datapath;
        end
        
        
        function GetSize(this)
            props = properties(this);
            totSize = 0;
            
            for ii=1:length(props)
                currentProperty = getfield(this, char(props(ii)));
                s = whos('currentProperty');
                totSize = totSize + s.bytes;
            end
            
            fprintf(1, '%d bytes\n', totSize);
        end
        
        %% Analysis methods
        function [condInds, condLogical] = trialSort(obj,directions,speeds,locations)
        % Find indices of all trials with direction  in directions, speed
        % in speeds, location in locations.
            if isnan(directions)
                dMask = true(size(obj.directions));
            else
                dMask = ismember(obj.directions,directions);
            end
            
            if isnan(speeds)
                sMask = true(size(obj.speeds));
            else
                sMask = ismember(obj.speeds,speeds);
            end
                                    
            condLogical = prod([dMask,sMask],2,'native');
            condInds = find(condLogical);
        end
                
        function obj = trainingTrials(obj,trialNs)
        % data from training paradigm
        
            % Add direction preference trials to data object
            files = dir(obj.datapath);
            
            % Determine the index of the first data file
            for fileInx = 1:length(files)
                if length(files(fileInx).name) >= length(obj.sname) && ...
                        strcmp(files(fileInx).name(1:length(obj.sname)),obj.sname)
                    break
                end
            end
            ind = 0;
                for ti = trialNs
                    
                    if length(files)-fileInx+1 > ti
                        
                        % Read file
                        file = readcxdata([obj.datapath '/' files(ti+fileInx-1).name]);
                        
                        trialname = file.trialname;
                        ind = ind+1;
                        % Update trial
                        obj.trialNumbers(ind,1) = ti;
                        obj.trialDataFiles{ind} = files(ti+fileInx-1).name;
                        
                        % Parse trial info
                        obj.directions(ind,1) = str2double(...
                            extractBefore(file.trialname,'-'));
                        obj.speeds(ind,1) = str2double(...
                            extractAfter(file.trialname,'-'));
                        obj.onOff(ind,:) = file.targets.on{1};
                        
                        % Add eye information
                        obj.eye(:,ind).hpos = (file.data(1,:))*obj.calib.posGain;
                        obj.eye(:,ind).vpos = (file.data(2,:))*obj.calib.posGain;
                        
                        obj.eye(:,ind).hvel = (file.data(3,:))*obj.calib.speedGain;
                        obj.eye(:,ind).vvel = (file.data(4,:))*obj.calib.speedGain;
                        
                        sacs = saccadeDetect(file.data(3,:)*obj.calib.speedGain,...
                            file.data(4,:)*obj.calib.speedGain,...
                            'accelerationThreshold',obj.calib.accThres,...
                            'windowSize',40);
                        obj.eye(:,ind).hvel(sacs) = NaN;
                        obj.eye(:,ind).vvel(sacs) = NaN;
                        obj.eye(:,ind).saccades = sacs;
                                                
                    end
                end
        end
        
        
        function [mE,steE,E] = MeanEyeSpeed(obj,condInds,varargin)
        % Plots mean eye speed for a set of trials specified in condLogical
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addRequired(Parser,'condInds')
            addParameter(Parser,'t',NaN)
            addParameter(Parser,'normalizer',1)
            
            parse(Parser,obj,condInds,varargin{:})
            
            obj = Parser.Results.obj;
            condInds = Parser.Results.condInds;
            t = Parser.Results.t;
            normalizer = Parser.Results.normalizer;
            
            E = (sqrt(vertcat(obj.eye(condInds).hvel).^2 + ...
                vertcat(obj.eye(condInds).vvel).^2 ))/normalizer;
            mE = nanmean(E,1);
            steE = sqrt(nanvar(E,[],1)/length(condInds));
            
        end
        
        %% Plotting methods
        function h = plotEyeTraces(obj,directions,speeds)
        % Plot the eye traces for each direction and speed specified by
        % directions and speeds
            
            h = figure;
            ind = 0;
            for di = 1:length(directions)
                for si = 1:length(speeds)
                    ind = ind+1;
                    subplot(length(directions),length(speeds),ind)
                    [condInds,~] = trialSort(obj,directions(di),speeds(si),NaN);
                    for i = 1:length(condInds)
                        plot(obj.eye(condInds(i)).vvel(obj.onOff(condInds(i),1):obj.onOff(condInds(i),2)),...
                            'Color',[1 0 0 5/length(condInds)])
                        hold on
                        plot(obj.eye(condInds(i)).hvel(obj.onOff(condInds(i),1):obj.onOff(condInds(i),2)),...
                            'Color',[0 0 0 5/length(condInds)])
                    end
                    xlim([0 1800])
                    ylim([-max(speeds)*1.25 max(speeds)*1.25])
                    lineProps.color = 'k';
                    lineProps.LineStyle = '--';
                    plotHorizontal(cosd(directions(di))*speeds(si),'lineProperties',lineProps);
                    lineProps.color = 'r';
                    plotHorizontal(sind(directions(di))*speeds(si),'lineProperties',lineProps);
                    
                    if di == length(directions)
                        xlabel('Time from motion onset (ms)')
                    end
                    if si == 1
                        ylabel('Eye velocity (deg/s)')
                    end
                    title(['Dir = ' num2str(directions(di)) ', Speed = ' num2str(speeds(si))])
                end
            end
        end
        
        function h = plotMeanEyeVelocity(obj,condLogical,varargin)
        % Plots mean eye velocity for a set of trials specified in condLogical
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addRequired(Parser,'condLogical')
            addParameter(Parser,'h',NaN)
            addParameter(Parser,'sh',NaN)
            addParameter(Parser,'t',NaN)
            addParameter(Parser,'color',NaN)
            
            parse(Parser,obj,condLogical,varargin{:})
            
            obj = Parser.Results.obj;
            condLogical = Parser.Results.condLogical;
            h = Parser.Results.h;
            sh = Parser.Results.sh;
            t = Parser.Results.t;
            color = Parser.Results.color;
            
            if ishandle(h)
                figure(h);
            else
                h = figure;
            end
            if ishandle(sh)
                subplot(sh)
            end
            E(:,:,1) = vertcat(obj.eye(~~condLogical).hvel);
            E(:,:,2) = vertcat(obj.eye(~~condLogical).vvel);
            mE = nanmean(E,1);
            steE = sqrt(nanvar(E,[],1)/sum(condLogical));
            
            if isnan(t)
                t = 0:size(mE,2)-1;
            end
            
            patchProps.FaceAlpha = 0.3;
            if ~any(isnan(color))
                patchProps.FaceColor = color;
            else
                color = [0 0 0];
                patchProps.FaceColor = color;
            end
            Etemp = mE(:,:,1);
            steEtemp = steE(:,:,1);
            myPatch(t(:),Etemp(:),steEtemp(:),'patchProperties',patchProps);
            hold on
            Etemp = mE(:,:,2);
            steEtemp = steE(:,:,2);
            myPatch(t(:),Etemp(:),steEtemp(:),'patchProperties',patchProps);
            plot(t,mE(:,:,1),'Color',color,'LineWidth',2)
            plot(t,mE(:,:,2),'--','Color',color,'LineWidth',2)
        end  
        
    end
end

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
        motionOnOff;
        eye;
        calib;
        rotation;
        reward;
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
        function [condInds, condLogical] = trialSort(obj,directions,speeds,undoRotation)
        % Find indices of all trials with direction  in directions, speed
        % in speeds, location in locations.
            if ~exist('undoRotation','var')
                undoRotation = false;
            end
            if isnan(directions)
                dMask = true(size(obj.directions));
            else
                if undoRotation
                    dMask = ismember(obj.directions+obj.rotation,directions);
                else
                    dMask = ismember(obj.directions,directions);
                end
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
                        obj.rotation(ind,1) = file.key.iVelTheta/1000;
                        obj.reward(ind,1) = file.key.iRewLen1;
                        obj.directions(ind,1) = str2double(...
                            extractBefore(file.trialname,'-'));
                        obj.speeds(ind,1) = str2double(...
                            extractAfter(file.trialname,'-'));
                        obj.onOff(ind,:) = file.targets.on{1};
                        patSpeed = sqrt(file.targets.patvelH.^2+file.targets.patvelV.^2);
                        tSpeed = sqrt(file.targets.hvel.^2 + file.targets.vvel.^2);
                        obj.motionOnOff(ind,:) = [find(patSpeed(1,:),1) find(tSpeed(1,:),1,'last')];
                        
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
        
        
        function [vE,hE,vEmean,hEmean,vEste,hEste] = MeanEyeSpeed(obj,varargin)
        % Computes mean eye speed for a set of trials specified in condInds
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'directions',NaN)
            addParameter(Parser,'speeds',NaN)
            addParameter(Parser,'onTime',-250)
            addParameter(Parser,'offTime',0)
            addParameter(Parser,'maxTrials',200)
            addParameter(Parser,'maxDir',3000)
            addParameter(Parser,'undoRotation',true)
            
            parse(Parser,obj,varargin{:})
            
            obj = Parser.Results.obj;
            directions = Parser.Results.directions;
            speeds = Parser.Results.speeds;
            onTime = Parser.Results.onTime;
            offTime = Parser.Results.offTime;
            maxTrials = Parser.Results.maxTrials;
            maxDir = Parser.Results.maxDir;
            undoRotation = Parser.Results.undoRotation;
            
            if any(isnan(directions))
                if undoRotation
                    directions = unique(obj.directions+obj.rotation);
                else
                    directions = unique(obj.directions);
                end
            end
            
            if any(isnan(speeds))
                speeds = unique(obj.speeds);
            end
            
            vE = nan(maxDir,maxTrials,length(directions),length(speeds));
            hE = nan(maxDir,maxTrials,length(directions),length(speeds));
            for si = 1:length(speeds)
                for di = 1:length(directions)
                    [condInds,~] = trialSort(obj,directions(di),speeds(si),undoRotation);
                    
                    for i = 1:length(condInds)
                        timePoints = length(obj.eye(condInds(i)).vvel( (obj.onOff(condInds(i),1)+onTime) : (obj.motionOnOff(condInds(i),2))+offTime) );
                        vE(1:timePoints,i,di,si) = obj.eye(condInds(i)).vvel( (obj.onOff(condInds(i),1)+onTime) : (obj.motionOnOff(condInds(i),2)+offTime) );
                        hE(1:timePoints,i,di,si) = obj.eye(condInds(i)).hvel( (obj.onOff(condInds(i),1)+onTime) : (obj.motionOnOff(condInds(i),2)+offTime) );
                    end
                    nTrials(1,1,di,si) = length(condInds);
                end
            end
            
            vEmean = nanmean(vE,2);
            vEste = sqrt(nanvar(vE,[],2)./repmat(nTrials,[maxDir,maxTrials,1,1]));
            
            hEmean = nanmean(hE,2);
            hEste = sqrt(nanvar(hE,[],2)./repmat(nTrials,[maxDir,maxTrials,1,1]));
            
        end
        
        %% Plotting methods
        function h = plotEyeTraces(obj,directions,speeds,varargin)
        % Plot the eye traces for each direction and speed specified by
        % directions and speeds
        
            % Pasrse inputs
            Parser = inputParser;
            
            addRequired(Parser,'obj')
            addRequired(Parser,'directions')
            addRequired(Parser,'speeds')
            addParameter(Parser,'NormalizeTime',false)
            addParameter(Parser,'sampleN',NaN)
            addParameter(Parser,'colorAlpha',NaN)
            addParameter(Parser,'undoRotation',true)
            
            parse(Parser,obj,directions,speeds,varargin{:})
            
            obj = Parser.Results.obj;
            directions = Parser.Results.directions;
            speeds = Parser.Results.speeds;
            NormalizeTime = Parser.Results.NormalizeTime;
            sampleN = Parser.Results.sampleN;
            colorAlpha = Parser.Results.colorAlpha;
            undoRotation = Parser.Results.undoRotation;
                        
            % Make figure
            h = figure('Name','Eye traces by direction and speed','Units','Normalized','Position',[0.1310 0.0340 0.4377 0.8597]);
            ind = 0;
            for si = 1:length(speeds)
                for di = 1:length(directions)
                    ind = ind+1;
                    subplot(length(speeds),length(directions),ind)
                    [condInds,~] = trialSort(obj,directions(di),speeds(si),undoRotation);
                    
                    if isnan(sampleN)
                        sampleN = length(condInds);
                    elseif sampleN < length(condInds)
                        condInds = randsample(condInds,sampleN);
                    end
                    
                    if isnan(colorAlpha)
                        colorAlpha = 1/length(condInds);
                    end
                    
                    for i = 1:length(condInds)
                        if NormalizeTime
                            plot(linspace(0,1,obj.motionOnOff(condInds(i),2)-obj.motionOnOff(condInds(i),1)+251),...
                                obj.eye(condInds(i)).vvel(obj.motionOnOff(condInds(i),1):obj.motionOnOff(condInds(i),2)+250),...
                                'Color',[1 0 0 colorAlpha])
                            hold on
                            plot(linspace(0,1,obj.motionOnOff(condInds(i),2)-obj.motionOnOff(condInds(i),1)+251),...
                                obj.eye(condInds(i)).hvel(obj.motionOnOff(condInds(i),1):obj.motionOnOff(condInds(i),2)+250),...
                                'Color',[0 0 1 colorAlpha])
                        else
                            plot(obj.eye(condInds(i)).vvel(obj.onOff(condInds(i),1):obj.onOff(condInds(i),2)),...
                                'Color',[1 0 0 colorAlpha])
                            hold on
                            plot(obj.eye(condInds(i)).hvel(obj.onOff(condInds(i),1):obj.onOff(condInds(i),2)),...
                                'Color',[0 0 1 colorAlpha])
                        end
                    end
                    if NormalizeTime
                        xlim([0 1])
                    else
                        xlim([0 2250])
                    end
                    ylim([-max(speeds)*1.25 max(speeds)*1.25])
                    lineProps.color = 'b';
                    lineProps.LineStyle = '--';
                    plotHorizontal(cosd(directions(di))*speeds(si),'lineProperties',lineProps);
                    lineProps.color = 'r';
                    plotHorizontal(sind(directions(di))*speeds(si),'lineProperties',lineProps);
                    
                    if di == floor(length(directions)/2) && si == length(speeds) && NormalizeTime
                        xlabel('Normalized time')
                    elseif di == floor(length(directions)/2) && si == length(speeds)
                        xlabel('Time from motion onset (ms)')
                    end
                    if di == 1
                        ylabel('Eye velocity (deg/s)')
                    end
                    title(['Dir = ' num2str(directions(di)) ', Speed = ' num2str(speeds(si))])
                end
            end
        end
        
        function h = plotMeanEyeVelocity(obj,varargin)
        % Plots mean eye velocity for a set of trials specified in condLogical
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'directions',NaN)
            addParameter(Parser,'speeds',NaN)
            addParameter(Parser,'onTime',-250)
            addParameter(Parser,'offTime',0)
            addParameter(Parser,'maxTrials',200)
            addParameter(Parser,'maxDir',3000)
            addParameter(Parser,'vEmean',NaN)
            addParameter(Parser,'hEmean',NaN)
            addParameter(Parser,'undoRotation',true)
            
            parse(Parser,obj,varargin{:})
            
            obj = Parser.Results.obj;
            directions = Parser.Results.directions;
            speeds = Parser.Results.speeds;
            onTime = Parser.Results.onTime;
            offTime = Parser.Results.offTime;
            maxTrials = Parser.Results.maxTrials;
            maxDir = Parser.Results.maxDir;
            vEmean = Parser.Results.vEmean;
            hEmean = Parser.Results.hEmean;
            undoRotation = Parser.Results.undoRotation;
            
            if any(isnan(directions))
                if undoRotation
                    directions = unique(obj.directions+obj.rotation);
                else
                    directions = unique(obj.directions);
                end
            end
            
            if any(isnan(speeds))
                speeds = unique(obj.speeds);
            end
            
            if any(isnan(vEmean(:))) | any(isnan(vHmean(:)))
                [~,~,vEmean,hEmean,~,~] = MeanEyeSpeed(obj,'directions',directions,'speeds',speeds,...
                    'onTime',onTime,'offTime',offTime,'maxTrials',maxTrials,'undoRotation',undoRotation);
            end
            t = onTime:size(vEmean,1)+onTime-1;
            
            figure('Name','Mean eye velocities','Units','Normalized','Position',[0.1310 0.0340 0.4377 0.8597])
            ind = 0;
            for si = 1:length(speeds)
                for di = 1:length(directions)
                    ind = ind+1;
                    subplot(length(speeds),length(directions),ind)
                    plot(t,hEmean(:,1,di,si),'b-')
                    hold on
                    plot(t,vEmean(:,1,di,si),'r-')
                    xlabel('Time (ms)')
                    ylabel('Velocity (deg/s)')
                    
                    xlim([onTime size(vEmean,1)+onTime])
                    ylim([-max(speeds)*1.25 max(speeds)*1.25])
                    lineProps.color = 'b';
                    lineProps.LineStyle = '--';
                    plotHorizontal(cosd(directions(di))*speeds(si),'lineProperties',lineProps);
                    lineProps.color = 'r';
                    plotHorizontal(sind(directions(di))*speeds(si),'lineProperties',lineProps);
                end
            end
        end  
        
    end
end

classdef MTAFet 


    properties (SetAccess = public)
        mode
        Features        
    end

    methods
        function Fet = MTAFet(Session,varargin)
            [label_mode,overwrite] = DefaultArgs(varargin,{'manual',0});
            if isa(Session,'MTASession'),
                Session = [Session.spath.analysis Session.filebase];
            end

            if ~exist([Session '.fet.' label_mode '.mat'],'file')||overwrite,
                Fet.mode = label_mode;
                Fet.Features  = {};
                Fet.save(Session,overwrite);
            else
                load([Session '.fet.' label_mode '.mat']);
            end
        end

        function save(Fet,Session,overwrite)
            if isa(Session,'MTASession'),
                Session = [Session.spath.analysis Session.filebase];
            end

            for i = 1:length(Fet.Features),
                if ~Fet.Features{i}.ifsave
                    Fet.Features{i}.feature = [];
                end
            end

            if ~exist([Session '.fet.' Fet.mode '.mat'],'file')
                save( [Session '.fet.' Fet.mode '.mat'],'Fet','-v7.3');
            elseif exist([Session '.fet.' Fet.mode '.mat'],'file')&&overwrite
                warning(['Overwriting: ' Session '.fet.' Fet.mode '.mat']);
                save( [Session '.fet.' Fet.mode '.mat'],'Fet','-v7.3');
            else
                warning(['File exists: ' Session '.fet.' Fet.mode '.mat', ' - flag the overwrite option  to save']);
            end
        end

        function Fet = addFeature(Fet,label,expression,ifsave)
            Fet.States(end+1) = {MTAFeature(label,expression,ifsave)};
        end

        function Feature = getFeature(Fet,label)
            for i = 1:length(Fet.Features),
                if strcmp(Fet.Features{i}.label,label), 
                    Feature = Fet.Features{i}; 
                    break;
                end
            end
        end


        %% make label(s) a cell array to be able to compute
        %% multiple features at once
        function Fet = computeFeatures(Fet,Session,labels)
            numOfLabels = length(labels);
            labelCount =  0;
            for i = 1:length(Fet.Features),
                if ~isempty(find(cellfun(@strcmp,labels,repmat({Fet.Features{i}.label},1,numOfLabels))==1))
                    Fet.Feature{i} = Fet.Features{i}.compute(Session); 
                    labelCount = labelCount+1;
                    if labelCount==numOfLabels,
                        break
                    end
                end
            end

        end


    end
end

classdef MTADang < MTAData
    properties 
        model
    end

    properties(Transient=true)
        data        % data
    end
    methods
        function Data = MTADang(varargin)
            [path,filename,data,sampleRate,model,type,ext] = ...
                DefaultArgs(varargin,{[],[],[],[],[],'TimeSeries','ang'});
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename = [filename '.' ext '.mat'];
                end
            end
            Data = Data@MTAData(path,filename,data,sampleRate,type,ext);
            Data.model = model;
        end
        
        function Data = create(Data,Session)
            ang = zeros(size(Session.xyz,1),Session.Model.N,Session.Model.N,5);
            diffMat = Session.markerDiffMatrix();
            for i=1:Session.Model.N,
                for j=1:Session.Model.N,
                    if i==j,
                        continue
                    end
                    [rz,~,direction] = rotZAxis(squeeze(diffMat(:,i,j,:)));
                    [ry,~,pitch ] = rotYAxis(rz);
                    ang(:,i,j,1) = direction;
                    ang(:,i,j,2) = pitch;
                    ang(:,i,j,3:5) = ry;
                end
            end
            Data.data = ang;
        end
        
        function Data = filter(Data,win)
        end
        function Data = resample(Data,newSampleRate,varargin)
            [interp_type] = DefaultArgs(varargin,{'linear'});
        end
        function Data = embed(Data,win,overlap)
        end

    end
    
end
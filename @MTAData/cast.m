function Data = cast(Data,type,varargin)
%function Data = cast(Data,type,varargin)
%
%
%  varargin:
%    [sampleRate,syncflag]
%
%    sampleRate: {numeric,MTAData}, modify the sampleRate of the
%                 returned object.
%   
%    syncflag:   string, Either 'relative' to the object's sync, 
%                zero starts at the beginig of the object's sync,
%                or absolute to the session, zero starts at the
%                begining of the recordings.
%
    
[sampleRate,syncflag] = DefaultArgs(varargin,{Data.sampleRate,'relative'});
nsData = {};
if isa(sampleRate,'MTAData'),    
    nsData = sampleRate;
    sampleRate = nsData.sampleRate;
    % Cast to size and sampleRate which matches nsData
    % - work in progress - never mind I guess I have to finish it now
end

oldType = Data.type;            
if ~strcmp(oldType,type)&&~isempty(oldType)
    switch type
      case 'TimePeriods'
        Data.data = ThreshCross(Data.data,0.5,0);
      case 'TimeSeries'

        if strcmp(Data.label,'sync')||isempty(Data.label),
            tmpdata = round(((Data.data)./Data.sampleRate+Data.origin).*sampleRate);
        else
            tmpdata = round((Data.data./Data.sampleRate+Data.sync.sync.data(1)).*sampleRate);
        end
        tmpdata(tmpdata==0) = 1;

% $$$                         data = false(round(Data.subsref(substruct('.','sync','()',...
% $$$                                     {prod(size(Data.subsref(substruct('.','sync','()',{':'}))))}))...
% $$$                                     .*sampleRate),1);
        
        if isa(Data.subsref(substruct('.','sync')),'MTAData'),
            data = false(round(Data.subsref(substruct('.','sync','.','data','()',...
                                                      {prod(size(Data.subsref(substruct('.','sync','()',{':'}))))}))...
                               .*sampleRate),1);
        else
            data = false(round(Data.subsref(substruct('.','sync','()',...
                                                      {prod(size(Data.subsref(substruct('.','sync','()',{':'}))))}))...
                               .*sampleRate),1);
        end                                

        % Flag each indice as included or excluded
        for j = 1:size(tmpdata,1),
            data(tmpdata(j,1):tmpdata(j,2)) = true;
        end         

        switch syncflag
          case 'relative'
            if isempty(nsData),
                try
                    data = [data(round(Data.sync.sync.data(1)*sampleRate)+1:round(Data.sync.sync.data(end)*sampleRate));0];
                catch err
                    data = cat(1,data,zeros(round(Data.sync.sync.data(end)*sampleRate)-round(Data.sync.data(end)*sampleRate),1));
                    data = data(round(Data.sync.sync.data(1)*sampleRate)+1:round(Data.sync.sync.data(end)*sampleRate));
                end
            else
                data = [data(round(Data.sync.sync.data(1)*sampleRate)+1:...
                             round(Data.sync.sync.data(1)*sampleRate)+nsData.size(1))];
            end
            
          case 'absolute' % Reserved for future versions
                          %data = [0;data];
        end
        %Data.data = logical(data);                        
        Data.data = data;
    end
end            
Data.type = type;
end

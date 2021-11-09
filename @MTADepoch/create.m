function Data = create(Data,Session,funcHandle,varargin)            
% function Data = create(Data,Session,funcHandle,varargin)            
% 
% Calls one of the heuristic fuctions for state segmentation
% 
% Note: This fuction is not ready for general use
%
[data,l,k] = feval(funcHandle,Session);
[label,key] = DefaultArgs(varargin,{l,k});
Data.path = Session.spath;
Data.filename = Session.filebase;
Data.data = data;
Data.sampleRate = data.sampleRate;
Data.label = label;
Data.key = key;
Data.sync = Session.sync.copy();
Data.origin = Session.sync.data(1);
Data.update_filename(Session);

end

 
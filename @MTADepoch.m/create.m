function Data = create(Data,Session,funcHandle,varargin)            
% function Data = create(Data,Session,funcHandle,varargin)            
% 
% Calls one of the heuristic fuctions for state segmentation
% 
% Note: This fuction is not ready for general use
%
[data,sampleRate,l,k] = feval(funcHandle,Session);
[label,key] = DefaultArgs(varargin,{l,k});
Data.path = Session.spath;
Data.filename = Session.filebase;
Data.data = data;
Data.sampleRate = sampleRate;
Data.label = label;
Data.key = key;
end

 
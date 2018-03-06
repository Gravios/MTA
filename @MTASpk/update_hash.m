function Data = update_hash(Data,varargin)
% function Data = update_hash(Data,varargin)
% 
% 
%
% varargin:
%    any number of arguments to be used to generate hash
%    

if isempty(varargin),
    modificationHash = [];
end

Data.hash = DataHash(varargin);
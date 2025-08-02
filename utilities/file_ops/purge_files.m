function purge_files(Trials,pattern)
%function purge_files(Trials,fileExt)
% remove files from project folder of listed trials which contain 
% the regular expression pattern
if isa(Trials,'MTASession'), Trials = {Trials}; end
cf(@(t)   MTATrial.validate(t), Trials);
cf(@(t,pat) delete(fullfile(t.spath,list_files(t.name,pat,'string'))),...
   Trials,repmat({pattern},size(Trials)));
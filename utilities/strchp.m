function strout = strchp(str,varargin)
[nchar,pat] = DefaultArgs(varargin,{20,'\s'});
strout = {};
if ~iscell(str)&&ischar(str),
    strout{1} = str;
else
    strout = str;
end


strl = cellfun(@numel,strout);
if all(strl<nchar),
    return,
else
    strlind = find(strl>nchar,1,'first');
    pinds = regexp(strout{strlind},pat);
    strpiece1 = {strout{strlind}(1:pinds(find(pinds<nchar,1,'last')))};
    strpiece2 = {strout{strlind}(pinds(find(pinds<nchar,1,'last'))+1:end)};
    strout = strchp([strout(1:strlind-1),strpiece1,strpiece2,strout(strlind+1:end)],nchar,pat);
end

   



    
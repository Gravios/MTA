function out = xmlget(z,varargin)
%xmlget(structure,'KinematicModel','MarkerSet','Markers',{'Marker','Attributes','NAME','RGB'})
    
zpart = z;
if ~strcmp(varargin(1),z.Name)
    error(['Top level: ' varargin(1) ' does not exist'])
end

out = {};
for i = 2:length(varargin),
    for j = 1:length(zpart),
        if  iscell(varargin{i}),
            for k = 3:length(varargin{i}),
                switch varargin{i}{2}
                  case 'Attributes'
                    out{end+1} = getChildAttribute(zpart.Children,varargin{i}{1},varargin{i}{k});
                end
            end
        else
            try,
                zpart = getChild(zpart.Children,varargin{i});
            catch
                error('child does not exist');
            end
        end
    end
end

if isempty(out)
   out = zpart;
end

end
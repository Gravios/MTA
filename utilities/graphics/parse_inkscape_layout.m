function layout = parse_inkscape_layout(filepath,elementBaseName)

data = textscan( fileread(filepath), ...
                 '%s%f%f%f%f'   ,...
                 'Delimiter'    , ',',   ...
                 'CollectOutput', true  );

layout = struct('body',    struct('width',[],'height',[]),...
                'subplots',[]);
for element = 1:numel(data{1})
    if ~isempty(regexp(data{1}{element}, elementBaseName))
        layout.subplots.(data{1}{element}) = struct('x',      data{2}(element,1),...
                                                    'y',      data{2}(element,2),...
                                                    'width',  data{2}(element,3),...
                                                    'height', data{2}(element,4));
        
    elseif ~isempty(regexp(data{1}{element}, '^body$'))
        layout.body.width = data{2}(element,3);
        layout.body.height = data{2}(element,4);
        
    end
    
end

    
    
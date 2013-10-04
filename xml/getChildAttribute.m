function chattrib = getChildAttribute(children,name,attribute)
    chattrib ={};
    for i = 1:length(children)
        if strcmp(children(i).Name,name)
            for j = 1:length(children(i).Attributes)
                if strcmp(children(i).Attributes(j).Name,attribute)
                    chattrib{end+1} = children(i).Attributes(j);
                end
            end
        end
    end
end

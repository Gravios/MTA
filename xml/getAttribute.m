function attrib = getAttribute(a,name)
    for i = 1:length(a)
        if strcmp(a(i).Name,name),
            attrib = a(i);
        end
    end
end

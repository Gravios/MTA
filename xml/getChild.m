function child = getChild(c,name)
    for i = 1:length(c)
        if strcmp(c(i).Name,name),
            child = c(i);
        end
    end
end

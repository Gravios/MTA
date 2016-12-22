function indicies = select_ranges(figureHandle)
fprintf(['\n***************************************\n',...
           '* Select each range:                  *\n',...
           '*   Enter:   Activate selection tool  *\n',...
           '*   Escape:  Quit                     *\n',...
           '***************************************\n\n']); 
indicies=[];
ylm = ylim;
bflag = false;
while ~strcmp(get(figureHandle,'CurrentCharacter'),char(27));
    set(figureHandle,'CurrentCharacter',char(32));        
    while ~strcmp(get(figureHandle,'CurrentCharacter'),char(13))
        waitforbuttonpress
        if strcmp(get(figureHandle,'CurrentCharacter'),char(27))
            bflag = true;
            break
        end
    end
    if bflag, break,end
    rect = getrect(figureHandle);
    patch([rect(1),rect(1),rect(1)+rect(3),rect(1)+rect(3)],...
          [ylm(1),ylm(2),ylm(2),ylm(1)],...
          [-1,-1,-1,-1],[.9,.9,.9]);
    indicies = cat(1,indicies,[round(rect(1)):round(rect(1)+rect(3))]');
end

indicies = unique(indicies);

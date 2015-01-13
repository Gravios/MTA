function index = figure_controls(hfig,index,varargin)
%function index = figure_controls(hfig,index,varargin)
%[indmap,autoIncr] = DefaultArgs(varargin,{[],false});

[indmap,autoIncr,figname] = DefaultArgs(varargin,{[],false,[]});


if autoIncr,
    whatkey = 'n';
    if ~isempty(indmap),
        if index==indmap(end),index=-1;
            return;
        end
    end
else
    B = waitforbuttonpress;
    whatkey = get(hfig,'CurrentCharacter');
    if ~B,return,end
end

switch double(whatkey)
  case double('i')
    ii = index;
    index = input('Enter index #: ');
    if ~isempty(indmap)
        if isempty(find(indmap==index))
            warning(['Index: ' num2str(index) ', does not exist'])
            index = ii;
        end
    end


  case double('n')
    if ~isempty(indmap)
        ii = find(indmap==index);
        if (ii+1)>numel(indmap),
            index=indmap(1);
        else
            index = indmap(ii+1);
        end
    else
        index = index+1;
    end
    
  case double('p')
    if ~isempty(indmap)
        ii = find(indmap==index);
        if (ii-1)<1,
            index=indmap(end);
        else
             index = indmap(ii-1);
        end
    else
        index=index-1;
    end
    
  case double('s')
    saveas(hfig,figname,'png');
    
  case double('q')
    index = -1;
end

function index = figure_controls(hfig,index,varargin)
%function index = figure_controls(hfig,index,varargin)
%[indmap,autoIncr] = DefaultArgs(varargin,{[],false});

[indmap,autoIncr,figname,flags] = DefaultArgs(varargin,{[],false,[],''});


if autoIncr,                              % ITERATE automatically over indmap
    currentChar = 'n';
    if ~isempty(indmap),
        if index==indmap(end),index=-1;
            return;
        end
    end
else                                      % ITERATE manually over indmap
    B = waitforbuttonpress();
    currentChar = hfig.CurrentCharacter;  % GET the last key pressed
    if ~B,return,end
end

switch double(currentChar)                % DO some action
  case double('i')
    ii = index;
    index = input('Enter index #: ');
    if ~isempty(indmap)
        if isempty(find(indmap==index))
            warning(['Index: ' num2str(index) ', does not exist'])
            index = ii;
        end
    end

  case double('n')                        % GET next index
    if ~isempty(indmap)                   % GET next index within indmap
        ii = find(indmap==index);
        if (ii+1)>numel(indmap),          % RETURN to begining of indmap            
            index=indmap(1);
        else                              % GET next index
            index = indmap(ii+1);
        end
    else                                  % GET next index        
        index = index+1;
    end
    
  case double('p')                        % GET previous index    
    if ~isempty(indmap)                   % GET previous index within indmap        
        ii = find(indmap==index);
        if (ii-1)<1,                      % MOVE to the end of indmap
            index=indmap(end);
        else                              % GET previous index                
             index = indmap(ii-1);
        end
    else                                  % GET previous index            
        index=index-1;
    end
    
  case double('s')                        % PRINT figure
    saveas(hfig,figname,'png');
    
  case double('q')                        % QUIT figure    
    index = -1;
end


for f = flags,
    switch f
      case '-'
      case 'v'
        sprintf('Figure Id: %d\nIndex: %i\n',[hfig,index])
      otherwise
        warning(['flag: ' f ' does not exist in figure_controls.m, ' ...
                            'will ignore.']);
    end
end

    

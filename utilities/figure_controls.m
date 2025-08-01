function index = figure_controls(hfig,index,varargin)
%function index = figure_controls(hfig,index,varargin)
%[map,autoIncr] = DefaultArgs(varargin,{[],false});

[map,autoIncr,figname,flags,groups] = DefaultArgs(varargin,{[],false,[],'',[]});

map_length        = numel( map );
map_is_not_empty  = logical( map_length );
current_map_index = find( index == map );

if autoIncr,                              % ITERATE automatically over map
    key = 'n';
    if map_is_not_empty,
        if index==map(end)
            index = -1;
            return;
        end
    end
else                                      % ITERATE manually over map
    B = waitforbuttonpress();
    key = hfig.CurrentCharacter;  % GET the last key pressed
    if ~B,return,end
end

switch double(key) 
  case double('i')
    old_map_index = map(find( index == map ));
    new_index = input('Enter index #: ');
    new_map_index = find( new_index == map );
    
    if map_is_not_empty && new_map_index
        index = map( new_map_index);
    else
        disp(['[INFO] index: ' num2str(new_index) ' does not exist.']);
    end
    
% MOVE to next index        
  case double('n')
    if map_is_not_empty 
        new_index = find( map == index) + 1;
        if new_index > map_length % LOOP back to the begining
            index = map(1);
        else 
            index = map( new_index );
        end
    else
        index = index + 1;
    end
    
% MOVE to previous index
  case double('p')                        
    if map_is_not_empty
        previous_index = current_map_index - 1;
        if previous_index < 1 % LOOP to the end of the map
            index = map( end );
        else
            index = map( previous_index );
        end
    else                         
        index = index - 1;
    end

% APPEND index to a group
  case double('m')
    disp(sprintf('[INFO] append to group: '))
    B = waitforbuttonpress();
    gkey = hfig.CurrentCharacter;
    key_exists_in_group = ismember( gkey, groups.keys);
    if isempty(B)
        hfig.CurrentCharacter = 'm';
    else
        if key_exists_in_group
            disp(sprintf('%s\n\n', gkey));
            groups(gkey) =                                                   ...
                setfield(                                                   ...
                    groups(gkey),                                            ...
                    'data',                                                 ...
                    [groups(gkey).data, index]                               ...
                );
        else
            hfig.CurrentCharacter = 'm';
            disp(sprintf('[INFO] append to group: '))
        end
    end

% PRINT figure
  case double('s')                        
    saveas( hfig, [figname '-' num2str(index) '.png'], 'png');

% QUIT figure    
  case double('q')
    index = -1;

  case double('h')
    disp(sprintf('\n\n'));
    disp(sprintf('+-- figure_controlls - Help Section ----------------------------------+'));
    disp(sprintf('|                                                                     |'));
    disp(sprintf('| Keyboard Shortcuts:                                                 |'));
    disp(sprintf('|                                                                     |'));
    disp(sprintf('|   Navigation:                                                       |'));
    disp(sprintf('|     q : quit figure_controlls                                       |'));
    disp(sprintf('|     n : move to next index                                          |'));
    disp(sprintf('|     p : move to previous index                                      |'));
    disp(sprintf('|     i : select an index                                             |'));
    disp(sprintf('|                                                                     |'));
    disp(sprintf('|   Groups Labeling:                                                  |'));
    gkeys = groups.keys;
    for gid = 1:double(groups.Count)
        gstring = ['|    ', gkeys{gid}, ' : ',groups(gkeys{gid}).label];
        if length(gstring) < 73
            gstring = [ gstring, repmat(' ',[1, 70 - length(gstring)]), '|' ];
        end
        disp(gstring);
    end
    disp(sprintf('|                                                                     |'));
    disp(sprintf('+---------------------------------------------------------------------+\n'));

end

for f = flags,
    switch f
      case '-'
      case 'v'
        disp(sprintf('Figure Id: %d \nIndex: %i\n', [hfig.Number, index]))
      otherwise
        sprintf('[WARING] MTA:utilities:figure_controls - Invalid flag \n')
    end
end

    

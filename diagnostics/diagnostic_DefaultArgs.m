function diagnostic_DefaultArgs(varargin)

    fprintf('Diagnostic for DefaultArgs\n')    
    
    fprintf('Test for empty varargin\n')
    varargin = {};
    [tv_int,tv_str,tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}})
    varargin = {2,'bob',{'top',2}}
    fprintf('check tv_int: %i\n',tv_int   == varargin{1});
    fprintf('check tv_str: %i\n',tv_str   == varargin{2});
    fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1});
    fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2});
        
    fprintf('\n\n')
    
    fprintf('Test for named fields in order\n')
    varargin = {'tv_int',2,'tv_str','bob','tv_cell',{'top',2}}
    [tv_int,tv_str,tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}});
    varargin = {2,'bob',{'top',2}}
    fprintf('check tv_int: %i\n',tv_int   == varargin{1})
    fprintf('check tv_str: %i\n',tv_str   == varargin{2})
    fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1})
    fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2})
    
    fprintf('\n\n')
    
    fprintf('Test for named fields random order [2,1,3]\n')
    varargin = {'tv_int',2,'tv_str','bob','tv_cell',{'top',2}}
    [tv_int,tv_str,tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}});
    varargin = {2,'bob',{'top',2}}
    fprintf('check tv_int: %i\n',tv_int   == varargin{1})
    fprintf('check tv_str: %i\n',tv_str   == varargin{2})
    fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1})
    fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2})

    fprintf('\n\n')
    
    fprintf('Test for named fields random order [3,2,1]\n')
    varargin = {'tv_cell',{'top',2},'tv_str','bob','tv_int',2}
    [tv_int,tv_str,tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}})
    varargin = {2,'bob',{'top',2}}
    fprintf('check tv_int: %i\n',tv_int   == varargin{1})
    fprintf('check tv_str: %i\n',tv_str   == varargin{2})
    fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1})
    fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2})

    fprintf('\n\n')
    
    fprintf(['Test for 1 normal order and 2 randomly ordered\n'...
             'named fields random order [1,3,2]\n']);
    varargin = {2,'tv_cell',{'top',2},'tv_str','bob'}
    [tv_int,tv_str,tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}});
    varargin = {2,'bob',{'top',2}}
    fprintf('check tv_int: %i\n',tv_int   == varargin{1});
    fprintf('check tv_str: %i\n',tv_str   == varargin{2});
    fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1});
    fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2});
    
    fprintf('\n\n')
    
    fprintf(['Elipses with multiple assignment\n'...
             'named fields random order [1,3,2]\n']);
    varargin = {2,'tv_cell',{'top',2},'tv_str','robert'}
    [tv_int,...
     tv_str,... test comment
     ...
     tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}});
    varargin = {2,'robert',{'top',2}}
    fprintf('check tv_int: %i\n',tv_int   == varargin{1});
    fprintf('check tv_str: %i\n',tv_str   == varargin{2});
    fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1});
    fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2});
    
    fprintf('\n\n')
    
    fprintf('Test single variabel assignment with empty varargin\n');
    varargin = {};
    tv_str = DefaultArgs(varargin,{'bob'});
    varargin = {'bob'}
    fprintf('check tv_str: %i\n',tv_str   == varargin{1});
    
    fprintf('\n\n')
    
    fprintf('Test single variabel assignment with empty varargin and elpsis\n');
    varargin = {};
    tv_str =...
        DefaultArgs(varargin,{'bob'});
    varargin = {'bob'}
    fprintf('check tv_str: %i\n',tv_str   == varargin{1});

end
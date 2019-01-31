function diagnostic_DefaultArgs()

    fprintf('Diagnostic for DefaultArgs\n')    
    
    fprintf('Test deep function defargs override with global var AP \n');        
    global AP;
    AP.MTADiagnostic_helper_deep_functions.ron = 'can';    
    MTADiagnostic_helper_deep_functions({});
    fprintf('complete\n\n');        

    fprintf('Test shallow defargs override with global var AP \n');        
    global AP;
    AP.diagnostic_DefaultArgs.ron = 'bot';    
    fprintf('complete\n\n');        
    defArgs = struct('bob','cool',...
                     'ron','not',...
                     'peter','cant');    
    varargin = {};
    [bob,ron,peter] = DefaultArgs(varargin,defArgs,'--struct');
    fprintf('Control var, bob( default = cool ): %s\n',bob);
    fprintf('OVERRIDE var, ron( default = not ): %s\n',ron);    
    fprintf('Control var, peter( default = cant ): %s\n',peter);
    fprintf('complete\n\n');            

    
    
    
    
    fprintf('Test no struct shallow defargs override with global var AP \n');        
    global AP;
    AP.diagnostic_DefaultArgs.ron = 'blast';        
    fprintf('complete\n\n');        
    varargin = {};
    [bob,ron,peter] = DefaultArgs(varargin,{'cool','not','cant'});
    fprintf('Control var, bob( default = cool ): %s\n',bob);
    fprintf('OVERRIDE var, ron( default = not ): %s\n',ron);    
    fprintf('Control var, peter( default = cant ): %s\n',peter);
    fprintf('complete\n\n');            
    
    
    
% $$$     
% $$$ 
% $$$     fprintf('Test for empty varargin\n')
% $$$     varargin = {};
% $$$     [tv_int,tv_str,tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}})
% $$$     varargin = {2,'bob',{'top',2}}
% $$$     fprintf('check tv_int: %i\n',tv_int   == varargin{1});
% $$$     fprintf('check tv_str: %i\n',tv_str   == varargin{2});
% $$$     fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1});
% $$$     fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2});
% $$$ 
% $$$     fprintf('\n\n')
% $$$ 
% $$$     fprintf('Test for named fields in order\n')
% $$$     varargin = {'tv_int',2,'tv_str','bob','tv_cell',{'top',2}}
% $$$     [tv_int,tv_str,tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}});
% $$$     varargin = {2,'bob',{'top',2}}
% $$$     fprintf('check tv_int: %i\n',tv_int   == varargin{1})
% $$$     fprintf('check tv_str: %i\n',tv_str   == varargin{2})
% $$$     fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1})
% $$$     fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2})
% $$$ 
% $$$     fprintf('\n\n')
% $$$ 
% $$$     fprintf('Test for named fields random order [2,1,3]\n')
% $$$     varargin = {'tv_int',2,'tv_str','bob','tv_cell',{'top',2}}
% $$$     [tv_int,tv_str,tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}});
% $$$     varargin = {2,'bob',{'top',2}}
% $$$     fprintf('check tv_int: %i\n',tv_int   == varargin{1})
% $$$     fprintf('check tv_str: %i\n',tv_str   == varargin{2})
% $$$     fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1})
% $$$     fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2})
% $$$ 
% $$$     fprintf('\n\n')
% $$$ 
% $$$     fprintf('Test for named fields random order [3,2,1]\n')
% $$$     varargin = {'tv_cell',{'top',2},'tv_str','bob','tv_int',2}
% $$$     [tv_int,tv_str,tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}})
% $$$     varargin = {2,'bob',{'top',2}}
% $$$     fprintf('check tv_int: %i\n',tv_int   == varargin{1})
% $$$     fprintf('check tv_str: %i\n',tv_str   == varargin{2})
% $$$     fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1})
% $$$     fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2})
% $$$ 
% $$$     fprintf('\n\n')
% $$$ 
% $$$     fprintf(['Test for 1 normal order and 2 randomly ordered\n'...
% $$$              'named fields random order [1,3,2]\n']);
% $$$     varargin = {2,'tv_cell',{'top',2},'tv_str','bob'}
% $$$     [tv_int,tv_str,tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}});
% $$$     varargin = {2,'bob',{'top',2}}
% $$$     fprintf('check tv_int: %i\n',tv_int   == varargin{1});
% $$$     fprintf('check tv_str: %i\n',tv_str   == varargin{2});
% $$$     fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1});
% $$$     fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2});
% $$$ 
% $$$     fprintf('\n\n')
% $$$ 
% $$$     fprintf(['Elipses with multiple assignment\n'...
% $$$              'named fields random order [1,3,2]\n']);
% $$$     varargin = {2,'tv_cell',{'top',2},'tv_str','robert'}
% $$$     [tv_int,...
% $$$      tv_str,... test comment
% $$$      ...
% $$$      tv_cell] = DefaultArgs(varargin,{2,'bob',{'top',2}});
% $$$     varargin = {2,'robert',{'top',2}}
% $$$     fprintf('check tv_int: %i\n',tv_int   == varargin{1});
% $$$     fprintf('check tv_str: %i\n',tv_str   == varargin{2});
% $$$     fprintf('check tv_cell{1}: %i\n',tv_cell{1}   == varargin{3}{1});
% $$$     fprintf('check tv_cell{2}: %i\n',tv_cell{2}   == varargin{3}{2});
% $$$ 
% $$$     fprintf('\n\n')
% $$$ 
% $$$     fprintf('Test single variabel assignment with empty varargin\n');
% $$$     varargin = {};
% $$$     tv_str = DefaultArgs(varargin,{'bob'});
% $$$     varargin = {'bob'}
% $$$     fprintf('check tv_str: %i\n',tv_str   == varargin{1});
% $$$ 
% $$$     fprintf('\n\n')
% $$$ 
% $$$     fprintf('Test single variabel assignment with empty varargin and elpsis\n');
% $$$     varargin = {};
% $$$     tv_str =...
% $$$         DefaultArgs(varargin,{'bob'});
% $$$     varargin = {'bob'}
% $$$     fprintf('check tv_str: %i\n',tv_str   == varargin{1});
% $$$ 



function MTADiagnostic_helper_deep_functions(varargin)

    defArgs = struct('bob','cool',...
                     'ron','not',...
                     'peter','cant');
    
    [bob,ron,peter] = DefaultArgs(varargin,defArgs,'--struct');

    fprintf('Control var, bob( default = cool ): %s\n',bob);
    fprintf('OVERRIDE var, ron( default = not ): %s\n',ron);    
    fprintf('Control var, peter( default = cant ): %s\n',peter);
    
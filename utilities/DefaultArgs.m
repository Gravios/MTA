function varargout = DefaultArgs(Args, DefArgs, varargin)
%% DefaultArgs(Args, DefArgs, FuncPath, NLineExcerpt ,NLineBackSeek)
% auxillary function to replace argument check in the beginning and def. 
% args assigmentsets the absent or empty values of the Args (cell array, 
% usually varargin)to their default values from the cell array DefArgs. 
% Output should contain the actuall names of arguments that you use in 
% the function
%
% e.g. : in function MyFunction(Arg1,Arg2,... , varargin)
% calling [SampleRate, BinSize] = DefaultArgs(varargin, {20000, 20});
% will assign the defualt values to SampleRate and BinSize arguments if they
% are empty or absent in the varargin cell list 
% (not passed to a function or passed empty)
%
% Update: The output variable names are now parsed directly from the
% function which calls DefaultArgs, via dbstack. All variable names are 
% then compaired to all string variables in Args. Matching string inputs
% will have the subsequent variable set as its value.
%
% Modes:
%    Static
%  
%
% NOTE: Variable names are matched case insensitive
% 
% NOTE: Always use comma separated multiple assignment syntax
%
% NOTE: The varargin in variables must be in order up to the first matching
%       string matching an output variables name.
%
% example (implementation):
%
%   Function Path is not specified, uses dbstack, to determin which
%      ExampleFunction(somearguments , varargin)
%      [SampleRate,BinSize] = DefaultArgs(varargin,{2000,20});
%
%   Specified function path, 
%      ExampleFunction(somearguments , varargin)
%      [SampleRate,BinSize] = DefaultArgs(varargin,{2000,20},which('ExampleFunction'));
%
% example (Usage): ExampleFunction(somearguments , SampleRate, BinSize, Flag, Overwrite)
%
%      ExampleFunction(somearguments, {'SampleRate',30000,'BinSize',30});
% 
%      ExampleFunction(somearguments, {30000,'BinSize',30,'Flag','-comp'});
%
%      ExampleFunction(somearguments,30000,30,'overwrite',true,'Flag','-comp'});
%
%      ExampleFunction(somearguments,'Overwrite',1)
%
% See also struct2varargin

%% Set DefaultArgs' default args


% IF varargin is empty, setup path to the m-file calling, hereafter known as CMF, DefaultsArgs.
if numel(varargin)>0,    

    switch varargin{1}
      
      case '--struct'
        
% ASSIGN default aguments fields in order to varargout cell array
        varargout = {};
        defaultArgFields = fields(DefArgs);
        assert( nargout == numel(defaultArgFields), ...
                'MTA:utilities:set_default_args:NArginDoesNotEqualNArgout',...
                'Count of outputs does not match the number of default arguments');
        for field = defaultArgFields', varargout{end+1} = DefArgs.(field{1});  end
        
% QUERY default argument variables from input struct
        numInputArgs = numel(Args);
        index = 1;
        nameValueMode = false;
        while index <= numInputArgs
            
            if ischar(Args{index}) && isfield(DefArgs,Args{index}),
% ASSIGN next postion in Args to varargout based on field tag in previous position in Args
                outputIndex = ~cell2mat(cf(@isempty,regexp(Args{index},cf(@(x) ['^',x,'$'],defaultArgFields'))));
                varargout(outputIndex) = Args(index+1);
                index = index+2;
                nameValueMode = true;
                
            elseif nameValueMode, error('MTA:utilities:set_default_args:MissingNameValuePair');
            
            else,
% ASSIGN next position in Args to varargout assuming correct order
                if ~isempty(Args{index})
                    varargout(index) = Args(index);
                end
                index = index+1;
            end
        end
        return;

      otherwise
        FuncPath = varargin{1};        
    end
else
    FuncPath = dbstack('-completenames'); 
end

% THROW error if output count doesn't equal the number of default args
assert( nargout==numel(DefArgs), ...
        'MTA:utilities:set_default_args:NArginDoesNotEqualNArgout',...
        'Count of outputs does not match the number of default arguments');


% assume generic array or logical
% 1. single element array is a flag to not attempt parameter value paring 
% 2. Character array with a total number of elements greater than 1
%    will be treated as the path for the CMF.
if numel(FuncPath)>1,
    % The total number of lines read from the CMF.
    NLineExcerpt = 15;
    if numel(varargin)>1,
        if ~isempty(varargin{2}),
            NLineExcerpt  = varargin{2};
        end
    end

% DEFINED line position for reading the CMF relative to the of the DefaultArgs call in the CMF.
    NLineBackSeek = 4;
    if numel(varargin)>2,
        if ~isempty(varargin{3}),
            NLineBackSeek = varargin{3};
        end
    end

    if ~iscell(DefArgs)
% WRAP DefArgs in cell
        DefArgs = {DefArgs};
    end

% Don't remember what this does 
    skiphead = 0;
    if isstruct(FuncPath),
        skiphead = FuncPath(2).line-NLineBackSeek;
        if skiphead<0, skiphead=0; end
        FuncPath = FuncPath(2).file;
    end

% OPEN matfile of function, copy code where DefaultArgs is located (default is 15 lines)
    fid = fopen(FuncPath);
    FSC = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '\b','HeaderLines',skiphead);
    fclose(fid);
    FSC = FSC{1};
    if numel(FSC)>NLineExcerpt,FSC = FSC(1:NLineExcerpt);end


% CLEAN function code excerpt 
    for i = 1:numel(FSC),
% REMOVE white spaces
        FSC{i}(ismember(FSC{i},' ')) = [];
        if  isempty(FSC{i}),continue,end
% FLAG lines without code for deletion (i.e. comments)
        if strcmp(FSC{i}(1),'%'),
            FSC{i} = [];
        end
    end
% DELETE flagged lines
    cind = find(cellfun(@isempty,FSC));
    if ~isempty(cind), 
        FSC(cind)=[];
    end

% PARSE variable names from excerpt -------------------------------------------------------
% FIND the line where DefaultArgs exists
    dfa_index = find(~cellfun(@isempty,regexpi(FSC,'defaultargs','once')),1);

    outind = [];
    tind = dfa_index;
    tdfas = '';
    while isempty(outind),
        tdfas = strcat(FSC{tind},tdfas);

        % Seach for the begining of variable assignments
        outind = find(ismember(tdfas(1),'[')==1,1,'first');    

        % Search backwards for ellipses 
        if tind-1==0,
            prvln = [];
        else
            prvln = strfind(FSC{tind-1},'...');        
        end
        
        if ~isempty(prvln)
            tind = tind-1;
            FSC{tind}(prvln:numel(FSC{tind}))=[];
        else
            if isempty(outind)
                outind = -1;
            end
        end
    end

    tdfas = tdfas(1:strfind(tdfas,'='));
    if outind~=-1,
        tdfas = tdfas(strfind(tdfas,'[')+1:strfind(tdfas,']')-1);
    else
        tdfas = tdfas(1:numel(tdfas)-1);
    end

    ParsedArgs = str2cell(tdfas,' ,');
    RPargs = strcat('(^',strjoin(ParsedArgs(:)','$)|(^'),'$)');
    strArgsInd = cellfun(@isstr,Args);


    nDefArgs = length(DefArgs);
    if (nargout~=nDefArgs)
        error('number of defaults is different from assigned');
    end

    NewMAI = -1;
    if numel(strArgsInd)>0,
        MatchedArgs = cell(size(DefArgs));
        MatchedArgsInd = false(size(strArgsInd));
        MatchedArgsInd(strArgsInd) = ~cellfun(@isempty,regexpi(Args(strArgsInd),RPargs));
        MatchedArgsInd = find(MatchedArgsInd);
        
        if ~isempty(MatchedArgsInd)        
            mai = MatchedArgsInd;
            while ~isempty(mai)
                ArgInd = find(~cellfun(@isempty,regexpi(ParsedArgs,{['^',Args{mai(1)},'$']})));
                NewMAI(numel(NewMAI)+1) = ArgInd;
                MatchedArgs{ArgInd} = Args{mai(1)+1};
                mai(1) = [];
            end
            if MatchedArgsInd(1)>1
                tArgs = Args;
                Args = cell(size(ParsedArgs));
                Args(1:MatchedArgsInd(1)-1) = tArgs(1:MatchedArgsInd(1)-1);
            else
                Args = cell(size(ParsedArgs));
            end
        end
    end

    %% Ensure the Args is a cell with at least one value
    if isempty(Args)
        Args = cell(nDefArgs,1);
    elseif numel(Args) ~= nDefArgs,
        tArgs = Args;
        Args = cell(nDefArgs,1);
        Args(1:numel(tArgs)) = tArgs;
    end

    varargout = cell(nDefArgs,1);

    for i=1:nDefArgs
        if isempty(Args{i})&& ~ismember(i,NewMAI)
            varargout(i) = DefArgs(i);
        elseif ~ismember(i,NewMAI)
            varargout(i) = Args(i);
        else
            varargout(i) = MatchedArgs(i);
        end
    end
    
else
    if isempty(Args),
        Args ={[]}; 
    end
    if ~iscell(DefArgs),
        DefArgs = {DefArgs};
    end
    nDefArgs = length(DefArgs);
    numInputArgs = length(Args);
    if (nargout~=nDefArgs), ...
                 error('number of defaults is different from assigned'); 
    end
    
    for i = 1:nDefArgs,
        if(i>numInputArgs || isempty(Args{i}))
            varargout(i) = {DefArgs{i}};
        else
            varargout(i) = {Args{i}};
        end
    end
end

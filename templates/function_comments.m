function function_comments(varargin)
% function function_comments(varargin)
% 
% varargin:
%    var1 - numericArray: data provided for transformation
%    var2 - numeric:      transformation parameter
%


% DEFARGS ------------------------------------------------------------------------------------------

% EXAMPLE 
defargs = struct('var1',                          'data provided for transform',                 ...
                 'var2',                          'some parameter',                              ...
                 'var3',                          'another parameter'                            ...
);
[var1,var2,var3] = DefaultArgs(varargin,defargs,'--struct');

% TEMPLATE  
defargs = struct(...
);
[] = DefaultArgs(varargin,defargs,'--struct');

%---------------------------------------------------------------------------------------------------


% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag

% EXAMPLE 
if isempty(tag),
    tag = DataHash(struct(var2,var3));
end

% TEMPLATE  
if isempty(tag),
    tag = DataHash(struct());
end
%---------------------------------------------------------------------------------------------------


% ASSERTIONS ---------------------------------------------------------------------------------------
assert(regexp(mode,'^COMPUTE$'),'MTA:transformations:decompose_xy_motion_wrt:UnknownMode');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------



% UPDATE hash property of Data object
Data.update_hash(DataHash(struct('order',order,'numApplications',numApplications)));



% END MAIN -----------------------------------------------------------------------------------------


% AUX METHODS --------------------------------------------------------------------------------------
% END AUX METHODS ----------------------------------------------------------------------------------



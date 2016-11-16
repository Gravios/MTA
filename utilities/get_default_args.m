function defargs = get_default_args(varargin)
%function defargs = get_default_args(anaysisTag,funcName,outputType)
%
% loads default arguments for a function from an m-file define by
% the user.
%
% the m-file should be stored in the project folder with format:
%   ['get_default_args_',anaysisTag,'.m']
%
% varargin:
%   analysisTag:    string, some short tag for a given analysis or sub analysis
%   funcName:       string, any function enumerated in the target m-file
%   outputType:     string, valid values - 'struct' and 'cell'
global MTA_CURRENT_PROJECT
global MTA_PROJECT_PATH

[anaysisTag,funcName,outputType] = DefaultArgs(varargin,{MTA_CURRENT_PROJECT,'','struct'},true);

funHandle = function_handle(fullfile(MTA_PROJECT_PATH,...
                                     ['get_default_args_',anaysisTag,'.m'])...
                            );

defargs = funHandle(funcName,outputType);
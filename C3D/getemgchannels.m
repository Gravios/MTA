function EMG = getemgchannels(itf, chanstart, chanend, type, index1, index2)
% GETEMGCHANNELS - returns emg data in raw or processed forms.
% 
%   USAGE:  EMG = getemgchannels(itf, chanstart, chanend, type*, index1*, index2*)
%           * = not a necessary input
%   INPUTS:
%   itf         = variable name used for COM object
%   chanstart   = first emg channel
%   chanend     = last emg channel
%   type        = 0 or no input for raw, 1 for FWR, 2 for LE (6Hz single-pass BW)
%   index1      = start frame index, all frames if not used as an argument
%   index2      = end frame index, all frames if not used as an argument
%   OUTPUTS:
%   EMG         = stucture with channel fields

%   C3D directory contains C3DServer activation and wrapper Matlab functions.
%   This function written by:
%   Matthew R. Walker, MSc. <matthewwalker_1@hotmail.com>
%   Michael J. Rainbow, BS. <Michael_Rainbow@brown.edu>
%   Motion Analysis Lab, Shriners Hospitals for Children, Erie, PA, USA
%   Questions and/or comments are most welcome.  
%   Last Updated: April 21, 2006
%   Created in: MATLAB Version 7.0.1.24704 (R14) Service Pack 1
%               O/S: MS Windows XP Version 5.1 (Build 2600: Service Pack 2)
%   
%   Please retain the author names, and give acknowledgement where necessary.  
%   DISCLAIMER: The use of these functions is at your own risk.  
%   The authors do not assume any responsibility related to the use 
%   of this code, and do not guarantee its correctness. 


% CHECK INPUTS
if nargin == 4, 
    index1 = itf.GetVideoFrame(0); % frame start
    index2 = itf.GetVideoFrame(1); % frame end
elseif nargin == 3,
    type = 0; % raw EMG
    index1 = itf.GetVideoFrame(0); % frame start
    index2 = itf.GetVideoFrame(1); % frame end
elseif nargin == 1 | nargin == 2 | nargin == 5,
    disp('Error: wrong number of inputs.');
    help getemgchannels; return
end
% CHECK ANALOG FORMAT
 fIndex = itf.GetParameterIndex('ANALOG', 'FORMAT');
 if fIndex == -1,
     disp('Sorry, this function does not support your analog data format at this time.');
     return
 else,
    format = itf.GetParameterValue(fIndex, 0);
    if upper(format(1)) == 'S',
        disp('Error: only supports VICON unsigned analog data at this time');
    return
    else, 
        disp('C3D file contains unsigned analog data, retrieving data...');
    end
 end

% CHECK EMG LABELS
 nlabelIndex = itf.GetParameterIndex('ANALOG', 'LABELS');
 label1 = itf.GetParameterValue(nlabelIndex, chanstart-1);
if upper(label1(2)) ~= 'E',
    disp('Error: the emg channels are labelled incorrectly')
    help getemgchannels; return
end
 
% GET EMG DATA
 nIndex = itf.GetParameterIndex('ANALOG', 'LABELS');
 unitIndex = itf.GetParameterIndex('ANALOG', 'UNITS');
 rateIndex = itf.GetParameterIndex('ANALOG', 'RATE');
 srate = itf.GetParameterValue(rateIndex, 0);
for i = chanstart : chanend,
    channel_name = itf.GetParameterValue(nIndex, i-1);
    newstring = channel_name(1:min(findstr(channel_name, ' '))-1);
    if strmatch(newstring, [], 'exact'),
        newstring = channel_name;
    end
    EMG.(newstring) = ...
        itf.GetAnalogDataEx(i-1,index1,index2,'1',0,0,'0');
    EMG.(newstring) = cell2mat(EMG.(newstring));
    if type == 1, 
        EMG.(newstring) = abs(EMG.(newstring));
    elseif type == 2,
        EMG.(newstring) = singlebutt(srate, 6, (abs(EMG.(newstring))), 0);
    end
    EMG.units.(newstring) = itf.GetParameterValue(unitIndex, i-1);
end
%--------------------------------------------------------------------------
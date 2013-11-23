function truncated = analog2videoframes(itf, signals, index1, index2)
% ANALOG2VIDEOFRAMES - truncates input signals to video frame rate.
% 
%   USAGE:  truncated = analog2videoframes(itf, signals, index1*, index2*)
%           * = not a necessary input
%   INPUTS:
%   itf     = variable name used for COM object
%   signals = array to be truncated
%   index1  = video start frame index, all frames if not used as an argument
%   index2  = video end frame index, all frames if not used as an argument
%   OUTPUTS:
%   truncated  = truncated signals array
   
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


if nargin == 2, 
    index1 = itf.GetVideoFrame(0); % frame start
    index2 = itf.GetVideoFrame(1); % frame end
elseif nargin == 1 | nargin == 3, 
    disp('Error: incorrect number of inputs.');
    help analog2videoframes;
    return;
end
 mult = itf.GetAnalogVideoRatio;
 truncated = signals((1+((index1-1)*mult)):mult:(1+((index2-1)*mult)), :); 
%--------------------------------------------------------------------------
    
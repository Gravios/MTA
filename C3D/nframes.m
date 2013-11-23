function frames = nframes(itf)
% NFRAMES - returns integer value of the number of video frames
% contained in an open C3D file.
% 
% USAGE: frames = nframes(itf)
% INPUTS:
%   itf  = variable name used for COM object
% OUTPUTS:
%   frames = integer value of number of video frames in c3d file

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


 frames = itf.GetVideoFrame(1) - itf.GetVideoFrame(0) + 1;
%--------------------------------------------------------------------------
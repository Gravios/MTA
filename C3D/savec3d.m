function savec3d(itf, newfname, filetype)
% SAVEC3D - saves C3D file.
% 
% USAGE: savec3d(itf, newfname*, filetype*)
%        * = not a necessary input
% INPUTS:
%   itf      = variable name used for COM object
%   newfname = new filename to save, '' or no input to save to open file
%   filetype = Intel(1), DEC(2), SGI(3), existing type (-1)
%
%   Note: 'newfname' and 'filetype' default to the name of the currently 
%   open file and the existing filetype if only 'itf' is entered as input.
%   savec3d will handle 'newfname' with or without the '.c3d' extension.
%   'newfname' can be entered with directory path (to save in a directory 
%   different from your current Matlab path).

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


if nargin == 3,
    dotind = findstr(newfname,'.');
    if isempty(dotind) == 1, newfname = [newfname '.c3d']; end
elseif nargin == 2,
    filetype = -1;
    dotind = findstr(newfname,'.');
    if isempty(dotind) == 1, newfname = [newfname '.c3d']; end
elseif nargin == 1,
    newfname = ''; filetype = -1;
end

 pRet = itf.SaveFile(newfname, filetype); 
if pRet == 1 & nargin ~=1, 
    disp([newfname ' has been saved.']); 
elseif pRet == 1 & nargin == 1,
    disp('Existing file has been saved.');
end
%--------------------------------------------------------------------------
 function XYZPOS = get3dtargets(itf, residual, index1, index2)
% GET3DTARGETS - returns structure containing all X,Y,Z trajectory data and
% residuals if chosen.
% 
%   USAGE:  XYZPOS = get3dtargets(itf, residual*, index1*, index2*) 
%           * = not a necessary input
%   INPUTS:
%   itf        = variable name used for COM object
%   residual   = Return matrix with point residual in column 4.  
%                0 or no 3rd argument = false (returns nx3 with XYZ data only)
%                1 = true (returns nx4 with XYZ and residuals) 
%   index1     = start frame index, all frames if not used as an argument
%   index2     = end frame index, all frames if not used as an argument
%   OUTPUTS:
%   XYZPOS     = structure with target fields of X, Y, Z, and/or residual as columns
   
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
     

if nargin == 1, 
    residual = 0;
    index1 = itf.GetVideoFrame(0); % frame start
    index2 = itf.GetVideoFrame(1); % frame end
elseif nargin == 2, 
    index1 = itf.GetVideoFrame(0); 
    index2 = itf.GetVideoFrame(1); 
end
 nIndex = itf.GetParameterIndex('POINT', 'LABELS');
 nItems = itf.GetParameterLength(nIndex);
 unitIndex = itf.GetParameterIndex('POINT', 'UNITS');
    
for i = 1 : nItems,
    target_name = itf.GetParameterValue(nIndex, i-1);
    newstring = target_name(1:min(findstr(target_name, ' '))-1);
    if strmatch(newstring, [], 'exact'),
        newstring = target_name;
    end
    if findstr('-', newstring) >= 1, 
        slashind = findstr('-', newstring); 
        newstring = [newstring(1:slashind-1) newstring(slashind+1:end)];
    end
    if strcmpi(newstring(1), '*'), newstring = ['T' newstring(2:end)]; end
    XYZPOS.(newstring) = ...
            [itf.GetPointDataEx(i-1,0,index1,index2,'1'), ...
             itf.GetPointDataEx(i-1,1,index1,index2,'1'), ...
             itf.GetPointDataEx(i-1,2,index1,index2,'1')];
    RESIDS = itf.GetPointResidualEx(i-1,index1,index2);
    XYZPOS.(newstring) = cell2mat(XYZPOS.(newstring));
    RESIDS = cell2mat(RESIDS);
    residindex = find(RESIDS == -1);
    XYZPOS.(newstring)(residindex, :) = NaN;
    XYZPOS.units = itf.GetParameterValue(unitIndex, 0);
    if residual == 1,
    XYZPOS.(newstring) = [XYZPOS.(newstring), RESIDS];
    end 
end 
%--------------------------------------------------------------------------
function [XYZPOS,MARKERS] = getAll3dTargets(itf, index1, index2)
% GET3DTARGET - returns nx4 array containing X,Y,Z trajectory data and
% residual.
% 
%   USAGE:  XYZPOS = get3dtarget(itf, signalname, residual*, index1*, index2*) 
%           * = not a necessary input
%   INPUTS:
%   itf        = variable name used for COM object
%   signalname = string name of desired marker
%   residual   = Return matrix with point residual in column 4.  
%                0 or no 3rd argument = false (returns nx3 with XYZ data only)
%                1 = true (returns nx4 with XYZ and residuals) 
%   index1     = start frame index, all frames if not used as an argument
%   index2     = end frame index, all frames if not used as an argument
%   OUTPUTS:
%   XYZPOS     = nx3/4 matrix with n frames and X, Y, Z, and/or residual as columns
   
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
  
    index1 = itf.GetVideoFrame(0); % frame start
    index2 = itf.GetVideoFrame(1); % frame end

elseif nargin == 0,
    disp('Error: wrong number of inputs.');
    help getanalogchannel;
    return
end
 nIndex = itf.GetParameterIndex('POINT', 'LABELS');
 if nIndex==-1,
     XYZPOS=zeros(index2,8,3);
     MARKERS={'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_left','head_front','head_right'};
     return
 end
 nItems = itf.GetParameterLength(nIndex);

    
for i = 1 : nItems,
    target_name = itf.GetParameterValue(nIndex, i-1);
    marker_name = target_name(min(findstr(target_name, ':'))+1:end);
    
   %fprintf('%s\n',marker_name);
  
    MARKERS{i} = marker_name;
    XYZPOS(:,i,1) = cell2mat(itf.GetPointDataEx(i-1,0,index1,index2,'1'));
    XYZPOS(:,i,2) = cell2mat(itf.GetPointDataEx(i-1,1,index1,index2,'1'));
    XYZPOS(:,i,3) = cell2mat(itf.GetPointDataEx(i-1,2,index1,index2,'1'));
   
end


function createc3d(itf,name,vfrt,nfram,nmkr,avr,achn,ftype,dtype,pscal)
% CREATEC3D - creates C3D file with the input parameters.
% 
% USAGE: createc3d(itf,name,vfrt*,nfram*,nmkr*,avr*,achn*,ftype*,dtype*,pscal*)
%        * = not a necessary input
% INPUTS:
%   itf      = variable name used for COM object
%   name     = name given to C3D file
%   vfrt     = video frame rate, must be > 0, default 120
%   nfram    = number of frames, default 500
%   nmkr     = number of markers, not < 0, default 20
%   avr      = analog to video ratio, default 13
%   achn     = number of analog channels, not < 0, default 28 
%   ftype    = Intel(1), DEC(2), SGI(3), default 1
%   dtype    = Integer(1), Floating Point(2), default 2
%   pscal    = point data scaling factor, must be > 0, default 0.1 
%
%   Note: Only requires 2 inputs ('itf' and 'name') if user accepts
%   defaults for other parameters.  createc3d will handle 'name' with or 
%   without the '.c3d' extension.  'name' can be entered with
%   directory path (to create in a directory different from your 
%   current Matlab path).
    
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


if nargin==9,pscal=0.1;
elseif nargin==8,pscal=0.1;dtype=2;
elseif nargin==7,pscal=0.1;dtype=2;ftype=1;
elseif nargin==6,pscal=0.1;dtype=2;ftype=1;achn=28;
elseif nargin==5;pscal=0.1;dtype=2;ftype=1;achn=28;avr=13;
elseif nargin==4;pscal=0.1;dtype=2;ftype=1;achn=28;avr=13;nmkr=20;
elseif nargin==3;pscal=0.1;dtype=2;ftype=1;achn=28;avr=13;nmkr=20;nfram=500;
elseif nargin==2;pscal=0.1;dtype=2;ftype=1;achn=28;avr=13;nmkr=20;nfram=500;vfrt=120;
end
 dotind = findstr(name,'.');
if isempty(dotind) == 1, name = [name '.c3d']; end
pRet = itf.NewFile(name,ftype,dtype,achn,avr,nmkr,vfrt,pscal,nfram);
if pRet == 1, disp([name ' has been created.']); end
%--------------------------------------------------------------------------
    
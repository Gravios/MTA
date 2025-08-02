

fid = fopen('/storage/gravio/data/raw/vicon/Vicon_Down/Top Level 2/jg05/jg05-20120317/cof/Trial001_0001.c3d','r');
fid = fopen('/storage/gravio/data/raw/motive/ER11-20150520-cof3d.c3d','r')

fseek(fid,0,-1)
c3d.hexP2ParmSection = fread(fid,2,'uint8');
c3d.numOfStoredTrajs = fread(fid,1,'uint16');
c3d.numOfAnalogMperframe = fread(fid,1,'uint16');
c3d.numOfFirstFrame = fread(fid,1,'uint16');
c3d.numOfLastFrame = fread(fid,1,'uint16');
c3d.MaxInterpGap = fread(fid,1,'uint16');
c3d.scaleFactor = fread(fid,1,'float32');
c3d.dataStart = fread(fid,1,'uint16');
c3d.numOfAnalogSamples = fread(fid,1,'int16');
c3d.sampleRate = fread(fid,1,'float32');


fseek(fid,512,-1)
c3d.pramRes = fread(fid,2,'uint8');
c3d.numOfParmBlocks = fread(fid,1,'uint8');
c3d.procType = fread(fid,1,'uint8');

fseek(fid,516,-1)
fseek(fid,512,-1)
h = '';
for i = 1:512,
    %fread(fid,1,'uint8')
    h(i) = fread(fid,1,'*char');
end

% Parameter Blocks are 512 bytes

fseek(fid,1024,-1)
parm.b1 = fread(fid,1,'uint8');
parm.b2 = fread(fid,1,'uint8');
parm.b3 = fread(fid,1,'uint8');
parm.b4 = fread(fid,1,'uint8');

parm.b1 = fread(fid,1,'*char');
parm.b2 = fread(fid,1,'*char');
parm.b3 = fread(fid,1,'*char');
parm.b4 = fread(fid,1,'*char');


h1 = '';
for i = 1:512,
    %fread(fid,1,'uint8')
    h1(i) = fread(fid,1,'*char');
end


fclose(fid);



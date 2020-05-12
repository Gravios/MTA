

map = 

  struct with fields:

        chanMap: [16x1 double]
        xcoords: [16x1 double]
        ycoords: [16x1 double]
    chanMap0ind: [16x1 double]
      connected: [16x1 logical]
        kcoords: [16x1 double]
           name: 'Linear 16x1'


% MAP setup
clear map           

% CHANMAP 
map.chanMap = [1:96];

% XCOORDS 
buz64x = reshape(repmat(0:200:200*7,[8,1]),[],1);
lin32x = [buz64x(end)+500].*ones([32,1]);
map.xcoords = [buz64x;lin32x];

% YCOORDS 
buz64y = repmat(0:20:(20*7),[1,8])';
lin32y = [0:50:(50*31)]';
map.ycoords = [buz64y;lin32y];

           
% CHANMAP0IND 
map.chanMap0ind = map.chanMap-1;

% CONNECTED 
map.connected  = ones([96,1]);

buz64k = reshape(repmat(1:8,[8,1]),[],1);
lin32k = ones([32,1]).*(max(buz64k)+1);
map.kcoords = [buz64k;lin32k];
map.name = 'buz64lin32';

save('/storage/share/matlab/Kilosort2/configFiles/buz64lin32_kilosortChanMap.mat','-struct','map');





ds = load('/storage/gravio/code/Kilosort2/configFiles/neuropixPhase3B1_kilosortChanMap.mat');


% MAP setup
clear map           

% CHANMAP 
map.chanMap = [1:64];

% XCOORDS 
buz64x = reshape(repmat(0:500:500*7,[8,1]),[],1);
map.xcoords = [buz64x];

% YCOORDS 
buz64y = repmat(0:20:(20*7),[1,8])';
map.ycoords = [buz64y];

           
% CHANMAP0IND 
map.chanMap0ind = map.chanMap-1;

% CONNECTED 
map.connected  = true([numel(map.chanMap),1]);

buz64k = reshape(repmat(0:7,[8,1]),[],1);
map.kcoords = [buz64k];
map.name = 'buz64';
map.fs = 32552;

save('/storage/gravio/code/Kilosort2/configFiles/buz64_kilosortChanMap.mat','-struct','map');

ds = load('/storage/share/matlab/Kilosort2/configFiles/buz64lin32_kilosortChanMap.mat');


%% testing 8 chan -----------------------------------------------------------------------------------------------------
% MAP setup
clear map           

% CHANMAP 
map.chanMap = [57:64];

% XCOORDS 
map.xcoords = [repmat(0,[8,1])];

% YCOORDS 
map.ycoords = [0:20:(20*7)]';

% CHANMAP0IND 
map.chanMap0ind = map.chanMap-1;

% CONNECTED 
map.connected  = true([numel(map.chanMap),1]);

map.kcoords = zeros([8,1]);
map.name = 'buz8';
map.fs = 32552;

save('/data/gravio/code/Kilosort2/configFiles/buz64_kilosortChanMap.mat','-struct','map');

ds = load('/storage/gravio/code/Kilosort2/configFiles/buz64lin32_kilosortChanMap.mat');



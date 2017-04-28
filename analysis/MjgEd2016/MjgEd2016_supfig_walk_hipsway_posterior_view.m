
%% supfig Lateral sway of the hips during walk
%videofile = ['/storage/share/exchange/gravio/optitrack/Session\ 2017-01-12/Take\ 2017-01-12\ 02.52.14\ PM-Camera\ 13029.avi'];
%videofile = '/storage/gravio/data/project/general/JM11-20170112/Trial001_C13029.avi';
videofile_side = '/storage/gravio/data/project/general/JM11-20170112/Trial001_C12881.avi';
%vdr = VideoReader(videofile);
vdr_side = VideoReader(videofile_side);

% Plot Timepoints
hfig = figure(20170131);clf
set(hfig,'Units','centimeters',...
         'PaperPositionMode', 'auto',...
         'Position',[0,0,21,30])
index = 74400;


vdw = VideoWriter(['/storage/gravio/ownCloud/videos/'...'
                   'rat_walk_hip_sway.avi'],...
                  'Uncompressed AVI');
vdw.FrameRate = 180;
vdw = VideoWriter(['/storage/gravio/ownCloud/videos/'...'
                   'rat_walk_hip_sway_slow_3x.avi'],...
                  'Uncompressed AVI');
vdw.FrameRate = 60;

vdw = VideoWriter(['/storage/gravio/ownCloud/videos/'...'
                   'rat_walk_hip_sway_slow_6x.avi'],...
                  'Uncompressed AVI');
vdw.FrameRate = 30;

vdw.open;

% $$$ timePoints = [0,15,45,75,85];
%imageCloseUpBox = [195,600;170,340]
for i = 1:460;
    vdr_side.CurrentTime = (index+i)/vdr_side.FrameRate;
    currentFrame = rgb2gray(readFrame(vdr_side));
    writeVideo(vdw,currentFrame);
% $$$     
% $$$     hax = axes('Units','centimeters',...
% $$$                'Position',[2+3*(i-1)+i/5,20,3,1.515]);
% $$$     imagesc(rot90(vs(imageCloseUpBox(2,1):imageCloseUpBox(2,2),...
% $$$                      imageCloseUpBox(1,1):imageCloseUpBox(1,2),1))');
% $$$     hax.XTickLabel = {};
% $$$     hax.YTickLabel = {};
end

vdw.close;








vdw = VideoWriter(['/storage/gravio/ownCloud/videos/'...'
                   'rat_walk_hip_sway_small.avi'],...
                  'Uncompressed AVI');
vdw.FrameRate = 180;

vdw = VideoWriter(['/storage/gravio/ownCloud/videos/'...'
                   'rat_walk_hip_sway_slow_small_3x.avi'],...
                  'Uncompressed AVI');
vdw.FrameRate = 60;

vdw = VideoWriter(['/storage/gravio/ownCloud/videos/'...'
                   'rat_walk_hip_sway_slow_small_6x.avi'],...
                  'Uncompressed AVI');
vdw.FrameRate = 30;

vdw.open;

% $$$ timePoints = [0,15,45,75,85];
%imageCloseUpBox = [195,600;170,340]
for i = 1:460;
    vdr_side.CurrentTime = (index+i)/vdr_side.FrameRate;
    currentFrame = rgb2gray(readFrame(vdr_side));
    writeVideo(vdw,currentFrame(150:350,1:200));
% $$$     
% $$$     hax = axes('Units','centimeters',...
% $$$                'Position',[2+3*(i-1)+i/5,20,3,1.515]);
% $$$     imagesc(rot90(vs(imageCloseUpBox(2,1):imageCloseUpBox(2,2),...
% $$$                      imageCloseUpBox(1,1):imageCloseUpBox(1,2),1))');
% $$$     hax.XTickLabel = {};
% $$$     hax.YTickLabel = {};
end

vdw.close;
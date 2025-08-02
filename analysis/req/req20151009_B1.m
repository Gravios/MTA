function req20151009_B1(Trial,varargin)
[hfig,state,stcColor] = DefaultArgs(varargin,{gcf,'rear','r'});

% FXYZ Filtered Marker Positions {Low Pass 2.5 Hz}
fxyz = Trial.load('xyz');
fxyz.filter('ButFilter',3,2.5,'low');

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);

%% Inter animal JPDF labeled contours

fet = fxyz.copy;
fet.data = [fang(:,1,4,2),[log10(abs(diff(fang(:,3,4,2))));0]];
edgs    = {linspace(0,1.4,75)};
edgs(2) = {linspace(-6,-1,75)};
edc = edgs;
[edc{:}] = get_histBinCenters(edc);
[X,Y] = meshgrid(edc{:});

figure(hfig);

% Feature and states of jg05-20120317 
if strcmp(Trial.name,'jg05-20120317'),
    % JPDF - Head/Body speed
    %ind = Trial.stc{'a'};
    ind = nniz(fet);
    b = fet(ind,:);
    hist2(b,edgs{1},edgs{2});
    xlabel('log10 body pitch (radians)');
    ylabel('log10(abs(d(BMBU_p_i_t_c_h)/dt)) log10(rad/sec)');
    title({'JPDF of log10(abs(d(BMBU_p_i_t_c_h)/dt)) VS BMBU_p_i_t_c_h',...
           [Trial.filebase ': overlayed with jg05-20120317 labeled states']});
    stc = Trial.load('stc','hand_labeled_rev2');
else
    stc = Trial.load('stc','LGR-hand_labeled_rev2-wrnpms');
end

hold on,
b = fet(stc{state},:);
o = hist2(b,edgs{1},edgs{2});
F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
o = conv2(o,F,'same');
contour(X,Y,o',[10,10],'linewidth',1.5,'Color',stcColor)




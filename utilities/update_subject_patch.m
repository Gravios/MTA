function subject = update_subject_patch(subject,appendage,headBodyAngIndex,overlay,hbaBinEdg,hbaBinCtr)
% SET rotation and arc boundaries
    if ~isempty(headBodyAngIndex);
        rotHba = -hbaBinCtr(headBodyAngIndex);
        rotBnd = -linspace([hbaBinEdg(headBodyAngIndex:headBodyAngIndex+1),31])';
    else,
        rotHba = [];
        rotBnd = [];
    end

% DESIGNATE midline
    subject.(appendage).midline.vert = ...
        {[0,subject.(appendage).patch.vert{1}(end-floor(numel(subject.(appendage).patch.vert{1})/2))],...
         [0,subject.(appendage).patch.vert{2}(end-floor(numel(subject.(appendage).patch.vert{2})/2))]};
    subject.(appendage).midline.color = 'r';
    
% HEAD - hba range overlay
    if ~isempty(rotBnd) && overlay, 
        verts = [subject.(appendage).midline.vert{1}(2);    ...
                 subject.(appendage).midline.vert{2}(2)]    ...
                .* 1.2; % scaled
        subject.(appendage).overlay.vert{1} = ...
            [0,multiprod([cos(rotBnd),-sin(rotBnd)], verts, 2, 1)'];
        subject.(appendage).overlay.vert{2} = ...
            [0,multiprod([sin(rotBnd), cos(rotBnd)], verts, [1,2], 1)'];
    end
    
% COMPUTE rotated patch verticies for appendage
    if ~isempty(headBodyAngIndex)
        verts = cat(1,subject.(appendage).patch.vert{:});
        subject.(appendage).patch.vert{1} = ...
            sq(multiprod([cos(rotHba),-sin(rotHba)], verts, 2, 1));
        subject.(appendage).patch.vert{2} = ...
            sq(multiprod([sin(rotHba), cos(rotHba)], verts, 2, 1));
    end
    
% UPDATE midline
    subject.(appendage).midline.vert = ...
        {[0,subject.(appendage).patch.vert{1}(end-floor(numel(subject.(appendage).patch.vert{1})/2))],...
         [0,subject.(appendage).patch.vert{2}(end-floor(numel(subject.(appendage).patch.vert{2})/2))]};
    subject.(appendage).midline.color = 'r';

end



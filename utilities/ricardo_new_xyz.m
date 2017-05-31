function ricardo_new_xyz(Session)

xyz = Session.load('xyz');

xyz.addMarker('head_front',...     Name
    [.7,0,.7],...  Color
    {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
     {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
     {'head_right','hcom',[0,0,255]},... new marker to skeleton
     {'head_FL','hcom',[0,0,255]},...
     {'head_FR','hcom',[0,0,255]}},... 
        mean(xyz(:,{'head_FL','head_FR'},:),2));

xyz.save
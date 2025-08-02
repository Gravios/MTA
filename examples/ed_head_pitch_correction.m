
RefTrial = MTATrial.validate('jg05-20120317.cof.all');


%T = {'Ed03-20140624.cof.all',...}

Trial = MTATrial.validate('Ed03-20140624.cof.all');
fet = fet_head_pitch(Trial);
fet.map_to_reference_session(Trial,RefTrial);







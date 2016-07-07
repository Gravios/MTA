% Need to correct marker shift error in sessions {'er01-20110719','Ed01-20140707'}
% marker shift error occurs when a marker must be reattached mid-session
% 
% MOD : MTA/@MTADfet/map_to_reference_session.m

Trial = MTATrial('Ed01-20140707');
Trial.load('stc','hand_labeled_rev2_Ed');
xyz = Trial.load('xyz');
xyz = Trial.load('xyz','seh');
ang = create(MTADang,Trial,xyz);

%% Checking ANG
eds = log10(linspace(40,80,200));
ind = Trial.stc{'w'};
figure,bar(eds,histc(log10(ang(ind,1,2,3)),eds),'histc');
% Okay

eds = log10(linspace(80,130,200));
ind = Trial.stc{'w'};
figure,bar(eds,histc(log10(ang(ind,1,3,3)),eds),'histc');
% Okay

eds = log10(linspace(40,80,200));
ind = Trial.stc{'w'};
figure,bar(eds,histc(log10(ang(ind,2,3,3)),eds),'histc');
% Okay

eds = log10(linspace(40,80,200));
ind = Trial.stc{'w'};
figure,bar(eds,histc(log10(ang(ind,3,4,3)),eds),'histc');
% Okay


eds = linspace(-pi/2,pi/2,200);
ind = Trial.stc{'w'};
figure,bar(eds,histc(ang(ind,1,2,2),eds),'histc');
% Okay
eds = linspace(-pi/2,pi/2,200);
ind = Trial.stc{'w'};
figure,bar(eds,histc(ang(ind,1,3,2),eds),'histc');
% Okay

eds = linspace(-pi/2,pi/2,200);
ind = Trial.stc{'w'};
figure,bar(eds,histc(ang(ind,2,3,2),eds),'histc');
% Okay

eds = linspace(-pi/2,pi/2,200);
ind = Trial.stc{'w'};
figure,bar(eds,histc(ang(ind,3,4,2),eds),'histc');
% Okay




%% Checking XYZ

ind = Trial.stc{'w'};

eds = log10(linspace(10,100,200));
figure,bar(eds,histc(log10(xyz(ind,1,3)),eds),'histc');
% Okay

eds = log10(linspace(50,150,200));
figure,bar(eds,histc(log10(xyz(ind,2,3)),eds),'histc');
% Okay

eds = log10(linspace(50,150,200));
figure,bar(eds,histc(log10(xyz(ind,4,3)),eds),'histc');
% Okay

eds = log10(linspace(20,150,200));
figure,bar(eds,histc(log10(xyz(ind,5,3)),eds),'histc');
% Okay



function jpdf_bv_hz(Sessions)

if isempty(Sessions)

elseif ischar(Sessions),


Sessions = {{'er01-20110719',     'all', 'bach'},... CA3
               {'er01-20110721',     'all', 'bach'},... CA3
               {'er06-20130612',     'all',  'cin'},... CA1
               {'er06-20130613', 'all-cof',  'cin'},... CA1
               {'er06-20130614', 'all-cof',  'cin'},... CA1
               {'jg04-20120129',     'all', 'bach'},... CA1
               {'jg04-20120130',     'all', 'bach'},... CA1
               {'jg05-20120309',     'all', 'bach'},... CA1
               {'jg05-20120310',     'all', 'bach'},... CA1
               {'jg05-20120311',     'all', 'bach'},... CA1
               {'jg05-20120315',     'all', 'bach'},... CA1
               {'jg05-20120317',     'all', 'bach'},... CA2???
               {'jg05-20120324',     'all', 'bach'}}; % CA3

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.5,xyz.sampleRate));
bv = MTADxyz('data',log10([0;Trial.vel('spine_lower',[1,2])]),'sampleRate',xyz.sampleRate);
hz = MTADxyz('data',log10(xyz(:,7,3)),'sampleRate',xyz.sampleRate);

plot_distrbs(Sessions,
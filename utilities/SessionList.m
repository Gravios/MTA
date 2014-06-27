function Sessions = SessionList(varargin)


path = load('MTAPaths.mat');

switch varargin{1}
  case 'all'
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

end
path.data
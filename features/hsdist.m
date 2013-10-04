function feature=hsdist(Session,method,varargin)
Session = Session.load_ang;
feature = Session.ang(:,4,5,3);
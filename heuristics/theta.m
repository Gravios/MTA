function [data,sampleRate,label,key] = theta(Session,varargin)
[mode,label,key] = DefaultArgs(varargin,{'sts2epoch','theta','t'});


switch mode
% $$$   case 'calculate'
% $$$         s = load([Trial.name '.eegseg.mat']);
% $$$         [Periods,frin,frout,ch2use] = DefaultArgs(varargin,{ [s.t(1) s.t(end)], [5 11], [1 5; 11 15],1});
% $$$ 
% $$$         %find freq. indexes within and outsize the band of interest
% $$$         thfin = logical(WithinRanges(s.f,frin));
% $$$         thfout = logical(WithinRanges(s.f,frout));
% $$$ 
% $$$         timein = logical(WithinRanges(s.t,Periods));
% $$$         timeout = ~timein; %setdiff([1:length(s.t)],timein);
% $$$         thratio = log(mean(sq(s.y(timein,thfin,ch2use)),2))-log(mean(sq(s.y(timein,thfout,ch2use)),2));
% $$$         nStates =2;
% $$$ 
% $$$         % fit gaussian mixture and to HMM - experimental version .. uses only
% $$$         % thetaratio 
% $$$         [theState thhmm thdec] = gausshmm(thratio,nStates,1,0);
% $$$ 
% $$$         %find which state is theta  and make step vector from 0 to 1 (outside to
% $$$         %inside)
% $$$         for ii=1:nStates 
% $$$             thratio_st(ii) = mean(thratio(find(theState==ii)));
% $$$         end
% $$$         [dummy TheInd] = max(thratio_st);
% $$$         InTh = zeros(length(s.t),1);
% $$$         InTh(timein) = (theState==TheInd);
% $$$         InTh(timeout) = 0; % that puts all times outside the periods as not theta
% $$$ 
% $$$         %detect beg and end of the theta runs
% $$$         thebeg= SchmittTrigger(InTh,0.9, 0.9);
% $$$         theend= SchmittTrigger(-InTh,-0.9, -0.9);
% $$$ 
% $$$         if thebeg(1)>theend(1)
% $$$             theend =theend(2:end);
% $$$         end
% $$$ 
% $$$         if thebeg(end)>theend(end)
% $$$             thebeg =thebeg(1:end-1);
% $$$         end
% $$$         %theend = theend-1; %shif tend back as the it is detected on donw/up
% $$$ 
% $$$         thrun =round([s.t(thebeg) s.t(theend)]*Trial.lfp.sampleRate);
% $$$ 
% $$$         if nargout<1
% $$$             figure
% $$$             myf= [min([frin(:);frout(:)]) max([frin(:);frout(:)])];
% $$$             myfi = find(f>myf(1) & f<myf(2));
% $$$             imagesc(t,f,log(sq(y(:,myfi,1)))');
% $$$             axis xy
% $$$             ca = caxis; caxis([ca(1) 0.8*ca(2)]);
% $$$             hold on
% $$$             h1 = Lines(thrun(:,1)/1250,[],'b');
% $$$             h2 = Lines(thrun(:,2)/1250,[],'r');
% $$$             h3 = Lines(Periods(:),[],'k');
% $$$             set(h1,'LineWidth',2);
% $$$             set(h2,'LineWidth',2);
% $$$             set(h3,'LineWidth',2);
% $$$             TimeBrowse(100,100);
% $$$             %   msave([FileBase '.thrun'],'thrun');
% $$$         end
% $$$ 
% $$$         % nSeg = size(thebeg,1);
% $$$         % for ii=1:nSeg
% $$$         %     thsegm(ii) = mean(thratio([thebeg(ii):theend(ii)]));
% $$$         %     thsegv(ii) = std(thratio([thebeg(ii):theend(ii)]));
% $$$         % end
% $$$         %     
% $$$         % pow = log(sq(y(:,myfi,1)));
% $$$         % mpow = mean(pow(find(InTh),:));
% $$$         % for ii=1:nSeg
% $$$         %     md = mahal(mpow,pow([thebeg(ii):theend(ii)],:));
% $$$         %     
% $$$         %     
% $$$         %     
% $$$         return


  case 'sts2double'
        sampleRate = Session.lfp.sampleRate;
        
        data = load(fullfile(Session.spath,[Session.name '.sts.theta']));
        
  
  case 'sts2epoch'
      data = load(fullfile(Session.spath,[Session.name '.sts.theta']));
      sync = Session.lfp.sync.copy;
        lsync = sync.sync.copy;
        lsync.resample(Session.lfp.sampleRate);
        data = IntersectRanges(lsync.data,data)-lsync.data(1)+1;
        data = MTADepoch(Session.spath,...
                        Session.filebase,...
                        data,...
                        Session.lfp.sampleRate,...
                        sync,...
                        sync(1),...
                        label,key);
end

end
    
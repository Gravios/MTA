function [rhm,varargout] = fet_rhm_exp(Trial,varargin)
% [rhm,fs,ts] = fet_rhm(Trial,varargin)
% [sampleRate,mode,windowSize] = DefaultArgs(varargin,{Trial.xyz.sampleRate,'spectral',1});
% Need to update the spectral window size to adapt to the xyz.sampleRate
%
% 
%

parspec = empty_spec;
xyz = Trial.load('xyz');
%xyz = Trial.load('xyz','seh');
%xyz.data(isnan(xyz.data(:)))=0;
varargout = cell([1,nargout-1]);

[sampleRate,mode,wsig,defspec,overwrite,newSR,type] = DefaultArgs(varargin,{'xyz','mta',1,def_spec_parm(xyz),false,[],'mta'});

fs = []; ts = [];

% create a ridgid body model
rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
% find the center of mass of the model
hcom = xyz.com(rb);
% add coordinates of the model's center of mass to the xyz object
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,3,[1]./(xyz.sampleRate/2),'low'));
xyz.addMarker('hcom',[128,255,128],{{'head_back','head_front',[0,0,1]}},hcom);

offset = 450;

nz = -cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_left',:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
nm = nz.*offset+hcom;
xyz.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

ny = cross(xyz(:,'htx',:)-hcom,xyz(:,'head_back',:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny.*offset+hcom;
xyz.addMarker('hrx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

nx = cross(xyz(:,'hrx',:)-hcom,xyz(:,'htx',:)-hcom);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx.*offset+hcom;    
xyz.addMarker('hbx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);


nx = cross(xyz(:,'hrx',:)-hcom,xyz(:,'htx',:)-hcom);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx.*offset*sqrt(2)+hcom;    
xyz.addMarker('hbxe',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);
xyz.addMarker('hbte',[128,255,128],{{'head_back','head_front',[0,0,1]}},...
                  genRotatedMarker(xyz,'hbxe',45,{'hbx','htx'}));
xyz.addMarker('fhbte',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);
xyz.data(nniz(xyz(:,'hbte',:)),end,:) = ButFilter(xyz(nniz(xyz(:,'hbte',:)),'hbte',:),...
                                                 3,[1]./(xyz.sampleRate/2),'low');


rb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle'});
bcom = xyz.com(rb);
xyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(bcom,4,[1]./(xyz.sampleRate/2),'low'));
xyz.addMarker('bcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},bcom);
xyz.addMarker('fsl',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_lower',:),4,[1.5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('fsm',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_middle',:),4,[1.5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('fsu',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_upper',:),4,[1.5]./(xyz.sampleRate/2),'low'));

rb = xyz.model.rb({'pelvis_root','spine_middle','spine_upper'});
bcom = xyz.com(rb);
xyz.addMarker('fbucom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(bcom,4,[1]./(xyz.sampleRate/2),'low'));
xyz.addMarker('bucom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},bcom);
xyz.addMarker('fsl',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_lower',:),4,[1.5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('fsm',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_middle',:),4,[1.5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('fsu',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(xyz(:,'spine_upper',:),4,[1.5]./(xyz.sampleRate/2),'low'));


% if xyz sampling rat e is greater than 120 Hz then resample it to 120 Hz
if ischar(sampleRate), sampleRate = Trial.(sampleRate).sampleRate;end
if isa(sampleRate,'MTAData'),
    xyz.resample(sampleRate);    
elseif sampleRate > 120, 
    xyz.resample(120); 
end
xyz.filter('ButFilter',3,30,'low');


ang = create(MTADang,Trial,xyz);
ang.data(~nniz(xyz(:,1,1)),:,:,:) = 0;
hang = Trial.transform_origin;
hang.roll(~nniz(xyz)) = 0;
fet = xyz.copy;
%bang = ButFilter(diff(hang.roll),3,[2,30]./(xyz.sampleRate/2),'bandpass');
bang = ang(:,'hbx','fhcom',3);
%bang = ButFilter(ang(:,'fsm','spine_upper',3),4,[0.8,20]./(xyz.sampleRate/2),'bandpass');
%bang = ButFilter(ang(:,'fbcom','spine_upper',3)-ang(:,'spine_upper','head_back',3),4,[0.5,20]./(xyz.sampleRate/2),'bandpass');
%bang = ButFilter(ang(:,'fbcom','spine_upper',3),4,[0.5,20]./(xyz.sampleRate/2),'bandpass');
%bang = ButFilter(ang(:,'fbcom','spine_upper',3),4,[0.5,20]./(xyz.sampleRate/2),'bandpass');
%bang = ButFilter(ang(:,'fbcom','spine_upper',3),4,[0.8]./(xyz.sampleRate/2),'high');
%bang = ButFilter(ang(:,'bucom','fsl',3),4,[0.8]./(xyz.sampleRate/2),'high');
%bang = ang(:,'hbx','fhbte',3);
%bang = ang(:,'hbx','hcom',2);
%bang = ang(:,'hbt','fhcom',3);

%fet.data = bang;
%fet.data = [0;diff(bang)];
fet.data = [0;ButFilter(diff(bang),3,[.1,20]/(ang.sampleRate/2),'bandpass')];

switch mode
  case 'mta'
    rhm = MTADfet.encapsulate(Trial,...
                              fet.data,...
                              fet.sampleRate,...
                              'rhythmic head motion feature',...
                              'rhm',...
                              'r');
  case 'raw'
    rhm = fet.data;
  otherwise
    dsf = fieldnames(defspec);
    for i = 1:length(dsf),parspec.(dsf{i}) = defspec.(dsf{i});end
    data = zeros(fet.size);
    
    if wsig,
        try,load(fullfile(Trial.path.arm,[mfilename,'.arm.mat']));end

        if exist('ARmodel','var')||overwrite,
            data(nniz(fet.data),:) = WhitenSignal(fet.data(nniz(fet.data),:),...
                                                      [],...
                                                      true,...
                                                      ARmodel);
        else
            [data(nniz(fet.data),:),ARmodel] = WhitenSignal(fet.data(nniz(fet.data),:),...
                                                                [],...
                                                                true);
            save(fullfile(Trial.path.arm,[mfilename,'.arm.mat']),'ARmodel');
        end
    else
        data(nniz(fet.data),:) = fet.data(nniz(fet.data),:);
    end

     svout = cell([1,nargout-3]);
    [ys,fs,ts,svout{:}] = spec(str2func(mode),data,parspec);


    % Modify time stamps and spec; add padding (0's)
    ts = ts+(parspec.WinLength/2)/fet.sampleRate;
    ssr = 1/diff(ts(1:2));
    pad = round([ts(1),mod(fet.size(1)-round(parspec.WinLength/2),parspec.WinLength)/fet.sampleRate].*ssr)-[1,0];
    szy = size(ys);
    rhm = MTADlfp('data',cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)])),'sampleRate',ssr);

    ts = cat(2,([1:pad(1)]./ssr),ts',([pad(1)+size(ts,1)]+[1:pad(2)])./ssr)';

    if numel(svout)>0,
        for i = 1:numel(svout),
        svout{i} = cat(1,zeros([pad(1),size(svout{i},2),size(svout{i},3),size(svout{i},4)]),...
                         svout{i},...
                         zeros([pad(2),size(svout{i},2),size(svout{i},3),size(svout{i},4)]));
        end
    end

    if ~isempty(newSR),
        rhm.resample(newSR);
        temp_ts = Trial.xyz.copy;
        temp_ts.data = ts;
        temp_ts.resample(newSR);
        ts = temp_ts.data;
    end

    if strcmp(type,'raw'),
        rhm = rhm.data;
    end
    
    
    tvout = cat(2,fs,ts,svout);
    varargout = tvout(1:length(varargout));

end







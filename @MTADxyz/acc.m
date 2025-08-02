function a = acc(Data,varargin)
%a = acc(Data,varargin)
%calculate the acceleration of marker(s)
%[marker,dim] = DefaultArgs(varargin,{[1:Data.model.N],[1:size(Data.xyz,3)]});
[markers,dims,padded] = DefaultArgs(varargin,{1:Data.model.N, 1:Data.size(3),1});
vel = Data.vel(markers,dims);
a = diff(vel.data,1);
if padded==1
    a = MTADxyz('data',cat(1,a,a(end,:)),'sampleRate',Data.sampleRate);
end
a.model = vel.model.copy;
end

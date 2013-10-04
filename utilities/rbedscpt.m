
sl

hrb = s.xyz(:,s.Model.gmi('head_right'),:)-s.xyz(:,s.Model.gmi('head_back'),:);
hlb = s.xyz(:,s.Model.gmi('head_left'),:)-s.xyz(:,s.Model.gmi('head_back'),:);
crlb = sq(cross(hlb,hrb));
crlbt = reshape(repmat(sqrt(sum(crlb.^2,2)),3,1),[],3);
trans_xyzb = crlb./crlbt;

hrf = s.xyz(:,s.Model.gmi('head_right'),:)-s.xyz(:,s.Model.gmi('head_front'),:);
hlf = s.xyz(:,s.Model.gmi('head_left'),:)-s.xyz(:,s.Model.gmi('head_front'),:);
crlf = sq(cross(hlf,hrf));
crlft = reshape(repmat(sqrt(sum(crlf.^2,2)),3,1),[],3);
trans_xyzf = crlf./crlft;

hcom = s.com(s.Model.rb({'head_back','head_left','head_front','head_right'}));


rbinteg = sqrt(sum((trans_xyzb+trans_xyzf).^2,2));
plot(rbinteg)



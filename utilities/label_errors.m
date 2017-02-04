function [Errors2]=label_errors(Trial,period)

% Note that period is a value in seconds to add extra samples before and
% after the errors


%% Head errors
    if (nargin<2 || isempty(period)), period = 2; end

    errors= FindErrorPeriods(Trial);


    if isempty(Trial.xyz.data)
        Trial.load('xyz')
    end
    
    BinaryErrors = WithinRanges(1:length(Trial.xyz.data), errors);
    gauss = normpdf(linspace(-1,1,round(period * Trial.xyz.sampleRate)),0,1);

    ErrorsCorrected= (conv(double(BinaryErrors),gauss,'same')>.5);
    Errors1 = ErrorsCorrected;
    BinaryErrors=[];

    %% Body errors

    BodyRef = zeros(length(Trial.xyz.data),3,3);

    % Extracting vectors of the body 
    for ii=2:4; % working just with body markers including back of the head
        RefVector= bsxfun (@minus,Trial.xyz.data(:,ii,:),Trial.xyz.data(:,ii-1,:));
        RefUnit = bsxfun(@rdivide,RefVector,sqrt(sum(RefVector.^2,3))); % unit vector normalization
        distance=[];
        distance=(sqrt(nansum(RefVector(:,1,:).^2,3)));
        distance(distance==0)=nan; % note that we introduce nan from here to be consistent later
        distanceNorm= (distance-nanmean(distance))/nanvar(distance); 
        distanceNorm= (distance-nanmean(distance))/nanvar(distance); 

        [Az,Pit,xx]=cart2sph(RefUnit(:,1,1),RefUnit(:,1,2),RefUnit(:,1,3));

        BodyRef(:,:,ii-1)=[Az, Pit ,distanceNorm];
    end
    RefVector=[];
    RefUnit=[];
    xx=[];

    % calculating angular and distances differences for sections of the body
    AngDiffBody=zeros(length(Trial.xyz.data),3,2);
    for ii =2:3;

        AngDiffBody(:,1,ii-1)= circ_dist(BodyRef(:,1,ii-1),BodyRef(:,1,ii));
        AngDiffBody(:,2,ii-1)= circ_dist(BodyRef(:,2,ii-1),BodyRef(:,2,ii));
        AngDiffBody(:,3,ii-1)= nanvar([BodyRef(:,3,ii-1),BodyRef(:,3,ii)],0,2);% we are using just the variance between two body segments distances
    end

    %% Distance variance Normalized
    DistVar=zeros(length(Trial.xyz.data),1);

    DistVar =tiedrank(sum(AngDiffBody(:,3,:),3));

    % The normalization 
    DistVar = (DistVar/max(DistVar));

    %% Calculating the resultant length as a messure of variance of the angles
    % Azimuth
    X = sin(sq(AngDiffBody(:,1,:)));
    Y = cos(sq(AngDiffBody(:,1,:)));

    Xm = nanmean(X,2);
    Ym = nanmean(Y,2);

    AzimuthRL=zeros(length(Trial.xyz.data),1);
    AzimuthRL= sqrt(sum([Xm.^2 Ym.^2],2));
    AzimuthRL = AzimuthRL-DistVar;
    % Note that the resultant length value is weighted by the variance of the
    % distances (the general feature of body segments is that they show high angular clustering, that's why the minus sign)
    % mAng = atan2(Ym,Xm); % mean angle (not used at the moment)

    % Pitch
    Xp = sin(sq(AngDiffBody(:,2,:)));
    Yp = cos(sq(AngDiffBody(:,2,:)));

    Xmp = nanmean(Xp,2);
    Ymp = nanmean(Yp,2);

    PitchRL=zeros(length(Trial.xyz.data),1);
    PitchRL= sqrt(nansum([Xmp.^2 Ymp.^2],2));
    PitchRL = PitchRL-DistVar; % weighted R 
                               % mAngp = atan2(Ymp,Xmp);

    % To avoid detection of small changes
    medAzimuthRL=medfilt1(AzimuthRL,round(2*Trial.xyz.sampleRate));
    medPitchRL=medfilt1(PitchRL,round(2*Trial.xyz.sampleRate));

    % Adding the detection for Azimuth and Pitch (thresholds acoording with the distibutions observed in all Animals) 
    BodyErrors = sum([medAzimuthRL<-.15 medPitchRL<-.15],2)>0;

    gauss = normpdf(linspace(-1,1,round(period* Trial.xyz.sampleRate)),0,1);
    gauss=gauss./sum(gauss);


    BodyErrors = conv(double(BodyErrors),gauss,'same')>.5;
    
    Errors2 = ThreshCross((double(Errors1)+double(BodyErrors))>0,0.5,1);

    Trial.stc.states(Trial.stc.gsi('e')) = [];
    Trial.stc.addState(Trial.spath,...
                       Trial.filebase,...
                       Errors2,...
                       Trial.xyz.sampleRate,...
                       Trial.sync.copy,...
                       Trial.sync.data(1),...
                       'errors','e');
    Trial.stc{'e'}.save(1);

end




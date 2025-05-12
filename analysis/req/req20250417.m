MjgER2016_load_data();



UnitsNew = {};
Trial = MTATrial.validate('Ed10-20140820.cof.all');

Trial = MTATrial.validate('jg05-20120309.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',true);
Trial = MTATrial.validate('jg05-20120310.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',true);
Trial = MTATrial.validate('jg05-20120311.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',true);
Trial = MTATrial.validate('jg05-20120312.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',true);
Trial = MTATrial.validate('jg05-20120315.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',true);
Trial = MTATrial.validate('jg05-20120316.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',true);
Trial = MTATrial.validate('jg05-20120317.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',true);


s = MTASession.validate('Ed10-20140815.cof.all');
s.spk.create(s);
s.save();
Trial = MTATrial.validate('Ed10-20140815.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trial,'overwrite',true);
Trials{27} = Trial;

s = MTASession.validate('Ed10-20140816.cof.all');
s.spk.create(s);
s.save();
Trial = MTATrial.validate('Ed10-20140816.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trial,'overwrite',true);
Trials{6} = Trial;

s = MTASession.validate('Ed10-20140817.cof.all');
s.spk.create(s);
s.save();
Trial = MTATrial.validate('Ed10-20140817.cof.gnd');
pft = pfs_2d_theta(Trial,'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trial,'overwrite',true);
Trials{7} = Trial;

nq = NeuronQuality(Trial, [], 'overwrite', true);

s = MTASession.validate('jg05-20120311.cof.all');
s.spk.create(s);
s.save();
Trial = MTATrial.validate('jg05-20120311.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trial,'overwrite',true);
Trials{19} = Trial;


s = MTASession.validate('jg05-20120312.cof.all');
s.spk.create(s);
s.save();
Trial = MTATrial.validate('jg05-20120312.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trial,'overwrite',true);
Trials{20} = Trial;

s = MTASession.validate('jg05-20120309.cof.all');
s.spk.create(s);
s.save();
Trial = MTATrial.validate('jg05-20120309.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trial,'overwrite',true);
Trials{17} = Trial;


Trial = MTATrial.validate('jg05-20120315.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',false);
Trial = MTATrial.validate('jg05-20120316.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',false);


s = MTASession.validate('er01-20110719.cof.all');
s.spk.create(s);
s.save();
Trial = MTATrial.validate('er01-20110719.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trial,'overwrite',true);
Trials{1} = Trial;

s = MTASession.validate('er01-20110721.cof.all');
s.spk.create(s);
s.save();
Trial = MTATrial.validate('er01-20110721.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trial,'overwrite',true);
Trials{2} = Trial;


s = MTASession.validate('er01-20110722.cof.all');
s.spk.create(s);
s.save();
Trial = MTATrial.validate('er01-20110722.cof.all');
pft = pfs_2d_theta(Trial,'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trial,'overwrite',true);
Trials{28} = Trial;


states = {'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta','lpause&theta'};
Pft = cf(@(T) pfs_2d_theta(T,[],'overwrite',false,'purge',false), Trials);
Pfs = cf(@(T) pfs_2d_states(T,'all',[],states,'overwrite',true,'purge',true), Trials);


for tid = 24:33
    report_placefield_summary(Trials{tid})
end


units = Trials{2}.spk.get_unit_set(Trial,'placecells');
figure,
for unit = units(:)'
    plot(Pft{1}, unit, [],'text');
    title(num2str(unit));
    waitforbuttonpress();
end

    

% $$$ figure();
% $$$ clf();
% $$$ ucnt = numel(Trial.spk.map(:,1));
% $$$ for unit = 1:ucnt
% $$$     subplot(ceil(ucnt/15),15,unit);
% $$$     plot(pft,unit,[],'text');
% $$$     title(num2str(unit));
% $$$ end
% $$$ axs = flipud(get(gcf,'Children'));


for tid = 1:numel(Trials)
    disp(tid);
    ucnt = Trials{tid}.spk.map(:,1)';
    Trials{tid}.nq = [];
    Trials{tid}.load('nq');
    nq = Trials{tid}.nq;
    if numel(nq.eDist) ~= ucnt
        NeuronQuality(Trials{tid},'overwrite',true)
        Trials{tid}.nq = [];
        Trials{tid}.load('nq');
        nq = Trials{tid}.nq;
    end
    rmaps = af(@(u) plot(Pft{tid},u,[],'text'), ucnt);
    rmaps = cat(3, rmaps{:});
    mrate = pft{tid}.maxRate(ucnt)';
    mzero = sum((reshape(rmaps,[],size(rmaps,3))./mrate)<0.5,'omitnan');
    srate = sum(reshape(rmaps,[],size(rmaps,3)),'omitnan');
    pyrId = find(mrate>1 & mzero>500 & nq.SpkWidthR'>0.6 );
    UnitsNew{tid} = Trials{tid}.spk.map(pyrId,1);
    units = set_unit_set(Trials{tid}.spk, ...
                         Trials{tid},     ...
                         'placecells',    ...
                         UnitsNew{tid});
end





cf(@(T) report_placefield_summary(T), Trials(2:10));

for pyr = pyrId
    axs(pyr).YColor='red';
    axs(pyr).XColor='red';
end

for pyr = pyrId
    axs(pyr).YColor='black';
    axs(pyr).XColor='black';
end

xyz = Trial.load('xyz');
cmtch = {};


% >>> er01 >>> ------------------------------------------
% >>>    0719 - 0721   >>> ------------------------------

tid1 = 1;
tid2 = 2;

s = MTASession.validate(Trials{tid1}.filebase);
s.spk.create(s);
s.save();
Trials{tid1} = MTATrial.validate(Trials{tid1}.filebase);
pft = pfs_2d_theta(Trials{tid1},'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trials{tid1},'overwrite',true);

s = MTASession.validate(Trials{tid2}.filebase);
s.spk.create(s);
s.save();
Trials{tid2} = MTATrial.validate(Trials{tid2}.filebase);
pfs = pfs_2d_theta(Trials{tid2},'overwrite',false);
pfs.purge_savefile();
pfs = pfs_2d_theta(Trials{tid2},'overwrite',true);


cmtch = {};
shnk = 1;
cmtch{shnk} = [ ...

];
cmtch{shnk} = [0,0];

figure(1)
clf();
ucnt = sum(Trials{tid1}.spk.map(:,2)==shnk);
offset = sum(Trials{tid1}.spk.map(:,2)<shnk);
for k = 1:ucnt
    if ismember(k+offset, cmtch{shnk}(:,1)); continue; end;    
    subplot(ceil(ucnt/10),10,k);
    plot(pft,k+offset,[],'colorbar');
    title(num2str(k+offset));
end
figure(2);
clf();
ucnt = sum(Trials{tid2}.spk.map(:,2)==shnk);
offset = sum(Trials{tid2}.spk.map(:,2)<shnk);
for k = 1:ucnt
    if ismember(k+offset, cmtch{shnk}(:,2)); continue; end;
    subplot(ceil(ucnt/10),10,k);
    plot(pfs,k+offset,[],'colorbar');
    title(num2str(k+offset));
end



Trials{tid1}.spk.map(Trials{tid1}.spk.map(:,2)==shnk,:)
Trials{tid2}.spk.map(Trials{tid2}.spk.map(:,2)==shnk,:)
     1     1     2                 1     1     2
     2     1     3                 2     1     3
     3     1     4                 3     1     4
     4     1     5                 4     1     5
     5     1     6                 5     1     6
     6     1     7                 6     1     7
     7     1     8                 7     1     8
     8     1     9                 8     1     9
     9     1    10                 9     1    10
    10     1    11                10     1    11
    11     1    12                11     1    12
    12     1    13                12     1    13
    13     1    14                13     1    14
    14     1    15                14     1    15
    15     1    16                15     1    16
    16     1    17                16     1    17
    17     1    18                17     1    18
    18     1    19                18     1    19
    19     1    20                19     1    20
    20     1    21                20     1    21
    21     1    22                21     1    22
    22     1    23                22     1    23
    23     1    24                23     1    24
    24     1    25                24     1    25
    25     1    26                25     1    26
    26     1    27                26     1    27
    27     1    28                27     1    28
    28     1    29                28     1    29
    29     1    30
    30     1    31
    31     1    32
    32     1    33


shnk = 2;
cmtch{shnk} = [ ...
    21, 16; ... ???
];
cmtch{shnk} = [0,0];

figure(1)
clf();
ucnt = sum(Trials{tid1}.spk.map(:,2)==shnk);
offset = sum(Trials{tid1}.spk.map(:,2)<shnk);
for k = 1:ucnt
    if ismember(k+offset, cmtch{shnk}(:,1)); continue; end;    
    subplot(ceil(ucnt/10),10,k);
    plot(pft,k+offset,[],'colorbar');
    title(num2str(k+offset));
end
figure(2);
clf();
ucnt = sum(Trials{tid2}.spk.map(:,2)==shnk);
offset = sum(Trials{tid2}.spk.map(:,2)<shnk);
for k = 1:ucnt
    if ismember(k+offset, cmtch{shnk}(:,2)); continue; end;
    subplot(ceil(ucnt/10),10,k);
    plot(pfs,k+offset,[],'colorbar');
    title(num2str(k+offset));
end

Trials{tid1}.spk.map(Trials{tid1}.spk.map(:,2)==shnk,:)
Trials{tid2}.spk.map(Trials{tid2}.spk.map(:,2)==shnk,:)

% <<< 0719 - 0721 <<< -----------------------------------
% >>> 0721 - 0722   >>> ------------------------------
cmtch = {};
shnk = 1; cmtch{shnk} = [ ];
shnk = 2;
cmtch{shnk} = [ ...
    31, 37; ...
    53, 56; ...
];
shnk = 3;
cmtch{shnk} = [ ...
    65, 91; ... p ?
    75, 74; ... P
    57, 92; ... P S
    59, 76; ... I    
];
shnk = 4;
cmtch{shnk} = [ ...
    80, 111; ...
];
shnk = 5; cmtch{shnk} = [ ];
shnk = 6; cmtch{shnk} = [ ];
shnk = 7; cmtch{shnk} = [ ];
shnk = 8; cmtch{shnk} = [ ];



% <<< 0721 - 0722 <<< -----------------------------------
% <<< er01 <<< ------------------------------------------
% >>> ER06 >>> ------------------------------------------
% >>> I 0612 - 0613   >>> -------------------------------
figure,plot(xyz(:,1,1),xyz(:,1,2),'.');

tid1 = 2;
tid2 = 2;

s = MTASession.validate(Trials{tid1}.filebase);
s.spk.create(s);
s.save();
Trials{tid1} = MTATrial.validate(Trials{tid1}.filebase);
pft = pfs_2d_theta(Trials{tid1},'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trials{tid1},'overwrite',true);

s = MTASession.validate(Trials{tid2}.filebase);
s.spk.create(s);
s.save();
Trials{tid2} = MTATrial.validate(Trials{tid2}.filebase);
pfs = pfs_2d_theta(Trials{tid2},[],'thetarc-sit-groom','overwrite',false);
pfs.purge_savefile();
pfs = pfs_2d_theta(Trials{tid2},[],'thetarc-sit-groom','overwrite',true);


cmtch = {};

shnk = 8;
cmtch{shnk} = [0,0];

figure(1)
clf();
ucnt = sum(Trials{tid1}.spk.map(:,2)==shnk);
offset = sum(Trials{tid1}.spk.map(:,2)<shnk);
for k = 1:ucnt
    if ismember(k+offset, cmtch{shnk}(:,1)); continue; end;    
    subplot(ceil(ucnt/10),10,k);
    plot(pft,k+offset,[],'colorbar');
    title(num2str(k+offset));
end
figure(2);
clf();
ucnt = sum(Trials{tid2}.spk.map(:,2)==shnk);
offset = sum(Trials{tid2}.spk.map(:,2)<shnk);
for k = 1:ucnt
    if ismember(k+offset, cmtch{shnk}(:,2)); continue; end;
    subplot(ceil(ucnt/10),10,k);
    plot(pfs,k+offset,[],'colorbar');
    title(num2str(k+offset));
end

Trials{tid1}.spk.map(Trials{tid1}.spk.map(:,2)==shnk,:)
Trials{tid2}.spk.map(Trials{tid2}.spk.map(:,2)==shnk,:)

     1     1     2             1     1     2
     2     1     3             2     1     3
     3     1     4             3     1     4
     4     1     5             4     1     5
     5     1     6             5     1     6
     6     1     7             6     1     7
     7     1     8             7     1     8
     8     1     9             8     1     9
     9     1    10             9     1    10
    10     1    11            10     1    11
    11     1    12            11     1    12
    12     1    13            12     1    13
    13     1    14            13     1    14
    14     1    15            14     1    15
    15     1    16            15     1    16
    16     1    17            16     1    17
    17     1    18            17     1    18
    18     1    19            18     1    19
    19     1    20            19     1    20
    20     1    21            20     1    21
    21     1    22            21     1    22
    22     1    23            22     1    23
    23     1    24            23     1    24
    24     1    25            24     1    25
    25     1    26            25     1    26
    26     1    27            26     1    27
    27     1    28            27     1    28
    28     1    29            28     1    29
    29     1    30            29     1    30
    30     1    31            30     1    31
    31     1    32            31     1    32
    32     1    33            32     1    33
    33     1    34            33     1    34
    34     1    35            34     1    35
    35     1    36            35     1    36
    36     1    37            36     1    37
    37     1    38            37     1    38
    38     1    39            38     1    39
    39     1    40            39     1    40
    40     1    41            40     1    41
    41     1    42            41     1    42
    42     1    43            42     1    43
    43     1    44            43     1    44
    44     1    45            44     1    45
                              45     1    46
                              46     1    47
                              47     1    48
                              48     1    49
% <<< 0612 - 0613 <<< -----------------------------------
% >>> I 0624 - 0625   >>> -------------------------------
tid1 = 2;
tid2 = 2;

s = MTASession.validate(Trials{tid1}.filebase);
s.spk.create(s);
s.save();
Trials{tid1} = MTATrial.validate(Trials{tid1}.filebase);
pft = pfs_2d_theta(Trials{tid1},'overwrite',false);
pft.purge_savefile();
pft = pfs_2d_theta(Trials{tid1},'overwrite',true);

s = MTASession.validate(Trials{tid2}.filebase);
s.spk.create(s);
s.save();
Trials{tid2} = MTATrial.validate(Trials{tid2}.filebase);
pfs = pfs_2d_theta(Trials{tid2},[],'thetarc-sit-groom','overwrite',false);
pfs.purge_savefile();
pfs = pfs_2d_theta(Trials{tid2},[],'thetarc-sit-groom','overwrite',true);


cmtch = {};

shnk = 8;
cmtch{shnk} = [0,0];

figure(1)
clf();
ucnt = sum(Trials{tid1}.spk.map(:,2)==shnk);
offset = sum(Trials{tid1}.spk.map(:,2)<shnk);
for k = 1:ucnt
    if ismember(k+offset, cmtch{shnk}(:,1)); continue; end;    
    subplot(ceil(ucnt/10),10,k);
    plot(pft,k+offset,[],'colorbar');
    title(num2str(k+offset));
end
figure(2);
clf();
ucnt = sum(Trials{tid2}.spk.map(:,2)==shnk);
offset = sum(Trials{tid2}.spk.map(:,2)<shnk);
for k = 1:ucnt
    if ismember(k+offset, cmtch{shnk}(:,2)); continue; end;
    subplot(ceil(ucnt/10),10,k);
    plot(pfs,k+offset,[],'colorbar');
    title(num2str(k+offset));
end

Trials{tid1}.spk.map(Trials{tid1}.spk.map(:,2)==shnk,:)
Trials{tid2}.spk.map(Trials{tid2}.spk.map(:,2)==shnk,:)











% <<< 0624 - 0625 <<< -----------------------------------
% <<< ER06 <<< ------------------------------------------
% >>> Ed10 >>> ------------------------------------------
% >>> I 0813 - 0814   >>> -------------------------------
% <<<   0813 - 0814   <<< -------------------------------
% >>> I 0814 - 0815   >>> -------------------------------
% <<<   0814 - 0815   <<< -------------------------------
% >>>   0815 - 0816   >>> -------------------------------
cmtch = {};
shnk = 1;
cmtch{shnk} = [ 9,  5];
shnk = 2;
cmtch{shnk} = [ 16, 35; ...
                17, 29; ...
                20, 30];
shnk = 2;
cmtch{shnk} = [ 16, 35; ...
                17, 29; ...
                20,  30];
shnk = 3;
cmtch{shnk} = [ 22, 52; ...
                25, 47; ...
                26, 53; ...
                29, 48; ... ???
                29, 41; ... ???
                32, 66; ...
                40, 64; ... ???
                51, 67; ... ???
              ];
shnk = 4;
cmtch{shnk} = [ 63, 82; ...
                81, 93; ... ???
                92, 89; ... ???
              ]; 
% <<< 0815 - 0816 <<< -----------------------------------
% >>>   0816 - 0817   >>> -------------------------------
cmtch = {};
shnk = 1; cmtch{shnk} = [];
shnk = 2;
cmtch{shnk} = [ 28; 28; ... S
                30, 31; ... I
                36; 32; ... S
                50, 39];
shnk = 3; 
cmtch{shnk} = [ 41, 51; ... 
                45, 41; ...
                49, 46; ...
                48, 47; ... 
                52, 52; ...
                53, 48; ...
                61, 71; ...
                62, 54; ...
                63, 50; ... ???
                64, 46; ...
                65, 70; ...
                66, 63; ...
                67, 69; ...
                75, 60; ...
                76, 61; ... 
                ];
shnk = 4;
cmtch{shnk} = [  ];
% <<< 0816 - 0817 <<<------------------------------------
% <<< Ed10 <<< ------------------------------------------
% >>> jg05 >>> ------------------------------------------
% >>>   0309 - 0310   >>> -------------------------------
tid = 17;
cmtch{tid} = {};

shnk = 1;  cmtch{tid}{shnk} = [];
shnk = 2;  cmtch{tid}{shnk} = [];
shnk = 3;  cmtch{tid}{shnk} = [];
shnk = 4;  cmtch{tid}{shnk} = [];
shnk = 5;  cmtch{tid}{shnk} = [];

shnk = 6;  
cmtch{tid}{shnk} = [ ...
    8, 18;
    10,  7; ...
    11,  5; ...
    15,  6; ...
    20, 11; ...
];


shnk = 7; 
cmtch{tid}{shnk} = [ ...
    26, 27; ...
    27, 28; ...
    30, 44; ... ???
    34, 25; ...
    33, 43; ...                
    35, 40; ...
    43, 39; ...
    48, 36; ...
    51, 23; ...
    52, 21; ... ???
    55, 34; ...
];

shnk = 8; 
cmtch{tid}{shnk} = [ ...
    59, 49; ...
    61, 60; ...
    62, 50; ...
    64, 59; ...
    66, 71; ...
    68, 47; ...
    69, 56; ...
    70, 74; ... ???
    73, 72; ... ???
    74, 57; ...
    78, 52; ...
    80, 75; ...
    82, 80; ...
    85, 84; ...
    87, 46; ...
    89, 66; ...
    91, 64; ...                
];


% <<< 0309 - 0310 <<< -----------------------------------
% >>>   0310 - 0311   >>> -------------------------------

tid = 18;
cmtch{tid} = {};

shnk = 1;  cmtch{tid}{shnk} = [];
shnk = 2;  cmtch{tid}{shnk} = [];
shnk = 3;  cmtch{tid}{shnk} = [];
shnk = 4;  cmtch{tid}{shnk} = [];
shnk = 5;  cmtch{tid}{shnk} = [];

shnk = 6;  
cmtch{tid}{shnk} = [ ...
    5, 69; ...
    6, 68; ...
    7, 70; ...
    9, 74; ...
    10, 71; ...
    12, 87; ...
    14, 78; ...
    16, 89; ...
    17, 88; ...
    18, 83; ...
    19, 72; ...
];

shnk = 7; 
cmtch{tid}{shnk} = [ ...
    21, 112;...
    28, 107; ...
    23, 103; ...
    24, 123; ...
    33, 124; ...
    43, 116; ...
];

shnk = 8; 
cmtch{tid}{shnk} = [ ...
    46,148; ...
    49,133; ...
    50,141; ...
    51,135; ...
    56,139; ...
    57,150; ...
    59,132; ...
    60,134; ...
    61,142; ...
    66,144; ...
    71,130; ...
    74,149; ...
    76,136; ...
    79,143; ...
    82,146  ...
];

shnk = 9;  cmtch{tid}{shnk} = [];
shnk = 10;  cmtch{tid}{shnk} = [];
shnk = 11;  cmtch{tid}{shnk} = [];
shnk = 12;  cmtch{tid}{shnk} = [];
% <<< 0310 - 0311 <<< -----------------------------------
% >>>   0310 - 0312   >>> -------------------------------

s = MTASession.validate('jg05-20120310.cof.all');
s.spk.create(s);
s.save();
Trials{18} = MTATrial.validate('jg05-20120310.cof.all');
%pft.purge_savefile();
pft = pfs_2d_theta(Trials{18},'overwrite',true);

s = MTASession.validate('jg05-20120312.cof.all');
s.spk.create(s);
s.save();
Trials{20} = MTATrial.validate('jg05-20120312.cof.all');
pfs.purge_savefile();
pfs = pfs_2d_theta(Trials{20},'overwrite',true);

cmtch = {};

shnk = 1;  cmtch{shnk} = [];
shnk = 2;  cmtch{shnk} = [];
shnk = 3;  cmtch{shnk} = [];
shnk = 4;  cmtch{shnk} = [];
shnk = 5;  cmtch{shnk} = [];
shnk = 5;  cmtch{shnk} = [];
shnk = 6;
cmtch{shnk} = [15, 42; ...
               17, 57;...
               16, 61;];
shnk = 7;  cmtch{shnk} = [0,0];
cmtch{shnk} = [ 27, 74; ...
                28, 75; ...
                38, 91; ...
                39, 81; ...
                42, 83; ...
                43, 90; ...
              ];

shnk = 8;  cmtch{shnk} = [0,0];
cmtch{shnk} = [ 49, 128; ...
                52, 137; ...
                54, 129; ...
                56, 109; ... ???
                57, 124; ...
                59, 117; ...
                60, 115; ... ???
                61, 102; ...
                63, 111; ... ???
                71, 104; ...
                81, 106; ...
                84, 142; ... ???
              ];

umerge0312 = [120,122];
rclust0312 = {[118,132], ...
              [8, 19, 33];

% <<< 0310 <-> 0312 <<< ---------------------------------
% >>>   0311 - 0312   >>> -------------------------------
tid = 19;
cmtch{tid} = {};

shnk = 1;  cmtch{tid}{shnk} = [];
shnk = 2;  cmtch{tid}{shnk} = [];
shnk = 3;  cmtch{tid}{shnk} = [];
shnk = 4;
cmtch{tid}{shnk} = [ ...
    7,   7; ...
    14,  8; ...
    29, 12 ...
];

shnk = 5;
cmtch{tid}{shnk} = [ ...
    64, 23; ...
    65, 29; ...
    49, 18; ...
    46, 33; ...
    33, 17; ...
    67, 19; ...
    60, 32 ...
];
shnk = 6;
cmtch{tid}{shnk} = [ ...
    68, 41; ...
    69, 43; ...
    70, 48; ...
    88, 57; ...
    89, 61; ...
    90, 39; ...
    74, 50; ...
    73, 36; ...
    72, 42; ...
    93, 59; ...
    93, 58; ...
    86, 55; ...
    94, 35 ...
];

shnk = 7;
cmtch{tid}{shnk} = [ ...
    97, 87 ...
];

shnk = 8;
cmtch{tid}{shnk} = [ ...
    129, 130;...  
    133, 128;...
    134, 115;...
    137, 106;...
    138, 144;...
    142, 102;...
    139, 114;...
    140, 136;...
    141, 105;...
    143, 118;...
    144, 145;...
    146, 122;...
    147, 132;...
    148, 124;...
    149, 103;...
    150, 124 ... ???
];

shnk = 9;
cmtch{tid}{shnk} = [];
shnk = 10;
cmtch{tid}{shnk} = [ ...
    169, 159;...
    173, 171 ...
];
shnk = 11;
cmtch{tid}{shnk} = [ ...
    178, 176;...
    179, 177;...
    180, 178 ...
];
shnk = 12;
cmtch{tid}{shnk} = [];

% <<< 0311 <-> 0312 <<< ---------------------------------
% >>>   0312 - 0315   >>> -------------------------------

shnk = 1; cmtch{shnk} = [];
shnk = 2; cmtch{shnk} = [];
shnk = 3; cmtch{shnk} = [];
shnk = 4; cmtch{shnk} = [];
shnk = 5; cmtch{shnk} = [];
shnk = 6; cmtch{shnk} = [];

shnk = 7;
cmtch{shnk} = [0,0];
cmtch{shnk} = [ 70, 35; ...
                75, 34; ...
                76, 31; ... M
                92, 33; ... M
                95, 33; ... M
                98, 27; ... M
                96, 25 ... M
              ];
shnk = 8;
cmtch{shnk} = [0,0];
cmtch{shnk} = [101, 45;...
               106, 63;...
               111, 75;...
               128, 77;... M
               134, 65;...
              ];

% <<< 0312 - 0315 <<< -----------------------------------
% >>>   0315 - 0316   >>> -------------------------------

cmtch{21} = {};
cmtch{21}{6} = [];
cmtch{21}{6}(end+1,:) = [ 3, 7];% GUESS PYR
cmtch{21}{6}(end+1,:) = [ 4,14];% MATCH PYR
cmtch{21}{6}(end+1,:) = [ 5,17];% MATCH INT
cmtch{21}{6}(end+1,:) = [ 6,13];% MATCH PYR good
cmtch{21}{6}(end+1,:) = [ 9,15];% GUESS PYR 10-EMPTY
cmtch{21}{6}(end+1,:) = [12,19];% MATCH PYR LLL
cmtch{21}{6}(end+1,:) = [22,16];% Partial MATCH PYR
cmtch{21}{6}(end+1,:) = [23, 6];% GUESS PYR LowRate

cmtch{21}{7} = [];
cmtch{21}{7}(end+1,:) = [24,42];% GUESS PYR good Mid Loc-lloc
cmtch{21}{7}(end+1,:) = [26,31];% GUESS PYR EMPTY
cmtch{21}{7}(end+1,:) = [27,41];% MATCH PYR Good??
cmtch{21}{7}(end+1,:) = [29,29];% MATCH INT
cmtch{21}{7}(end+1,:) = [30,27];% MATCH INT
cmtch{21}{7}(end+1,:) = [32,38];% MATCH PYR Good Edg Loc-All
cmtch{21}{7}(end+1,:) = [33,30];% MATCH PYR Good Mid Loc-Loc
cmtch{21}{7}(end+1,:) = [34,28];% MATCH INT
cmtch{21}{7}(end+1,:) = [40,23];% GUESS PYR Good Mid Rear-Rear
cmtch{21}{7}(end+1,:) = [42,36];% GUESS PYR Prob Edg lloc-lloc

cmtch{21}{8} = [];
cmtch{21}{8}(end+1,:) = [45,56];% GUESS PYR
cmtch{21}{8}(end+1,:) = [48,59];% MATCH PYR good rear
cmtch{21}{8}(end+1,:) = [51,46];% GUESS PYR unkn
cmtch{21}{8}(end+1,:) = [53,62];% GUESS PYRunkn
cmtch{21}{8}(end+1,:) = [57,64];% GUESS PYR unkn
cmtch{21}{8}(end+1,:) = [58,49];% MATCH INT
cmtch{21}{8}(end+1,:) = [59,60];% GUESS PYR
cmtch{21}{8}(end+1,:) = [60,55];% MATCH INT
cmtch{21}{8}(end+1,:) = [61,48];% MATCH PYR good loc loc
cmtch{21}{8}(end+1,:) = [62,54];% GUESS PYR
cmtch{21}{8}(end+1,:) = [63,65];% MATCH PYR good
cmtch{21}{8}(end+1,:) = [65,47];% GUESS PYR remapping or no match
cmtch{21}{8}(end+1,:) = [70,53];% MATCH PYR unkn S S ?
cmtch{21}{8}(end+1,:) = [72,57];% MATCH PYR unkn S R
cmtch{21}{8}(end+1,:) = [73,51];% MATCH PYR good R R
cmtch{21}{8}(end+1,:) = [74,50];% GUESS PYR good S L
cmtch{21}{8}(end+1,:) = [75,58];% MATCH PYR unkn S S??
cmtch{21}{8}(end+1,:) = [77,61];% GUESS PYR good L L


% <<< 0315 - 0316 <<< -----------------------------------
% >>>   0316 - 0317   >>> -------------------------------


cmtch{22} = {}
cmtch{22}{6} = [];
cmtch{22}{6}(end+1,:) = [3,11]; % Maybe
%cmtch{22}{6}(end+1,:) = [4,11]; % Maybe
cmtch{22}{6}(end+1,:) = [5,10]; % Probably
cmtch{22}{6}(end+1,:) = [12,12];% Maybe PYR
%cmtch{22}{6}(end+1,:) = [12,20];% Maybe PYR
cmtch{22}{6}(end+1,:) = [13,29];% MATCH PYR
cmtch{22}{6}(end+1,:) = [18,26];% MATCH PYR
cmtch{22}{6}(end+1,:) = [19,22];% MATCH PYR
%cmtch{22}{6}(end+1,:) = [19,8]; % maybe PYR
cmtch{22}{6}(end+1,:) = [21,18];% MATCH PYR
cmtch{22}{6}(end+1,:) = [22,31];% MATCH PYR
%cmtch{22}{6}(end+1,:) = [22,30];% MATCH PYR

cmtch{22}{7} = [];
cmtch{22}{7}(end+1,:) = [23,58];% MATCH PYR
cmtch{22}{7}(end+1,:) = [27,49];% maybe INT 
cmtch{22}{7}(end+1,:) = [28,52];% MATCH INT
cmtch{22}{7}(end+1,:) = [29,44];% MATCH INT
cmtch{22}{7}(end+1,:) = [30,54];% MATCH PYR  56 xcorr 41
cmtch{22}{7}(end+1,:) = [32,60];% maybe PYR sleep/immobile?
cmtch{22}{7}(end+1,:) = [38,48];% probably PYR
cmtch{22}{7}(end+1,:) = [40,42];% maybe PYR sleep/immobile?
cmtch{22}{7}(end+1,:) = [41,50];% MATCH PYR 
cmtch{22}{7}(end+1,:) = [42,51];% MATCH PYR 
cmtch{22}{7}(end+1,:) = [43,36];% MATCH PYR

cmtch{22}{8} = [];
cmtch{22}{8}(end+1,:) = [47,66];% MATCH PYR  
cmtch{22}{8}(end+1,:) = [48,72];% MATCH PYR  T
cmtch{22}{8}(end+1,:) = [50,70];% MATCH PYR  T
cmtch{22}{8}(end+1,:) = [51,67];% MATCH PYR  T
cmtch{22}{8}(end+1,:) = [59,71];% MATCH PYR  T
cmtch{22}{8}(end+1,:) = [65,69];% MATCH PYR  T
cmtch{22}{8}(end+1,:) = [61,63];% MATCH PYR? T
cmtch{22}{8}(end+1,:) = [64,75];% EMPTY PYR?

% <<< 0316 - 0317 <<< -----------------------------------
% >>> I 0324 - 0325   >>> -------------------------------
% <<< I 0324 - 0325   <<< -------------------------------
% <<< jg05 <<< ------------------------------------------

configure_default_args();

MjgER2016_figure_BhvPlacefields_args();


bfs = cf(                                             ...
    @(t,u)                                            ...
    compute_bhv_ratemaps( t, u, 'overwrite', true),   ...
    Trials,                                           ...
    UnitsNew([6,7,17,18,19,20])                       ...
);

session_name = 'ER06-20130625';
data_paths = struct( ...
    'xyz', '/storage/gravio/data/processed/xyz/ER06', ...
    'ephys', '/storage/evgeny/data/processed/ER06'    ...
);
link_session( session_name, data_paths);




figure,
for unit = units(:)'
for sts = 1:numel(Pfs{1})
    subplot(1,numel(Pfs{1}),sts)
    plot(Pfs{1}{sts},unit,[],'text');
end
waitforbuttonpress();
end
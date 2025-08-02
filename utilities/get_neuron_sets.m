function neuron_sets = get_neuron_sets(Trials)


%sessionList = get_session_list_v3('BehaviorPlaceCode');
%Trials = af(@(s) MTATrial.validate(s), sessionList);

cmtch = cell([numel(Trials), numel(Trials)]);
% >>> er01 >>> ------------------------------------------
% >>>  E 0719 - 0721   >>> ------------------------------

% $$$ tid1 = 1;
% $$$ tid2 = 2;
% $$$ 
% $$$ s = MTASession.validate(Trials{tid1}.filebase);
% $$$ s.spk.create(s);
% $$$ s.save();
% $$$ Trials{tid1} = MTATrial.validate(Trials{tid1}.filebase);
% $$$ pft = pfs_2d_theta(Trials{tid1},'overwrite',false);
% $$$ pft.purge_savefile();
% $$$ pft = pfs_2d_theta(Trials{tid1},'overwrite',true);
% $$$ 
% $$$ s = MTASession.validate(Trials{tid2}.filebase);
% $$$ s.spk.create(s);
% $$$ s.save();
% $$$ Trials{tid2} = MTATrial.validate(Trials{tid2}.filebase);
% $$$ pfs = pfs_2d_theta(Trials{tid2},'overwrite',false);
% $$$ pfs.purge_savefile();
% $$$ pfs = pfs_2d_theta(Trials{tid2},'overwrite',true);
% $$$ 
% $$$ tind1 = find_trial_index(Trials, 'er01-20110719.cof.all');
% $$$ tind2 = find_trial_index(Trials, 'er01-20110721.cof.all');
% $$$ cmtch{tind1,tind2} = {};
% $$$ 
% $$$ 
% $$$ figure(1)
% $$$ clf();
% $$$ ucnt = sum(Trials{tid1}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid1}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,1)); continue; end;    
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pft,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ figure(2);
% $$$ clf();
% $$$ ucnt = sum(Trials{tid2}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid2}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,2)); continue; end;
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pfs,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ Trials{tid1}.spk.map(Trials{tid1}.spk.map(:,2)==shnk,:)
% $$$ Trials{tid2}.spk.map(Trials{tid2}.spk.map(:,2)==shnk,:)
% $$$      1     1     2                 1     1     2
% $$$      2     1     3                 2     1     3
% $$$      3     1     4                 3     1     4
% $$$      4     1     5                 4     1     5
% $$$      5     1     6                 5     1     6
% $$$      6     1     7                 6     1     7
% $$$      7     1     8                 7     1     8
% $$$      8     1     9                 8     1     9
% $$$      9     1    10                 9     1    10
% $$$     10     1    11                10     1    11
% $$$     11     1    12                11     1    12
% $$$     12     1    13                12     1    13
% $$$     13     1    14                13     1    14
% $$$     14     1    15                14     1    15
% $$$     15     1    16                15     1    16
% $$$     16     1    17                16     1    17
% $$$     17     1    18                17     1    18
% $$$     18     1    19                18     1    19
% $$$     19     1    20                19     1    20
% $$$     20     1    21                20     1    21
% $$$     21     1    22                21     1    22
% $$$     22     1    23                22     1    23
% $$$     23     1    24                23     1    24
% $$$     24     1    25                24     1    25
% $$$     25     1    26                25     1    26
% $$$     26     1    27                26     1    27
% $$$     27     1    28                27     1    28
% $$$     28     1    29                28     1    29
% $$$     29     1    30
% $$$     30     1    31
% $$$     31     1    32
% $$$     32     1    33
% $$$ 
% $$$ 
% $$$ shnk = 2;
% $$$ cmtch{shnk} = [ ...
% $$$     21, 16; ... ???
% $$$ ];
% $$$ cmtch{shnk} = [0,0];
% $$$ 
% $$$ figure(1)
% $$$ clf();
% $$$ ucnt = sum(Trials{tid1}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid1}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,1)); continue; end;    
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pft,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ figure(2);
% $$$ clf();
% $$$ ucnt = sum(Trials{tid2}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid2}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,2)); continue; end;
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pfs,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ 
% $$$ Trials{tid1}.spk.map(Trials{tid1}.spk.map(:,2)==shnk,:)
% $$$ Trials{tid2}.spk.map(Trials{tid2}.spk.map(:,2)==shnk,:)

% <<<    0719 - 0721   <<< ------------------------------
% >>>    0721 - 0722   >>> ------------------------------
tind1 = find_trial_index(Trials, 'er01-20110721.cof.all');
tind2 = find_trial_index(Trials, 'er01-20110722.cof.all');
cmtch{tind1,tind2} = {};
cmtch{tind1,tind2} = [ ...
    31,  37, 1, 1; ...
    53,  56, 1, 1; ...
    65,  91, 1, 0; ... p ?
    75,  74, 1, 1; ... P
    57,  92, 1, 0; ... P S
    59,  76, 0, 1; ... I    
    80,  111, 1, 1; ...
];
% <<< 0721 - 0722 <<< -----------------------------------
% <<< er01 <<< ------------------------------------------
% >>> ER06 >>> ------------------------------------------
% >>> I 0612 - 0613   >>> -------------------------------

% $$$ figure,plot(xyz(:,1,1),xyz(:,1,2),'.');
% $$$ 
% $$$ tid1 = 2;
% $$$ tid2 = 2;
% $$$ 
% $$$ s = MTASession.validate(Trials{tid1}.filebase);
% $$$ s.spk.create(s);
% $$$ s.save();
% $$$ Trials{tid1} = MTATrial.validate(Trials{tid1}.filebase);
% $$$ pft = pfs_2d_theta(Trials{tid1},'overwrite',false);
% $$$ pft.purge_savefile();
% $$$ pft = pfs_2d_theta(Trials{tid1},'overwrite',true);
% $$$ 
% $$$ s = MTASession.validate(Trials{tid2}.filebase);
% $$$ s.spk.create(s);
% $$$ s.save();
% $$$ Trials{tid2} = MTATrial.validate(Trials{tid2}.filebase);
% $$$ pfs = pfs_2d_theta(Trials{tid2},[],'thetarc-sit-groom','overwrite',false);
% $$$ pfs.purge_savefile();
% $$$ pfs = pfs_2d_theta(Trials{tid2},[],'thetarc-sit-groom','overwrite',true);
% $$$ 
% $$$ 
% $$$ cmtch = {};
% $$$ 
% $$$ shnk = 8;
% $$$ cmtch{shnk} = [0,0];
% $$$ 
% $$$ figure(1)
% $$$ clf();
% $$$ ucnt = sum(Trials{tid1}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid1}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,1)); continue; end;    
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pft,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ figure(2);
% $$$ clf();
% $$$ ucnt = sum(Trials{tid2}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid2}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,2)); continue; end;
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pfs,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ 
% $$$ Trials{tid1}.spk.map(Trials{tid1}.spk.map(:,2)==shnk,:)
% $$$ Trials{tid2}.spk.map(Trials{tid2}.spk.map(:,2)==shnk,:)
% $$$ 
% $$$      1     1     2             1     1     2
% $$$      2     1     3             2     1     3
% $$$      3     1     4             3     1     4
% $$$      4     1     5             4     1     5
% $$$      5     1     6             5     1     6
% $$$      6     1     7             6     1     7
% $$$      7     1     8             7     1     8
% $$$      8     1     9             8     1     9
% $$$      9     1    10             9     1    10
% $$$     10     1    11            10     1    11
% $$$     11     1    12            11     1    12
% $$$     12     1    13            12     1    13
% $$$     13     1    14            13     1    14
% $$$     14     1    15            14     1    15
% $$$     15     1    16            15     1    16
% $$$     16     1    17            16     1    17
% $$$     17     1    18            17     1    18
% $$$     18     1    19            18     1    19
% $$$     19     1    20            19     1    20
% $$$     20     1    21            20     1    21
% $$$     21     1    22            21     1    22
% $$$     22     1    23            22     1    23
% $$$     23     1    24            23     1    24
% $$$     24     1    25            24     1    25
% $$$     25     1    26            25     1    26
% $$$     26     1    27            26     1    27
% $$$     27     1    28            27     1    28
% $$$     28     1    29            28     1    29
% $$$     29     1    30            29     1    30
% $$$     30     1    31            30     1    31
% $$$     31     1    32            31     1    32
% $$$     32     1    33            32     1    33
% $$$     33     1    34            33     1    34
% $$$     34     1    35            34     1    35
% $$$     35     1    36            35     1    36
% $$$     36     1    37            36     1    37
% $$$     37     1    38            37     1    38
% $$$     38     1    39            38     1    39
% $$$     39     1    40            39     1    40
% $$$     40     1    41            40     1    41
% $$$     41     1    42            41     1    42
% $$$     42     1    43            42     1    43
% $$$     43     1    44            43     1    44
% $$$     44     1    45            44     1    45
% $$$                               45     1    46
% $$$                               46     1    47
% $$$                               47     1    48
% $$$                               48     1    49

% <<< I 0612 - 0613   <<< -------------------------------
% >>> I 0613 - 0614   >>> -------------------------------
% <<<   0613 - 0614   <<< -------------------------------
% >>> I 0624 - 0625   >>> -------------------------------

% $$$ tid1 = 2;
% $$$ tid2 = 2;
% $$$ 
% $$$ s = MTASession.validate(Trials{tid1}.filebase);
% $$$ s.spk.create(s);
% $$$ s.save();
% $$$ Trials{tid1} = MTATrial.validate(Trials{tid1}.filebase);
% $$$ pft = pfs_2d_theta(Trials{tid1},'overwrite',false);
% $$$ pft.purge_savefile();
% $$$ pft = pfs_2d_theta(Trials{tid1},'overwrite',true);
% $$$ 
% $$$ s = MTASession.validate(Trials{tid2}.filebase);
% $$$ s.spk.create(s);
% $$$ s.save();
% $$$ Trials{tid2} = MTATrial.validate(Trials{tid2}.filebase);
% $$$ pfs = pfs_2d_theta(Trials{tid2},[],'thetarc-sit-groom','overwrite',false);
% $$$ pfs.purge_savefile();
% $$$ pfs = pfs_2d_theta(Trials{tid2},[],'thetarc-sit-groom','overwrite',true);
% $$$ 
% $$$ 
% $$$ cmtch = {};
% $$$ 
% $$$ shnk = 8;
% $$$ cmtch{shnk} = [0,0];
% $$$ 
% $$$ figure(1)
% $$$ clf();
% $$$ ucnt = sum(Trials{tid1}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid1}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,1)); continue; end;    
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pft,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ figure(2);
% $$$ clf();
% $$$ ucnt = sum(Trials{tid2}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid2}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,2)); continue; end;
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pfs,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ 
% $$$ Trials{tid1}.spk.map(Trials{tid1}.spk.map(:,2)==shnk,:)
% $$$ Trials{tid2}.spk.map(Trials{tid2}.spk.map(:,2)==shnk,:)

% <<< 0624 - 0625 <<< -----------------------------------
% <<< ER06 <<< ------------------------------------------
% >>> Ed10 >>> ------------------------------------------
% >>>   0813 - 0814   >>> -------------------------------
tind1 = find_trial_index(Trials, 'Ed10-20140813.cof.all');
tind2 = find_trial_index(Trials, 'Ed10-20140814.cof.all');

cmtch{tind1,tind2} = [ ...
    10, 40, 0, 1; ... Int
    11, 29, 0, 1; ... Int
    13, 63, 0, 1; ... Int
    24, 54, 1, 1; ... Pyr Slt
    25, 44, 1, 1; ... Pyr Slt
    33, 43, 1, 1; ... Pyr
    40, 70, 1, 0; ... Pyr ???
    41, 39, 1, 0; ... Pyr ???
];
% <<<   0813 - 0814   <<< -------------------------------
% >>>   0814 - 0815   >>> -------------------------------

tind1 = find_trial_index(Trials, 'Ed10-20140814.cof.all');
tind2 = find_trial_index(Trials, 'Ed10-20140815.cof.all');

cmtch{tind1,tind2} = [ ...
     4,  9, 1, 1; ... Pyr
     7, 21, 1, 1; ... Pyr
    10, 20, 0, 1; ... Int
    11, 25, 0, 1; ... Int
    16, 29, 0, 1; ... Int
    29, 26, 0, 1; ... Int
    30, 38, 0, 1; ... Int
    39, 46, 1, 1; ... Pyr
    43, 50, 1, 1; ... Pyr
    44, 40, 1, 1; ... Pyr
    50, 31, 0, 1; ... Int
    54, 48, 1, 1; ... Pyr
    63, 30, 0, 1; ... Int    
    70, 37, 1, 1; ... Pyr
    75, 47, 1, 1; ... Pyr
    77, 77, 1, 1; ... Pyr
    85, 63, 0, 1; ... Int
    86, 73, 1, 0; ... Pyr ???
    89, 81, 1, 1; ... Pyr
    94, 66, 0, 1; ... Int
    99, 60, 0, 1; ... Int
   103, 92, 1, 1; ... Pyr
   105, 57, 1, 1; ... Pyr
   107, 81, 1, 0; ... Pyr Slt
   108, 96, 1, 1; ... Pyr
   109, 68, 1, 0; ... Pyr Slt
   112, 86, 1, 1; ... Pyr
   113, 75, 1, 1; ... Pyr
   114, 74, 1, 0; ... Pyr Slt
   118, 92, 1, 0; ... Pyr Slt
   121, 88, 1, 1; ... Pyr
   123, 64, 0, 0; ... Int BAD
   124, 94, 1, 1; ... Pyr    
];
% merge 0814-4: 8,12
% merge 0814-4: 24, 29




% <<<   0814 - 0815   <<< -------------------------------
% >>>   0815 - 0816   >>> -------------------------------

tind1 = find_trial_index(Trials, 'Ed10-20140815.cof.all');
tind2 = find_trial_index(Trials, 'Ed10-20140816.cof.all');
cmtch{tind1,tind2} = [ ...
    9,  5, 1, 1; ... Pyr
    16, 35, 0, 1; ... Int
    17, 29, 0, 1; ... Int
    20, 30, 0, 1; ... Int
    23, 50, 0, 1; ... Int
    25, 41, 0, 1; ... Int
    26, 53, 0, 1; ... Int
    29, 47, 0, 1; ... Int
    30, 52, 0, 1; ... Int
    33, 62, 1, 0; ... Pyr Slt
    38, 47, 0, 1; ... Int ???
    32, 66, 1, 1; ... Pyr
    45, 32, 1, 0; ... Pyr ???
    47, 75, 1, 1; ... Pyr
    51, 67, 1, 1; ... Pyr
];

% <<<   0815 - 0816   <<< -------------------------------
% >>>   0816 - 0817   >>> -------------------------------

tind1 = find_trial_index(Trials, 'Ed10-20140816.cof.all');
tind2 = find_trial_index(Trials, 'Ed10-20140817.cof.gnd');

cmtch{ tind1, tind2} = [ ...
    28, 28, 1, 0; ... Pyr Slt
    30, 31, 0, 1; ... Int
    36, 32, 1, 0; ... Pyr Slt
    50, 39, 0, 1; ... Int
    41, 51, 0, 1; ... Int
    45, 41, 1, 1; ... Pyr
    49, 46, 1, 1; ... Pyr
    48, 47, 0, 1; ... Int
    52, 52, 0, 1; ... Int
    53, 48, 0, 0; ... Int
    61, 71, 1, 1; ... PYr
    62, 54, 1, 1; ... Pyr
    63, 50, 1, 0; ... Pyr ???
    64, 46, 1, 1; ... Pyr
    65, 70, 1, 1; ... Pyr
    66, 63, 1, 1; ... Pyr
    67, 69, 1, 1; ... Pyr
    75, 60, 1, 1; ... Pyr
    76, 61, 1, 1; ... Pyr
];
% <<<   0816 - 0817   <<<--------------------------------
% <<< Ed10 <<< ------------------------------------------
% >>> jg05 >>> ------------------------------------------
% >>>   0309 - 0310   >>> -------------------------------

tind1 = find_trial_index(Trials, 'jg05-20120309.cof.all');
tind2 = find_trial_index(Trials, 'jg05-20120310.cof.all');
cmtch{tind1,tind2} = [ ...
     8, 18, 1, 1; ... Pyr
    10,  7, 0, 1; ... Int
    11,  5, 0, 1; ... Int
    15,  6, 0, 1; ... Int
    20, 11, 1, 1; ... Pyr
    26, 27, 0, 1; ... Int
    27, 28, 0, 1; ... Int
    30, 44, 1, 0; ... Pyr ???
    34, 25, 1, 1; ... Pyr
    33, 43, 0, 1; ... Int
    35, 40, 1, 1; ... Pyr
    43, 39, 1, 1; ... Pyr
    48, 36, 1, 1; ... Pyr
    51, 23, 1, 1; ... Pyr
    52, 21, 1, 1; ... Pyr ???
    55, 34, 1, 1; ... Pyr
    59, 49, 1, 1; ... Pyr
    61, 60, 1, 1; ... Pyr
    62, 50, 1, 1; ... Pyr
    64, 59, 0, 1; ... Int
    66, 71, 0, 1; ... Int
    68, 47, 1, 1; ... Pyr
    69, 56, 1, 1; ... Pyr
    70, 74, 1, 1; ... Pyr
    73, 72, 1, 1; ... Pyr
    74, 57, 1, 1; ... Pyr
    78, 52, 1, 1; ... Pyr
    80, 75, 1, 1; ... Pyr
    82, 80, 1, 1; ... Pyr
    85, 84, 1, 1; ... Pyr
    87, 46, 1, 1; ... Pyr
    89, 66, 1, 1; ... Pyr
    91, 64, 1, 1; ... Pyr
];

% <<<   0309 - 0310   <<< -------------------------------
% >>>   0310 - 0311   >>> -------------------------------

tind1 = find_trial_index(Trials, 'jg05-20120310.cof.all');
tind2 = find_trial_index(Trials, 'jg05-20120311.cof.all');
cmtch{tind1,tind2} = [ ...
     5, 69, 0, 1; ... Int
     6, 68, 0, 1; ... Int
     7, 70, 0, 1; ... Int
     9, 74, 1, 1; ... Pyr
    10, 71, 1, 1; ... Pyr
    12, 87, 1, 1; ... Pyr
    14, 78, 1, 1; ... Pyr
    16, 89, 1, 1; ... Pyr
    17, 88, 1, 1; ... Pyr
    18, 83, 1, 1; ... Pyr
    19, 72, 1, 1; ... Pyr
    21, 112, 1, 1; ... Pyr
    28, 107, 0, 1; ... Int
    23, 103, 1, 1; ... Pyr
    24, 123, 1, 0; ... Pyr
    33, 124, 1, 1; ... Pyr
    43, 116, 0, 1; ... Int
    46, 148, 1, 1; ... Pyr
    49, 133, 1, 1; ... Pyr
    50, 141, 1, 1; ... Pyr
    51, 135, 1, 1; ... Pyr
    56, 139, 1, 1; ... Pyr
    57, 150, 1, 1; ... Pyr
    59, 132, 0, 1; ... Int
    60, 134, 1, 1; ... Pyr
    61, 142, 1, 1; ... Pyr
    66, 144, 1, 1; ... Pyr
    71, 130, 0, 1; ... Int
    74, 149, 1, 1; ... Pyr
    76, 136, 1, 1; ... Pyr
    79, 143, 1, 1; ... Pyr
    82, 146, 1, 1;  ... Pyr
];

% <<<   0310 - 0311   <<< -------------------------------
% >>>   0310 - 0312   >>> -------------------------------

tind1 = find_trial_index(Trials, 'jg05-20120310.cof.all');
tind2 = find_trial_index(Trials, 'jg05-20120312.cof.all');

cmtch{tind1,tind2} = [ ...
    15,  42, 1, 1; ... pyr
    17,  57, 1, 1; ... pyr
    16,  61, 1, 1; ... pyr
    27,  74, 0, 1; ... Int
    28,  75, 0, 1; ... Int
    38,  91, 1, 1; ... Pyr
    39,  81, 1, 1; ... Pyr
    42,  83, 1, 1; ... Pyr good
    43,  90, 0, 1; ... Int
    49, 128, 1, 1; ... pyr good
    52, 137, 1, 1; ... pyr good
    54, 129, 1, 1; ...
    56, 109, 1, 0; ... ???
    57, 124, 1, 1; ...
    59, 117, 0, 1; ...
    60, 115, 1, 0; ... ???
    61, 102, 1, 1; ...
    63, 111, 1, 0; ... ??? good?
    71, 104, 0, 1; ...
    81, 106, 1, 1; ...
    84, 142, 1, 0; ... ???
];

% $$$ umerge0312 = [120,122];
% $$$ rclust0312 = {[118,132], ...
% $$$               [8, 19, 33];

% <<<   0310 - 0312   <<< ---------------------------------
% >>>   0311 - 0312   >>> -------------------------------

tind1 = find_trial_index(Trials, 'jg05-20120311.cof.all');
tind2 = find_trial_index(Trials, 'jg05-20120312.cof.all');
cmtch{tind1,tind2} = [ ...
    7,   7, 0, 1; ... Int
    14,  8, 0, 1; ... Int
    29, 12, 1, 1; ... Pyr
    64, 23, 1, 1; ... Pyr
    65, 29, 1, 1; ... Pyr
    49, 18, 1, 1; ... Pyr
    46, 33, 1, 1; ... Pyr
    33, 17, 1, 1; ... Pyr
    67, 19, 1, 1; ... Pyr
    60, 32, 1, 1; ... Pyr
    68, 41, 0, 1; ... Int
    69, 43, 0, 1; ... Int
    70, 48, 0, 1; ... Int
    88, 57, 1, 1; ... Pyr
    89, 61, 1, 1; ... Pyr
    90, 39, 1, 1; ... Pyr
    74, 50, 1, 1; ... Pyr
    73, 36, 1, 1; ... Pyr
    72, 42, 1, 1; ... Pyr
    93, 59, 1, 1; ... Pyr
    93, 58, 1, 1; ... Pyr
    86, 55, 1, 1; ... Pyr
    94, 35, 1, 1; ... Pyr
    97, 87, 1, 1; ... Pyr
    129, 130, 1, 1; ... Pyr  
    133, 128, 1, 1; ... Pyr
    134, 115, 1, 1; ... Pyr
    137, 106, 1, 1; ... Pyr
    138, 144, 1, 1; ... Pyr
    142, 102, 1, 1; ... Pyr
    139, 114, 1, 1; ... Pyr
    140, 136, 1, 1; ... Pyr
    141, 105, 1, 1; ... Pyr
    143, 118, 1, 1; ... Pyr
    144, 145, 1, 1; ... Pyr
    146, 122, 1, 1; ... Pyr
    147, 132, 1, 1; ... Pyr
    148, 124, 1, 1; ... Pyr
    149, 103, 1, 1; ... Pyr
    150, 124, 1, 0; ... Pyr ???
];

% $$$ shnk = 10;
% $$$ cmtch{tid}{shnk} = [ ...
% $$$     169, 159, nan, 1;...
% $$$     173, 171, nan, 1; ...
% $$$ ];
% $$$ shnk = 11;
% $$$ cmtch{tid}{shnk} = [ ...
% $$$     178, 176, nan, 1;...
% $$$     179, 177, nan, 1;...
% $$$     180, 178, nan, 1; ...
% $$$ ];

% <<<   0311 - 0312   <<< -------------------------------
% >>>   0312 - 0315   >>> -------------------------------

% $$$ tind1 = find_trial_index(Trials, 'jg05-20120312.cof.all');
% $$$ tind2 = find_trial_index(Trials, 'jg05-20120315.cof.all');
% $$$ 
% $$$ cmtch{tind1,tind2} = [ ...
% $$$     70, 35, 1, 1; ... Pyr
% $$$     75, 34, 0, 1; ... Int
% $$$     76, 31, 1, 1; ... Pyr M
% $$$     92, 33, 1, 1; ... Pyr M
% $$$     95, 33, 1, 1; ... Pyr M
% $$$     98, 27, 1, 1; ... Pyr M
% $$$     96, 25, 1, 1; ... Pyr M
% $$$     101, 45, 1, 1; ... Pyr
% $$$     106, 63, 1, 1; ... Pyr
% $$$     111, 75, 1, 1; ... Pyr
% $$$ %    128, 77, 1, 1; ... Pyr M
% $$$     134, 65, 1, 1; ... Pyr
% $$$ ];

% <<<   0312 - 0315   <<< -------------------------------
% >>>   0315 - 0316   >>> -------------------------------

tind1 = find_trial_index(Trials, 'jg05-20120315.cof.all');
tind2 = find_trial_index(Trials, 'jg05-20120316.cof.all');
cmtch{tind1,tind2}  = [  ...
    3,  7, 0, 0; ... % PYR GUESS
    4, 14, 1, 1; ... % PYR MATCH
    5, 17, 1, 0; ... % INT MATCH
    6, 13, 1, 1; ... % PYR MATCH
    9, 15, 1, 0; ... % PYR GUESS
   12, 19, 1, 1; ... % PYR MATCH 
   22, 16, 1, 1; ... % PYR MATCH
   23,  6, 1, 0; ... % PYR GUESS
    24, 42, 1, 0; ... % PYR GUESS good Mid Loc-lloc
    26, 31, 1, 0; ... % PYR GUESS PYR EMPTY
    27, 41, 1, 1; ... % PYR MATCH Good??
    29, 29, 0, 1; ... % INT MATCH
    30, 27, 0, 1; ... % INT MATCH
    32, 38, 1, 1; ... % PYR MATCH Good Edg Loc-All
    33, 30, 1, 1; ... % PYR MATCH Good Mid Loc-Loc
    34, 28, 0, 1; ... % INT MATCH
    40, 23, 1, 1; ... % PYR GUESS Good Mid Rear-Rear
    42, 36, 1, 0; ... % PYR GUESS Prob Edg lloc-lloc
    45, 56, 1, 0; ... % PYR GUESS 
    48, 59, 1, 1; ... % PYR MATCH good rear
    51, 46, 1, 0; ... % PYR GUESS unkn
    53, 62, 1, 0; ... % PYR GUESS unkn
    57, 64, 1, 0; ... % PYR GUESS unkn
    58, 49, 0, 1; ... % INT MATCH
    59, 60, 1, 0; ... % PYR GUESS
    60, 55, 0, 1; ... % INT MATCH
    61, 48, 1, 1; ... % PYR MATCH good loc loc
    62, 54, 1, 0; ... % PYR GUESS
    63, 65, 1, 1; ... % PYR MATCH good
    65, 47, 1, 0; ... % PYR GUESS remapping or no match
    70, 53, 1, 0; ... % PYR MATCH Slt - unkn S S ?
    72, 57, 1, 0; ... % PYR MATCH Slt - unkn S R
    73, 51, 1, 1; ... % PYR MATCH good R R
    74, 50, 1, 0; ... % PYR GUESS Slt - good S L
    75, 58, 1, 0; ... % PYR MATCH Slt - unkn S S??
    77, 61, 1, 1; ... % PYR GUESS good L L
];    

% <<<   0315 - 0316   <<< -------------------------------
% >>>   0316 - 0317   >>> -------------------------------

tind1 = find_trial_index(Trials, 'jg05-20120316.cof.all');
tind2 = find_trial_index(Trials, 'jg05-20120317.cof.all');
cmtch{tind1,tind2} = [  ...
     3, 11, nan, nan; ... % ???  Maybe - [4,11]; % Maybe
     5, 10, nan, nan; ... % ??? Probably
    12, 12, 1, nan; ... % Maybe PYR  -[12,20];% Maybe PYR
    13, 29, 1, 1; ... % MATCH PYR
    18, 26, 1, 1; ... % MATCH PYR
    19, 22, 1, nan; ... % MATCH PYR%cmtch{22}{6}(end+1,:) = [19,8]; % maybe PYR
    21, 18, 1, 1;; ... % MATCH PYR
    22, 31, 1, nan; ... % MATCH PYR -[22,30];% MATCH PYR ??
    23, 58, 1, 1; ... % PYR MATCH
    27, 49, 0, 0; ... % INT GUESS 
    28, 52, 0, 1; ... % INT MATCH
    29, 44, 0, 1; ... % INT MATCH
    30, 54, 1, 1; ... % PYR MATCH   56 xcorr 41
    32, 60, 1, 0; ... % PYR GUESS  sleep/immobile?
    38, 48, 1, 0; ... % PYR GUESS
    40, 42, 1, 0; ... % PYR GUESS sleep/immobile?
    41, 50, 1, 1; ... % PYR MATCH  
    42, 51, 1, 1; ... % PYR MATCH 
    43, 36, 1, 1; ... % PYR MATCH
    47, 66, 1, 1; ... % PYR MATCH  
    48, 72, 1, 1; ... % PYR MATCH  T
    50, 70, 1, 1; ... % PYR MATCH  T
    51, 67, 1, 1; ... % PYR MATCH  T
    59, 71, 1, 1; ... % PYR MATCH  T
    65, 69, 1, 1; ... % PYR MATCH  T
    61, 63, 1, 1; ... % PYR MATCH  
    64, 75, 1, 1; ... % PYR MATCH Slt 
];

% <<< 0316 - 0317 <<< -----------------------------------
% >>> I 0323 - 0324   >>> -------------------------------
% $$$ tid1 = find_trial_index(Trials,'jg05-20120323.cof.all');
% $$$ tid2 = find_trial_index(Trials,'jg05-20120324.cof.all');
% $$$ 
% $$$ s = MTASession.validate(Trials{tid1}.filebase);
% $$$ s.spk.create(s);
% $$$ s.save();
% $$$ Trials{tid1} = MTATrial.validate(Trials{tid1}.filebase);
% $$$ pft = pfs_2d_theta(Trials{tid1},'overwrite',false);
% $$$ pft.purge_savefile();
% $$$ pft = pfs_2d_theta(Trials{tid1},'overwrite',true);
% $$$ 
% $$$ s = MTASession.validate(Trials{tid2}.filebase);
% $$$ s.spk.create(s);
% $$$ s.save();
% $$$ Trials{tid2} = MTATrial.validate(Trials{tid2}.filebase);
% $$$ pfs = pfs_2d_theta(Trials{tid2},'overwrite',false);
% $$$ pfs.purge_savefile();
% $$$ pfs = pfs_2d_theta(Trials{tid2},'overwrite',true);
% $$$ 
% $$$ cmtch = {};
% $$$ shnk = 4;
% $$$ cmtch{shnk} = [ ...
% $$$     6,  27; ... Int ???
% $$$    20,  34; ... Int ???
% $$$ ];
% $$$ figure(1)
% $$$ clf();
% $$$ ucnt = sum(Trials{tid1}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid1}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,1)); continue; end;    
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pft,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ figure(2);
% $$$ clf();
% $$$ ucnt = sum(Trials{tid2}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid2}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,2)); continue; end;
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pfs,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ 
% $$$ Trials{tid1}.spk.map(Trials{tid1}.spk.map(:,2)==shnk,:)
% $$$ Trials{tid2}.spk.map(Trials{tid2}.spk.map(:,2)==shnk,:)
% $$$ 
% $$$      5     4     2               23     4     2
% $$$      6     4     3               24     4     3
% $$$      7     4     4               25     4     4
% $$$      8     4     5               26     4     5
% $$$      9     4     6               27     4     6
% $$$     10     4     7               28     4     7
% $$$     11     4     8               29     4     8
% $$$     12     4     9               30     4     9
% $$$     13     4    10               31     4    10
% $$$     14     4    11               32     4    11
% $$$     15     4    12               33     4    12
% $$$     16     4    13               34     4    13
% $$$     17     4    14               35     4    14
% $$$     18     4    15               36     4    15
% $$$     19     4    16               37     4    16
% $$$     20     4    17               38     4    17
% $$$     21     4    18               39     4    18
% $$$     22     4    19               40     4    19
% $$$     23     4    20               41     4    20
% $$$     24     4    21               42     4    21
% $$$     25     4    22               43     4    22
% $$$     26     4    23               44     4    23
% $$$     27     4    24               45     4    24
% $$$     28     4    25               46     4    25
% $$$     29     4    26               47     4    26
% $$$                                  48     4    27



% <<<   0323 - 0324   <<< -------------------------------
% >>> I 0324 - 0325   >>> -------------------------------
% $$$ tid1 = find_trial_index(Trials,'jg05-20120323.cof.all');
% $$$ tid2 = find_trial_index(Trials,'jg05-20120324.cof.all');
% $$$ s = MTASession.validate(Trials{tid1}.filebase);
% $$$ s.spk.create(s);
% $$$ s.save();
% $$$ Trials{tid1} = MTATrial.validate(Trials{tid1}.filebase);
% $$$ pft = pfs_2d_theta(Trials{tid1},'overwrite',false);
% $$$ pft.purge_savefile();
% $$$ pft = pfs_2d_theta(Trials{tid1},'overwrite',true);
% $$$ s = MTASession.validate(Trials{tid2}.filebase);
% $$$ s.spk.create(s);
% $$$ s.save();
% $$$ Trials{tid2} = MTATrial.validate(Trials{tid2}.filebase);
% $$$ pfs = pfs_2d_theta(Trials{tid2},'overwrite',false);
% $$$ pfs.purge_savefile();
% $$$ pfs = pfs_2d_theta(Trials{tid2},'overwrite',true);
% $$$ 
% $$$ 
% $$$ shnk = 4
% $$$ cmtch{shnk} = [0,0];
% $$$ cmtch{shnk} = [ ...
% $$$     20, 24; ...
% $$$ ];
% $$$ 
% $$$ figure(1)
% $$$ clf();
% $$$ ucnt = sum(Trials{tid1}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid1}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,1)); continue; end;    
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pft,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ figure(2);
% $$$ clf();
% $$$ ucnt = sum(Trials{tid2}.spk.map(:,2)==shnk);
% $$$ offset = sum(Trials{tid2}.spk.map(:,2)<shnk);
% $$$ for k = 1:ucnt
% $$$     if ismember(k+offset, cmtch{shnk}(:,2)); continue; end;
% $$$     subplot(ceil(ucnt/10),10,k);
% $$$     plot(pfs,k+offset,[],'colorbar');
% $$$     title(num2str(k+offset));
% $$$ end
% $$$ 
% $$$ Trials{tid1}.spk.map(Trials{tid1}.spk.map(:,2)==shnk,:)
% $$$ Trials{tid2}.spk.map(Trials{tid2}.spk.map(:,2)==shnk,:)

% <<< I 0324 - 0325   <<< -------------------------------
% >>> I 0325 - 0326   >>> -------------------------------
% <<<   0325 - 0326   <<< -------------------------------
% <<< jg05 <<< ------------------------------------------

tmtch = cell([numel(Trials),numel(Trials)]);
% >>> Trial Name Pairs ----------------------------------
for tind1 = 1:numel(Trials)
    for tind2 = 1:numel(Trials)
        tmtch{tind1,tind2} = {Trials{tind1}.filebase,Trials{tind2}.filebase};
    end
end
% <<< Trial Name Pairs ----------------------------------

neuron_sets = generate_neuron_sets(tmtch, cmtch, false);



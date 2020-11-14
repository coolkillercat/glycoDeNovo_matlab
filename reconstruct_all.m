%% 20120327
clearvars filename gap minus2H;
datapath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20120327', filesep];
filename{1} = 'SLA.txt';

gap = 0;
minus2H = 0;
batch_reconstruct( datapath, filename, minus2H, gap );

%% 20130530
clearvars filename gap minus2H;
datapath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20130530', filesep];
filename{1} = 'Lewis B.EED.txt';
filename{2} = 'Lewis Y.EED.txt';
filename{3} = 'LNFP 1.EED.txt';
filename{4} = 'LNFP 2.EED.txt';
filename{5} = 'LNFP 3.EED.txt';

gap = [0 0 0 0 0];
minus2H = [0 0 1 0 0];
batch_reconstruct( datapath, filename, minus2H, gap );

%% 20130620
clearvars filename gap minus2H;
datapath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20130620', filesep];
filename{1} = 'Lewis B.CID.txt';
filename{2} = 'Lewis Y.CID.txt';

gap = [0 0];
minus2H = [0 0];
batch_reconstruct( datapath, filename, minus2H, gap );

%% 20140904
clearvars filename gap minus2H;
datapath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20140904', filesep];
filename{1} = 'A2F Reduced PM.txt';
filename{2} = 'Desialylated A2F O18 PM.txt';

gap = [0 0];
minus2H = [0 0];
batch_reconstruct( datapath, filename, minus2H, gap );

%% 20141029
clearvars filename gap minus2H;
datapath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20141029', filesep];
filename{1} = 'Man9_O18_PM.txt';

gap = 0;
minus2H = 0;
batch_reconstruct( datapath, filename, minus2H, gap );

%% 20151121
clearvars filename gap minus2H;
datapath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20151121', filesep];

id = 0;
id = id + 1; filename{id} = 'LNPF1 O18 PM Na.txt'; gap(id) = 1; minus2H(id) = 0;        % 1
id = id + 1; filename{id} = 'LNPF2 O18 PM Na.txt';  gap(id) = 0; minus2H(id) = 0;       % 2
id = id + 1; filename{id} = 'LNPF3 O18 PM Na.txt'; gap(id) = 0; minus2H(id) = 1;        % 3
id = id + 1; filename{id} = 'LewisB O18 PM Na.txt'; gap(id) = 0; minus2H(id) = 0;       % 4
id = id + 1; filename{id} = 'LewisY O18 PM Na.txt'; gap(id) = 0; minus2H(id) = 0;       % 5
id = id + 1; filename{id} = 'A2F Reduced PM Na.txt'; gap(id) = 0; minus2H(id) = 0;      % 6
id = id + 1; filename{id} = 'Desialyl-A2F O18 PM Na.txt'; gap(id) = 0; minus2H(id) = 0; % 7
id = id + 1; filename{id} = 'Man9 O18 PM Na.txt'; gap(id) = 0; minus2H(id) = 0;         % 8
id = id + 1; filename{id} = 'LNPF1 O18 PM Na CID.txt'; gap(id) = 1; minus2H(id) = 1;    % 9
id = id + 1; filename{id} = 'LNPF2 O18 PM Na CID.txt'; gap(id) = 1; minus2H(id) = 0;    % 10
id = id + 1; filename{id} = 'LNPF3 O18 PM Na CID.txt'; gap(id) = 1; minus2H(id) = 0;    % 11
id = id + 1; filename{id} = 'LewisB O18 PM Na CID.txt'; gap(id) = 0; minus2H(id) = 0;   % 12
id = id + 1; filename{id} = 'LewisY O18 PM Na CID.txt'; gap(id) = 0; minus2H(id) = 0;   % 13
id = id + 1; filename{id} = 'LNPF1 O18 PM Cs.txt'; gap(id) = 0; minus2H(id) = 0;        % 14
id = id + 1; filename{id} = 'LNPF2 O18 PM Cs.txt'; gap(id) = 0; minus2H(id) = 0;        % 15
id = id + 1; filename{id} = 'LNPF3 O18 PM Cs.txt'; gap(id) = 0; minus2H(id) = 0;        % 16
id = id + 1; filename{id} = 'LewisB O18 PM Cs.txt'; gap(id) = 0; minus2H(id) = 0;       % 17
id = id + 1; filename{id} = 'LewisY O18 PM Cs.txt'; gap(id) = 0; minus2H(id) = 0;       % 18

batch_reconstruct( datapath, filename, minus2H, gap );

%% 20160211
clear filename;
datapath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20160211', filesep];
filename{1} = 'N002 Reduced PM Na Deuterium.txt';
filename{2} = 'N003 Reduced PM Na Deuterium.txt';
filename{3} = 'N012 Reduced PM Na Deuterium.txt';
filename{4} = 'N013 Reduced PM Na Deuterium.txt';
filename{5} = 'N014 Reduced PM Na Deuterium.txt';
filename{6} = 'N222 Reduced PM Na Deuterium.txt';
filename{7} = 'N223 Reduced PM Na Deuterium.txt';
filename{8} = 'N233 Reduced PM Na Deuterium.txt';

gap = [0 0 0 0 1 1 0 0];
minus2H = [0 1 0 0 1 0 0 0];
batch_reconstruct( datapath, filename, minus2H, gap );

%% 20160410
clear filename;
datapath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20160410', filesep];
filename{1} = 'LewisB PM Na EED.20160410.txt';
filename{2} = 'LNnT PM Na EED.20160410.txt';
filename{3} = 'LNT PM Na EED.20160412.txt';

gap = [0 0 0];
minus2H = [0 0 0];
batch_reconstruct( datapath, filename, minus2H, gap );

%% 20170202
clear filename;
datapath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20170202', filesep];
filename{1} = 'LNFP-I.txt';
filename{2} = 'LNFP-II.txt';
filename{3} = 'LNFP-III.txt';
filename{4} = 'LNFP-V.txt';
filename{5} = 'LNFP-VI.txt';
filename{6} = 'Isomaltohexaose.txt';
filename{7} = 'Laminarihexaose.txt';
filename{8} = 'Maltohexaose.txt';

gap = [0 0 0 0 0 0 0 0];
minus2H = [0 0 0 0 0 0 0 0];
batch_reconstruct( datapath, filename, minus2H, gap );


%% 20170320
clearvars filename gap minus2H;
datapath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20170320', filesep];
filename{1} = 'Celhex PM Na.txt';
filename{2} = 'Lamhex O18 PM Na.txt';
filename{3} = 'LNnT PM.txt';
filename{4} = 'LNT PM.txt';
filename{5} = 'Malhex PM Na.txt';
filename{6} = 'SLX PM.txt';

gap = [0 0 0 0 0 0];
minus2H = [0 0 0 0 0 0];
batch_reconstruct( datapath, filename, minus2H, gap );

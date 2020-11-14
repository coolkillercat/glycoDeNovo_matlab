%% 2016/02/11
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

%%
specSet = CSpectrum.empty(length(filename), 0);
for f = 1 : length( filename )
    specSet{f} = CSpectrum.load( [datapath, filename{f}] );
    specSet{f}.protonate( specSet{f}.mPrecursor );
    specSet{f}.merge_peaks( 0.02 );
end

%% reconstruct topology
f = 6;
specU = specSet{f}.copy;
specU.add_complementary_ions();
specU.merge_peaks( 0.02 );

%%
disp( filename{f} );
reconstructor = CGlycanDeNovo( 0.02, [], minus2H(f), gap(f) );
reconstructor.mCheckMinusH = 0;
reconstructor.interpret_peaks( specU );
reconstructor.reconstruct_formulas();

%%
name = ['rec.' filename{f}];
specU.save_reconstruction( [datapath, 'results/', name], [], reconstructor.mCheckMinus2H*2, reconstructor.mCheckGap );
save( strrep( [datapath, 'results/', name], '.txt', '.mat' ), 'specU' );

%% Batch
batch_reconstruct( datapath, filename, minus2H, gap );
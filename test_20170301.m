%% 2017/03/01
clearvars filename gap minus2H;
datapath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20170301', filesep];
filename{1} = 'A2 unlabeled PM Cs.txt';
filename{2} = 'A2 unlabeled PM Na.txt';
filename{3} = 'LNnT unlabeled PM Cs.txt';
filename{4} = 'LNT unlabeled PM Cs.txt';
filename{5} = 'SLA unlabeled PM Cs.txt';
filename{6} = 'SLA unlabeled PM Na.txt';
filename{7} = 'SLX unlabeled PM Cs.txt';
filename{8} = 'SLX unlabeled PM Na.txt';

gap = [0 0 0 0 0 0 0 0];
minus2H = [1 0 0 0 0 0 0 0];

%%
specSet = CSpectrum.empty(length(filename), 0);
for f = 1 : length( filename )
    specSet{f} = CSpectrum.load( [datapath, filename{f}] );
    specSet{f}.protonate( specSet{f}.mPrecursor );
    specSet{f}.merge_peaks( 0.02 );
end

%% reconstruct topology
f = 4;
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
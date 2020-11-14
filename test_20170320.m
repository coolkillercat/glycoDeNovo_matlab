%% 2017/03/20
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

%%
specSet = CSpectrum.empty(length(filename), 0);
for f = 1 : length( filename )
    specSet{f} = CSpectrum.load( [datapath, filename{f}] );
    specSet{f}.protonate( [], specSet{f}.mPrecursor );
    specSet{f}.merge_peaks( 0.02 );
end

%% reconstruct topology
f = 1;
specU = specSet{f}.copy;
specU.add_complementary_ions();
specU.merge_peaks( 0.02 );

%%
disp( filename{f} );
reconstructor = CGlycanDeNovo( 0.02, [], minus2H(f), gap(f) );
reconstructor.mCheckMinusH = 0;
reconstructor.mIntensityThreshold = 0.002;
reconstructor.interpret_peaks( specU );

reconstructor.reconstruct_formulas();

%%
name = ['rec.' filename{f}];
specU.save_reconstruction( [datapath, 'results/', name], [], reconstructor.mCheckMinus2H*2, reconstructor.mCheckGap );
save( strrep( [datapath, 'results/', name], '.txt', '.mat' ), 'specU' );

%% Batch
batch_reconstruct( datapath, filename, minus2H, gap );
%% 207/02
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

%%
specSet = CSpectrum.empty(length(filename), 0);
for f = 1 : length( filename )
    specSet{f} = CSpectrum.load( [datapath, filename{f}] );
    specSet{f}.protonate( 1 );
    specSet{f}.merge_peaks( 0.001 );
end

%% reconstruct topology
f = 5;
specU = specSet{f}.copy;
specU.add_complementary_ions();
specU.merge_peaks( 0.001 );

%%
reconstructor = CGlycoDeNovo( 5, [], minus2H(f), gap(f) );
reconstructor.mCheckMinusH = 0;
reconstructor.mUseComplementOnly = 1;
reconstructor.interpret_peaks( specU );

%%
reconstructor.reconstruct_formulas();

%%
name = ['rec.' filename{f}];
specU.save_reconstruction( [datapath, 'results/', name], [], reconstructor.mCheckMinus2H*2, reconstructor.mCheckGap );
save( strrep( [datapath, 'results/', name], '.txt', '.mat' ), 'specU' );

%% Batch
batch_reconstruct( datapath, filename, minus2H, gap, 5, 0, 1, 1 );
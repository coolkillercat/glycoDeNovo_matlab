%% 2018/04/26
clearvars filename gap minus2H;
datapath = 'D:\Projects\Glycomics\data\20180426\';
filename{1} = 'spectrum1.txt';
filename{2} = 'spectrum2.txt';
filename{3} = 'spectrum3.txt';
filename{4} = 'spectrum4.txt';
filename{5} = 'spectrum5.txt';

gap = [0 0 0 0 0];
minus2H = [0 0 0 0 0];

%%
specSet = CSpectrum.empty(length(filename), 0);
for f = 1 : length( filename )
    specSet{f} = CSpectrum.load( [datapath, filename{f}] );
    specSet{f}.protonate( specSet{f}.mPrecursor );
    %specSet{f}.merge_peaks( 0.02 );
end

%% reconstruct topology
f = 5;
specU = specSet{f}.copy;
specU.add_complementary_ions();
specU.merge_peaks( 0.001 );

%%
disp( filename{f} );
reconstructor = CGlycoDeNovo( 5, [], minus2H(f), gap(f) ); % 5ppm accuracy
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



%% IonClassifier -- load training data
datapath = 'D:\Projects\Glycomics\data\20180426\';
spectraTrain = load_saved_spectra( '20170615' );
spectraTest = load_saved_spectra( '20180426' );

%% IonClassifier -- train classifier
ionClassifier = CIonClassifier;
ionClassifier.mMassAccuracy = 0.01;
ionClassifier.train( spectraTrain );
ionClassifier.rank_candidates( spectraTest );

%% IonClassifier -- save results
for k = 1 : length( spectraTest )
    name = strsplit( spectraTest(k).filename, '\' );
    name = strrep( name{end}, '.mat', '.txt' );
    name = strrep( name, 'rec.', 'IC.' );
    disp( name );
    if ~isempty( spectraTest(k).mPeaks(end).mInferredSuperSet )
        spectraTest(k).mPeaks(end).mInferredSuperSet.sort_topologies_by_score();
    end
    spectraTest(k).save_reconstruction( [datapath, 'results/', name], [] );
end


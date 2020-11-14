%% 2018/09/04
clear filename;
datapath = ['D:\Projects\Glycomics\data\', '20180822', filesep];
id = 1; gap = []; minus2H = [];
filename{id} = 'LNFP I.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'LNFP III.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 2
filename{id} = 'LNFP VI.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 3

%%
specSet = CSpectrum.empty(length(filename), 0);
for f = 1 : length( filename )
    specSet{f} = CSpectrum.load( [datapath, filename{f}] );
    specSet{f}.protonate( specSet{f}.mPrecursor );
    specSet{f}.merge_peaks( 0.001 );
end

%% reconstruct topology
for f = 1 : length( filename )
    specU = specSet{f}.copy;
    specU.add_complementary_ions();
    specU.merge_peaks( 0.001 );
    
    %
    reconstructor = CGlycoDeNovo( 5, [], minus2H(f), gap(f) ); % 5ppm accuracy
    reconstructor.mCheckMinusH = 0;
    reconstructor.interpret_peaks( specU );
    
    %
    reconstructor.reconstruct_formulas();
    
    %
    name = ['rec.' filename{f}];
    specU.save_reconstruction( [datapath, 'results', filesep, name], [], reconstructor.mCheckMinus2H*2, reconstructor.mCheckGap );
    save( strrep( [datapath, 'results', filesep, name], '.txt', '.mat' ), 'specU' );
    disp( filename{f} );
end

%% Batch
% batch_reconstruct( datapath, filename, minus2H, gap );
batch_reconstruct( datapath, filename, minus2H, gap, 5, 0, 1 );


%% IonClassifier -- load training data
datapath = 'D:\Projects\Glycomics\data\20180822\';
spectraTrain = load_saved_spectra( '20170615.2' ); % or use '20170615'
spectraTest = load_saved_spectra( '20180822' );

%% IonClassifier -- train classifier
ionClassifier = CIonClassifier;
ionClassifier.mMassAccuracy = 0.01;
ionClassifier.train( spectraTrain );
ionClassifier.rank_candidates( spectraTest );
disp( 'rank candidates done' );

%% IonClassifier -- save results
disp( 'Saving results ...' );
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
disp( 'done' );

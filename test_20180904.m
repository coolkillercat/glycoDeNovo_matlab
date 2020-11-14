%% 2018/09/04
clear filename;
datapath = ['D:', filesep, 'Projects', filesep, 'Glycomics', filesep, 'data', filesep, '20180904', filesep];
id = 1; gap = []; minus2H = [];
filename{id} = 'TIMS_1046.1.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'TIMS_1046.2.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'TIMS_1046.3.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 1

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
    %specU.merge_peaks( 0.001 );
    
    %
    reconstructor = CGlycoDeNovo( 5, [], minus2H(f), gap(f) ); % 5ppm accuracy
    reconstructor.mCheckMinusH = 0;
    reconstructor.mMaxBranchingNum = 3;
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
disp( 'Applying IonClassifier ...' );
datapath = 'D:\Projects\Glycomics\data\20180904\';
spectraTrain = load_saved_spectra( '20170615.2' ); % or use '20170615'
spectraTest = load_saved_spectra( '20180904' );

%%
ionClassifier = CIonClassifier;
ionClassifier.mMassAccuracy = 0.01;
ionClassifier.mUseOriginalPeaks = 0;
[trainVectors, trainData] = ionClassifier.train( spectraTrain );

ionClassifier.rank_candidates( spectraTest );

disp( 'Saving results ...' );
for k = 1 : length( spectraTest )
    name = strsplit( spectraTest(k).filename, '\' );
    name = strrep( name{end}, '.mat', '.txt' );
    if ionClassifier.mUseOriginalPeaks
        name = strrep( name, 'rec.', 'IC.OriginalPeak.' );
    else
        name = strrep( name, 'rec.', 'IC.' );
    end
    disp( name );
    if ~isempty( spectraTest(k).mPeaks(end).mInferredSuperSet )
        spectraTest(k).mPeaks(end).mInferredSuperSet.sort_topologies_by_score();
    end
    spectraTest(k).save_reconstruction( [datapath, 'results/', name], [] );
end
disp( 'done' );

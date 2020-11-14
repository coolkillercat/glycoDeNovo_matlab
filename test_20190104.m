%% 2018/09/04
clear filename;
datapath = ['D:\Projects\Glycomics\data\', '20190104', filesep];
id = 1; gap = []; minus2H = [];
filename{id} = 'Man4_Peak1_OLD.txt'; gap(id) = 1; minus2H(id) = 1; id = id + 1; % 1
filename{id} = 'Man4_Peak2_OLD.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 2
filename{id} = 'Man4_Peak3_OLD.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 3
filename{id} = 'Man4_Peak1_NEW.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 4
filename{id} = 'Man5_Peak2_OLD.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 5
filename{id} = 'Man6_Peak3_OLD.txt'; gap(id) = 0; minus2H(id) = 1; id = id + 1; % 6
filename{id} = 'Man4_Peak1_20190217.txt'; gap(id) = 0; minus2H(id) = 1; id = id + 1; % 7

%%
specSet = CSpectrum.empty(length(filename), 0);
for f = 1 : length( filename )
    specSet{f} = CSpectrum.load( [datapath, filename{f}] );
    specSet{f}.protonate( specSet{f}.mPrecursor );
    specSet{f}.merge_peaks( 0.001 );
end
disp('Done loading raw data ...');

%% reconstruct topology
for f = 7 : 7 %1 : length( filename )
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
datapath = 'D:\Projects\Glycomics\data\20190104\';
spectraTrain = load_saved_spectra( '20170615.2' ); % or use '20170615'
spectraTest = load_saved_spectra( '20190104' );

%% IonClassifier -- train classifier
ionClassifier = CIonClassifier;
ionClassifier.mMassAccuracy = 0.01;
ionClassifier.train( spectraTrain );
disp( 'Train IonClassifier done.' );

%%
ionClassifier.mEnableX15 = 1;
ionClassifier.rank_candidates( spectraTest );
disp( 'Test IonClassifier done.' );

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

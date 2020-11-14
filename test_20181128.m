%% 2018/11/28
clear filename;
if ismac
    datapath = [getenv('HOME'), '/Documents/Projects/Glycomics/data/', '20181128', filesep];
else
    datapath = ['D:\Projects\Glycomics\data\', '20190116', filesep];
end
id = 1; gap = []; minus2H = []; branchNum = [];
filename{id} = 'LNFP V.txt'; gap(id) = 1; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 1
filename{id} = 'SLeA.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 1
filename{id} = 'SLeX.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 3; id = id + 1; % 1

%%
specSet = CSpectrum.empty(length(filename), 0);
for f = 1 : length( filename )
    specSet{f} = CSpectrum.load( [datapath, filename{f}] );
    specSet{f}.protonate( 1 );
    % specSet{f}.merge_peaks( 0.001 );
end

%% reconstruct topology
for f = 1 : length( filename )
    specU = specSet{f}.replicate;
    specU.add_complementary_ions();
    %specU.merge_peaks( 0.001 );
    
    %
    reconstructor = CGlycoDeNovo( 5, [], minus2H(f), gap(f) ); % 5ppm accuracy
    reconstructor.mCheckMinusH = 0;
    reconstructor.mMaxBranchingNum = branchNum(f);
    reconstructor.interpret_peaks( specU );
    
    %
    reconstructor.reconstruct_formulas();
    
    %
    name = ['rec.' filename{f}];
    specU.save_reconstruction( [datapath, 'results', filesep, name], [], reconstructor.mCheckMinus2H*2, reconstructor.mCheckGap );
    save( strrep( [datapath, 'results', filesep, name], '.txt', '.mat' ), 'specU' );
    % specU.save_reconstruction_full( [datapath, 'results', filesep, strrep(name, '.txt', '.full.txt')], reconstructor.mCheckMinus2H*2, reconstructor.mCheckGap );
    % save( strrep( [datapath, 'results', filesep, name], '.txt', '.mat' ), 'specU' );
    disp( filename{f} );
end


%% Batch
% batch_reconstruct( datapath, filename, minus2H, gap );
batch_reconstruct( datapath, filename, minus2H, gap, 5, 0, 1 );

%% IonClassifier -- load training data
disp( 'Applying IonClassifier ...' );
datapath = ['D:', filesep, 'Projects', filesep, 'Glycomics', filesep, 'data', filesep, '20181128', filesep];
spectraTrain = load_saved_spectra( '20170615.2' ); % or use '20170615'
spectraTest = load_saved_spectra( '20181128' );

%%
ionClassifier = CIonClassifier;
ionClassifier.mMassAccuracy = 0.01;
ionClassifier.mUseOriginalPeaks = 0;
[trainVectors, trainData] = ionClassifier.train( spectraTrain );
disp('Done training.');

%%
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
    % spectraTest(k).save_reconstruction( [datapath, 'results/', name], [] );
    disp([datapath, 'results', filesep, strrep(name, '.txt', '.full.txt')]);
    spectraTest(k).save_reconstruction_full( [datapath, 'results', filesep, strrep(name, '.txt', '.full.txt')], minus2H(k)*2, gap(k) );
end
disp( 'done' );



%% For Cheng Lin App Note %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spectraTrain = load_saved_spectra( '20170615.2' ); % or use '20170615'
backup = spectraTrain;

%%
flag = ones(1, length(backup));
flag(end-1) = 0;
spectraTrain = backup(flag>0);
spectraTest = backup(flag == 0);

%%
ionClassifier = CIonClassifier;
ionClassifier.mMassAccuracy = 0.01;
ionClassifier.mUseOriginalPeaks = 0;
[trainVectors, trainData] = ionClassifier.train( spectraTrain );
disp('Done training.');

%%
resultpath = ['D:', filesep, 'Projects', filesep, 'Glycomics', filesep, 'data', filesep, '20181128', filesep];
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
    name = strrep( name, '.txt', '.NoReducing.txt' );
    disp( name );
    if ~isempty( spectraTest(k).mPeaks(end).mInferredSuperSet )
        spectraTest(k).mPeaks(end).mInferredSuperSet.sort_topologies_by_score();
    end
    spectraTest(k).save_reconstruction( [resultpath, 'results/', name], [] );
end
disp( 'done' );


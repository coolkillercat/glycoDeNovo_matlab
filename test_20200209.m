clear filename;
if ismac
    datapath = [getenv('HOME'), '/Documents/Projects/Glycomics/data/', '20200209', filesep];
else
    datapath = ['D:\Projects\Glycomics\data\', '20200209', filesep];
end
files = dir( [datapath, '*.txt'] );

%%
specSet = CSpectrum.empty(length(files), 0);
for f = 1 : length( files )
    filename = [datapath, files(f).name];
    disp( ['Loading ', filename] );
    specSet{f} = CSpectrum.load( filename );
end
disp('Done loading raw data ...');

%% activate mass2composition
m2c = CMass2Composition;
m2c.load( 'm2c.mat' );
m2c.mCheckMinus2H = 1;

%% reconstruct topology
debugMatchedSpectrum = 1;
for f = 1 : length( files )
    disp( '>> ' );
    disp( ['>> ', specSet{f}.filename] );
    aSpec = specSet{f}.replicate;
    aSpec.protonate();
    
    if aSpec.mMassAccuracy > 0
        m2c.mMassAccuracyPPM = aSpec.mMassAccuracy;
    else
        m2c.mMassAccuracyPPM = 5;
    end
    m2c.set_reducing_end_modification( aSpec.mReducingEndModification );
    m2c.set_permethylation( aSpec.mPermethylated );
    hypoSpecs = m2c.correct_spectrum( aSpec );
    
    trialNum = length(hypoSpecs);
    for s = 1 : trialNum
        disp( ['>>> Try ', num2str(s)] );
        reconstructor = CGlycoDeNovo( 5, find(hypoSpecs(s).mComposition > 0), 0, 0 );
        reconstructor.mCheckMinusH = 0;
        reconstructor.mMaxBranchingNum = 3;
        reconstructor.mPermethylated = aSpec.mPermethylated;
        reconstructor.mCompositionConstraint = hypoSpecs(s).mComposition;
        reconstructor.set_reducing_end_modification( hypoSpecs(s).mReducingEndModification );
        reconstructor.interpret_peaks( hypoSpecs(s) );
        reconstructor.reconstruct_formulas();
        if trialNum == 1
            name = strrep(files(f).name, '.txt', '.grec');
        else
            name = strrep(files(f).name, '.txt', ['.', num2str(s), '.grec']);
        end
        
        resultfile = [datapath, 'results', filesep, name];
        if isempty(hypoSpecs(s).mPeaks(end).mInferredFormulas)
            disp( '>>> Failed, try -2H' );
            reconstructor.mCheckMinus2H = 1;
            reconstructor.interpret_peaks( hypoSpecs(s) );
            reconstructor.reconstruct_formulas();
            
            if isempty(hypoSpecs(s).mPeaks(end).mInferredFormulas)
                debug_reconstruction(hypoSpecs(s), [datapath, 'results', filesep, strrep(name, '.grec', '.gdbg')], ...
                                    debugMatchedSpectrum );
            else
                disp( ['>>> Save reconstruction results in ', resultfile] );
                hypoSpecs(s).save_reconstruction( resultfile, [], ...
                                        reconstructor.mCheckMinus2H*2, reconstructor.mCheckGap );
            end
            reconstructor.mCheckMinus2H = 0;
        else
            disp( ['>>> Save reconstruction results in ', resultfile] );
            hypoSpecs(s).save_reconstruction( resultfile, [], ...
                                        reconstructor.mCheckMinus2H*2, reconstructor.mCheckGap );
        end
    end
    
    disp( ['>> Finished ', specSet{f}.filename] );
end
disp( 'Finished all.' );

%%
specU = specSet{1}.replicate;
specU.protonate();
specU.add_complementary_ions();
specU.merge_peaks( 0.001 );
reconstructor = CGlycoDeNovo( 5, [], 0, 0 );
reconstructor.mCheckMinusH = 0;
reconstructor.mMaxBranchingNum = 3;
reconstructor.interpret_peaks( specU );
reconstructor.reconstruct_formulas();
specU.save_reconstruction( 'test.txt', [], reconstructor.mCheckMinus2H*2, reconstructor.mCheckGap );


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
for k = 1 : 1 %length( spectraTest )
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


%%
id = 1; gap = []; minus2H = [];
filename{id} = 'G0.1b2.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 1
filename{id} = 'G0B.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 3; id = id + 1; % 2
filename{id} = 'G0F.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 3
filename{id} = 'G0FB.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 3; id = id + 1; % 4
filename{id} = 'G1.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 5
filename{id} = 'G1B.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 3; id = id + 1; % 6
filename{id} = 'G1F.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 7
filename{id} = 'G1FB.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 8
filename{id} = 'G2.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 9
filename{id} = 'G2F.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 10
filename{id} = 'G2FB.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 3; id = id + 1; % 11
filename{id} = 'G2S1.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 12
filename{id} = 'G2S1B.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 13
filename{id} = 'G2S1F.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 14
filename{id} = 'G2S1FB.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 15
filename{id} = 'G2S2.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 16
filename{id} = 'G2S2F.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 17
filename{id} = 'G2S2FB.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 18
filename{id} = 'G3S1F.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 19
filename{id} = 'G3S2B.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 3; id = id + 1; % 20
filename{id} = 'G3S3.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 21
filename{id} = 'G3S3F.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 22
filename{id} = 'M3G1S.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 23
filename{id} = 'M4G0.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 24
filename{id} = 'M4G0F.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 25
filename{id} = 'M4G1.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 26
filename{id} = 'M5G1.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 27
filename{id} = 'Man5.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 28
filename{id} = 'Man6.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 29
filename{id} = 'Man7.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 30
filename{id} = 'Man8.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 31
filename{id} = 'Man9.txt'; gap(id) = 0; minus2H(id) = 1; branchNum(id) = 2; id = id + 1; % 32

% This dataset use SNAP peak picking

clear filename;
dataset = '20170615.2';
if ismac
    datapath = [getenv('HOME'), '/Documents/Projects/Glycomics/data/', dataset, filesep];
else
    datapath = ['E:\Brandeis\data\', dataset, filesep];
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
    
    % Begin - 1: This is for copy_reconstruction and training ionclassifier
    aSpec.add_complementary_ions();
    aSpec.merge_peaks( 0.0025 );
    % End - 1

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
        reconstructor.mMaxBranchingNum = 2;
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
        
        copy_reconstruction( hypoSpecs(s), aSpec );
    end
    
    % save it for IonClassifier
    save( [datapath, 'results', filesep, strrep(files(f).name, '.txt', '.mat')], 'aSpec' );
    
    disp( ['>> Finished ', specSet{f}.filename] );
end
disp( 'Finished all.' );

%%
trainSpec = load_saved_spectra( '20190116_FTMS' );
allSpec = load_saved_spectra( '20200618' );

%% IonClassifier -- train classifier
ionClassifier = CIonClassifier2;
ionClassifier.mMassAccuracy = 0.01;

%%
disp( 'IonClassifier working ...' );
ionClassifier.train( trainSpec );

%%
for k = 1 : length( allSpec )
    aSpec = allSpec(k);
    if ~isempty( aSpec.mPeaks(end).mInferredFormulas )
        
        filename = strrep( aSpec.filename, '.mat', '.IC');
        disp( filename );
        
        ionClassifier.score_candidates( aSpec );
        aSpec.sort_topologies_by_score();
        aSpec.save_reconstruction( filename, [] );
    end
end
disp( 'done' );

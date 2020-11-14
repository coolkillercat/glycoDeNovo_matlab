%% 2018/09/04
clear filename;
clear filename;
if ismac
    datapath = [getenv('HOME'), '/Documents/Projects/Glycomics/data/', '20180816', filesep];
else
    datapath = ['D:\Projects\Glycomics\data\', '20180816', filesep];
end
id = 1; gap = []; minus2H = []; branchNum = [];
filename{id} = 'A2S2.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'HighMan_PGC.Man5_Peak1.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 2
filename{id} = 'HighMan_PGC.Man5_Peak2.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 3
filename{id} = 'HighMan_PGC.Man5_Peak3.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 4
filename{id} = 'HighMan_PGC.Man6_Peak1.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 5
filename{id} = 'HighMan_PGC.Man6_Peak2.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 6
filename{id} = 'HighMan_PGC.Man6_Peak3.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 7
filename{id} = 'HighMan_PGC.Man7_Peak1.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 8
filename{id} = 'HighMan_PGC.Man7_Peak2.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 9
filename{id} = 'HighMan_PGC.Man7_Peak3.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 10
filename{id} = 'HighMan_PGC.Man7_Peak4.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 11
filename{id} = 'HighMan_PGC.Man7_Peak5.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 12
filename{id} = 'HighMan_PGC.Man8_Peak1.txt'; gap(id) = 1; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 13
filename{id} = 'HighMan_PGC.Man8_Peak2.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 14
filename{id} = 'HighMan_PGC.Man8_Peak3.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 15
filename{id} = 'HighMan_PGC.Man9.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 16
filename{id} = 'Man5.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 17
filename{id} = 'Man6.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 18
filename{id} = 'NA2F.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 19
filename{id} = 'NG1A2F.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 20
filename{id} = 'NGA2F.txt'; gap(id) = 0; minus2H(id) = 0; branchNum(id) = 2; id = id + 1; % 21

%%
specSet = CSpectrum.empty(length(filename), 0);
for f = 1 : length( filename )
    specSet{f} = CSpectrum.load( [datapath, filename{f}] );
    specSet{f}.protonate( 1 );
    %specSet{f}.merge_peaks( 0.001 );
end

%% reconstruct topology
for f = 13 : 13 %1 : length( filename )
    specU = specSet{f}.replicate();
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
datapath = 'D:\Projects\Glycomics\data\20180816\';
spectraTrain = load_saved_spectra( '20170615.2' ); % or use '20170615'
spectraTest = load_saved_spectra( '20180816' );

%% IonClassifier -- train classifier
ionClassifier = CIonClassifier;
ionClassifier.mMassAccuracy = 0.01;
ionClassifier.train( spectraTrain );
ionClassifier.rank_candidates( spectraTest );

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

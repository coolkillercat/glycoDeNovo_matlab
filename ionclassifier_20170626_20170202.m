%%
% Use data with non-PRAGS reducing end moficiations to train IonClassifier, and apply
% the trained IonClassifier to the data with PRAGS reducing end moficiation.
%%
spectra = load_saved_spectra( {'20170615', '20170202'} );
IonClassifierResultPath = 'D:\Projects\Glycomics\data\20170202\results\';

%%
flagTrain = zeros(1, length(spectra));
flagTest = zeros(1, length(spectra));
for k = 1 : length(flagTrain)
    if strcmp( spectra(k).mExperimentMethod, 'EED' )
        if strcmp( spectra(k).mReducingEndModification, 'PRAGS' )
            flagTest(k) = 1;
        elseif ~isempty( spectra(k).mReducingEndModification )
            flagTrain(k) = 1;
        end
    end
end
spectraTrain = spectra(flagTrain > 0);
spectraTest = spectra(flagTest > 0);

%%
[vectors, data] = prepare_training_data( spectra, 0.01, 0, 0 );
trainData = data; trainVectors = vectors;
testData = data; testVectors = vectors;
ions = {'B', 'C', 'Y', 'Z', 'O'};
for k = 1 : length(ions)
    ion = ions{k};
    disp(ion);
    flagA = zeros(1, length(data.(ion)));
    flagB = zeros(1, length(data.(ion)));
    for m = 1 : length(data.(ion))
        flagA(m) = flagTrain( data.(ion)(m).glycanID );
        flagB(m) = flagTest( data.(ion)(m).glycanID );
    end
    trainData.(ion) = data.(ion)(flagA > 0);
    testData.(ion) = data.(ion)(flagB > 0);
    trainVectors.(ion) = vectors.(ion)(flagA > 0, :);
    testVectors.(ion) = vectors.(ion)(flagB > 0, :);
end

%%
ionClassifier = train_ionclassifier( trainVectors, trainData, 0 );
results = predict_by_ionclassifier( spectraTest, testData, testVectors, ionClassifier );
xlswrite( [IonClassifierResultPath, 'IonClassifier.20170906.xlsx'], results );

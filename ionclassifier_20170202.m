spectra = load_saved_spectra( {'20170202'} );

%% LOO
looData = spectra;
[vectors, data] = prepare_training_data( looData, 0.01, 0, 0 );
trainVectors.massFeatures = vectors.massFeatures;
testVectors.massFeatures = vectors.massFeatures;

%%
ions = ['B', 'C', 'Y', 'Z', 'O'];
CVResults = [];
clear trainData testData trainVectors testVectors;
trainVectors.massFeatures = vectors.massFeatures;
testVectors.massFeatures = vectors.massFeatures;
% for k = 1 : length(looData)
for k = 6 : 6
	flag = ones(1, length(looData));
	flag(k) = 0;

	for ion = ions
		flagSample = ones(1, length(data.(ion)));
		for m = 1 : length(data.(ion))
			flagSample(m) = flag( data.(ion)(m).glycanID );
		end
		trainData.(ion) = data.(ion)(flagSample > 0);
		testData.(ion) = data.(ion)(flagSample == 0);
		trainVectors.(ion) = vectors.(ion)(flagSample > 0, :);
		testVectors.(ion) = vectors.(ion)(flagSample == 0, :);
	end
	
	ionClassifier = train_ionclassifier( trainVectors, trainData, 0 );
	[aResult, icScores] = predict_by_ionclassifier( looData(k), testData, testVectors, ionClassifier );
    aResult{end,1} = k;
	if isempty( CVResults )
		CVResults = aResult;
	else % results(end+1,:)
		CVResults(end+1,:) = aResult(end,:);
	end
end
% IonClassifierResultPath = [get_project_path, 'Glycomics', filesep, 'denovo', filesep, 'data', filesep, '20170202', filesep, 'results'];

% xlswrite( [IonClassifierResultPath, 'IonClassifier.LOO_RD.20170202.xlsx'], CVResults );
% disp( 'Finished LOO' );

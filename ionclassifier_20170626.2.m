%%
% Use data with reducing end moficiations to train IonClassifier, and apply
% the trained IonClassifier to the data without reducing end moficiation.
%%
spectra = load_saved_spectra( '20170615.2' );
IonClassifierResultPath = 'D:\hong\Documents\Projects\Glycomics\denovo\data\20170615.2\results\';

%%
flagTrain = zeros(1, length(spectra));
for k = 1 : length(flagTrain)
    flagTrain(k) = strcmp( spectra(k).mExperimentMethod, 'EED' ) & ~isempty( spectra(k).mReducingEndModification );
end
spectraTrain = spectra(flagTrain > 0);
spectraTest = spectra(flagTrain == 0);

%%
[vectors, data] = prepare_training_data( spectra, 0.04, 0, 0 );
trainData = data; trainVectors = vectors;
testData = data; testVectors = vectors;
ions = {'B', 'C', 'Y', 'Z', 'O'};
for k = 1 : length(ions)
    ion = ions{k};
    disp(ion);
    flagA = ones(1, length(data.(ion)));
    for m = 1 : length(data.(ion))
        flagA(m) = flagTrain( data.(ion)(m).glycanID );
    end
    trainData.(ion) = data.(ion)(flagA > 0);
    testData.(ion) = data.(ion)(flagA == 0);
    trainVectors.(ion) = vectors.(ion)(flagA > 0, :);
    testVectors.(ion) = vectors.(ion)(flagA == 0, :);
end

%%

ionClassifier = train_ionclassifier( trainVectors, trainData, 0 );
results = predict_by_ionclassifier( spectraTest, testData, testVectors, ionClassifier );
xlswrite( [IonClassifierResultPath, 'IonClassifier.trained_RD.20170906.xlsx'], results );

%% LOO
looData = spectraTrain;
[vectors, data] = prepare_training_data( looData, 0.04, 0, 0 );
trainVectors.massFeatures = vectors.massFeatures;
testVectors.massFeatures = vectors.massFeatures;

ions = ['B', 'C', 'Y', 'Z', 'O'];
CVResults = [];
for k = 1 : length(looData)
% for k = 5 : 5
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
	aResult = predict_by_ionclassifier( looData(k), testData, testVectors, ionClassifier );
	if isempty( CVResults )
		CVResults = aResult;
	else % results(end+1,:)
		CVResults(end+1,:) = aResult(end,:);
    end
end
xlswrite( [IonClassifierResultPath, 'IonClassifier.LOO_RD.20170906.xlsx'], CVResults );

%% check features
clear selectedFeatures;
num = length( trainVectors.massFeatures );
fieldNames = fields( ionClassifier );
for k = 1 : length( fieldNames )
    aField = fieldNames{k};
    selectedFeatures.(aField).idx = ionClassifier.(aField).SelectedFeatureIdx;
    selectedFeatures.(aField).importance = ionClassifier.(aField).FeatureImportance(selectedFeatures.(aField).idx);
    idx = selectedFeatures.(aField).idx; 
    idx(idx > num) = idx(idx>num)-num;
    selectedFeatures.(aField).massFeatures = trainVectors.massFeatures(idx);
    [selectedFeatures.(aField).importance, idx] = sort( selectedFeatures.(aField).importance, 'descend' );
    selectedFeatures.(aField).idx = selectedFeatures.(aField).idx(idx);
    selectedFeatures.(aField).massFeatures = selectedFeatures.(aField).massFeatures(idx);
end

%% check A2F reduced
gtRank = spectra(2).get_rank( '[[Neu5Ac Hex HexNAc Hex] [Neu5Ac Hex HexNAc Hex] Hex HexNAc] [Fuc] HexNAc' );

%%
peak = spectra(2).mPeaks(end);
peakset = [];
for k = 1 : gtRank.idx-1
    peakset = [peakset, peak.mInferredSuperSet.mTopologies(k).mSupportPeaks];
    if length(peakset) > 1000
        peakset = unique(peakset);
    end
end
peakset = unique( peakset );
fakePeaks = setdiff( peakset, peak.mInferredSuperSet.mTopologies( gtRank.idx ).mSupportPeaks );
fakePeakMasses = [spectra(2).mPeaks.mMass];
fakePeakMasses = fakePeakMasses( fakePeaks);

%% Check A2F reduced (2)
% Peak 167: mass 406.207116
% Peak 252: mass 464.249016
c = 0;
for k = 1 : gtRank.idx-1
    if any( peak.mInferredSuperSet.mTopologies(k).mSupportPeaks == 167 )
        c = c + 1;
    end
end

%% Peaks, Interpretable Peak, Recontructed Peaks, Candidate
temp = [329	133	18	6	2
216	76	24	8	4
461	193	28	8	4
283	105	26	6	2
469	209	45	19	16
516	224	23	11	13
390	178	26	14	16
534	245	32	12	1
471	212	24	11	10
477	210	21	13	17
2389	1109	395	24	22
2532	1182	588	101	1870
2646	1222	597	151	990750
914	435	71	25	37
2320	1063	262	52	116290
1571	731	175	49	834
2683	1229	273	50	4619
2544	1179	351	48	2385
953	411	78	18	34
2674	1189	226	30	1577
2326	1078	234	33	1920
218	91	30	9	4	
317	126	21	7	5
270	105	23	9	5
459	195	48	17	14
333	125	55	18	22
412	166	47	11	11
468	207	58	18	22
];

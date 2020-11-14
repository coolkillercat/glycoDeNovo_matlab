spectraTrain = load_saved_spectra( '20170615.2' ); % or use '20170615'
ionClassifier = CIonClassifier;
ionClassifier.mMassAccuracy = 0.01;
[trainVectors, trainData] = ionClassifier.train( spectraTrain );

%%
fprintf('====\n');
dX15 = -27.9949;
idxes = find( abs(trainVectors.massFeatures - dX15) < 0.005 );
fprintf('Feature %f\n', dX15);
for ion = 'BCYZO'
    fprintf( '%s: %f\n', ion, sum(sum(trainVectors.(ion)(:,idxes))) / size(trainVectors.(ion),1) );
end
fprintf('Feature %f\n', -dX15);
idxes = find( abs(trainVectors.massFeatures + dX15) < 0.005 );
for ion = 'BCYZO'
    fprintf( '%s: %f\n', ion, sum(sum(trainVectors.(ion)(:,idxes))) / size(trainVectors.(ion),1) );
end
fprintf('----\n');
dX15 = -27.9949 - CMass.H2O;
idxes = find( abs(trainVectors.massFeatures - dX15) < 0.005 );
fprintf('Feature %f\n', dX15);
for ion = 'BCYZO'
    fprintf( '%s: %f\n', ion, sum(sum(trainVectors.(ion)(:,idxes))) / size(trainVectors.(ion),1) );
end
fprintf('Feature %f\n', -dX15);
idxes = find( abs(trainVectors.massFeatures + dX15) < 0.005 );
for ion = 'BCYZO'
    fprintf( '%s: %f\n', ion, sum(sum(trainVectors.(ion)(:,idxes))) / size(trainVectors.(ion),1) );
end
fprintf('====\n');

%%
spectraTest = load_saved_spectra( '20180816' );

%%
[testVectors, signals] = CIonClassifier.extract_data( spectraTest(13), ionClassifier.mMassFeatures, [], 0.01, 0 );
disp(sum(sum(testVectors.B(:,idxes))) / size(testVectors.B,1))

%%
ionClassifier.rank_candidates( spectraTest(13) );

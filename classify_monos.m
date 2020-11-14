spectra = load_saved_spectra( '20170615' );
[allVectors, allData] = CIonClassifier.prepare_training_data( spectra, [], [], 0.005, 0, 0 );
[glcData, glcVectors] = CIonClassifier.filter_BC_by_rootmono( allData, allVectors, {'Glc'});
[galData, galVectors] = CIonClassifier.filter_BC_by_rootmono( allData, allVectors, {'Gal'});
[manData, manVectors] = CIonClassifier.filter_BC_by_rootmono( allData, allVectors, {'Man'});

% save monoSpecificData manData manVectors glcData glcVectors galData galVectors;

%% 
% load monoSpecificData

%% 
ion = 'C';
% c1 = 'man'; c2 = 'gal';
c1 = 'gal'; c2 = 'glc';
if strcmp( c1, 'man' )
    if strcmp(c2, 'gal')
        X = [manVectors.(ion); galVectors.(ion)];
        Y = [ones(size(manVectors.(ion), 1), 1); zeros(size(galVectors.(ion), 1), 1)];
    elseif strcmp(c2, 'glc')
        X = [manVectors.(ion); glcVectors.(ion)];
        Y = [ones(size(manVectors.(ion), 1), 1); zeros(size(glcVectors.(ion), 1), 1)];
    end
elseif strcmp( c1, 'gal' )
    X = [galVectors.(ion); glcVectors.(ion)];
    Y = [ones(size(galVectors.(ion), 1), 1); zeros(size(glcVectors.(ion), 1), 1)];
end

%%
disp( ['LOO ' c1 ' vs ' c2] );
N = size(Y,1);
looresult = zeros(1, N);
for k = 1 : N
    flag = ones(N,1);
    flag(k) = 0; flag = flag > 0;
    trainX = X(flag,:);
    trainY = Y(flag,:);
    cc = fitensemble( trainX, trainY, 'AdaBoostM1', 5, 'Tree', 'Holdout', 0.2 );
    %cc = fitensemble( trainX, trainY, 'Bag', 5, 'Tree', 'Holdout', 0.2, 'Type', 'classification' );
    looresult(k) = cc.Trained{1}.predict(X(~flag,:));
    disp( [k, Y(k), looresult(k)] );
end
fprintf( 'LOO error %f\n', sum(abs(looresult - Y')) / N );

%%
if strcmp(c1, 'man') && strcmp(c2, 'gal')
    N = size(manVectors.(ion),1);
    for k = 1 : N
        fprintf( '%d\t%d\t%s\n', Y(k), looresult(k), manData.(ion)(k).ions(1).mFormula );
    end
    for k = 1 : size(galVectors.B,1)
        fprintf( '%d\t%d\t%s\n', Y(k+N), looresult(k+N), galData.(ion)(k).ions(1).mFormula );
    end
elseif strcmp(c1, 'man') && strcmp(c2, 'glc')
    N = size(manVectors.(ion),1);
    for k = 1 : N
        fprintf( '%d\t%d\t%s\n', Y(k), looresult(k), manData.(ion)(k).ions(1).mFormula );
    end
    for k = 1 : size(glcVectors.B,1)
        fprintf( '%d\t%d\t%s\n', Y(k+N), looresult(k+N), glcData.(ion)(k).ions(1).mFormula );
    end
elseif strcmp(c1, 'gal') && strcmp(c2, 'glc')
    N = size(galVectors.(ion),1);
    for k = 1 : N
        fprintf( '%d\t%d\t%s\n', Y(k), looresult(k), galData.(ion)(k).ions(1).mFormula );
    end
    for k = 1 : size(glcVectors.B,1)
        fprintf( '%d\t%d\t%s\n', Y(k+N), looresult(k+N), glcData.(ion)(k).ions(1).mFormula );
    end
end
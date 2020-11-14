function ionClassifier = train_ionclassifier( trainVectors, trainData, do_LOO )

if nargin < 3
    do_LOO = 0;
end

warning('off', 'stats:classreg:learning:modifier:AdaBoostM1:modify:Terminate');

boostNum = 100;
holdOut = 0.20;
%posSet = {'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C', 'Y', 'Y', 'Y', 'Y', 'Z', 'Z', 'Z', 'Z'}; % 16
%negSet = {'C', 'Y', 'Z', 'O', 'B', 'Y', 'Z', 'O', 'B', 'C', 'Z', 'O', 'C', 'B', 'Y', 'O'};
posSet = {'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C'}; % 8
negSet = {'C', 'Y', 'Z', 'O', 'B', 'Y', 'Z', 'O'};
ionClassifier = [];
for round = 1 : length(posSet)

    posIon = posSet{round};
    negIon = negSet{round}; % no enough data, put them in one bucket.
    % use_simple_model = {'B_v_Y', 'B_v_Z', 'C_v_Y'};
    use_simple_model = {};

    p_v_n = [posIon, '_v_', negIon];
    disp( p_v_n );
    
    posX = []; negX = [];
    
    % find samples that appear in both posX and negX
    posSource = [];
    negSource = [];
    
    for n = posIon
        posX = [posX; trainVectors.(n)];
        posSource = [posSource; [[trainData.(n).glycanID]; [trainData.(n).peakID]]'];
    end
    
    for n = negIon
        negX = [negX; trainVectors.(n)];
        negSource = [negSource; [[trainData.(n).glycanID]; [trainData.(n).peakID]]'];
    end
    
    X = [posX; negX];
    Y = [ones(size(posX,1), 1); zeros(size(negX,1), 1)];
    XInfo = [posSource; negSource];
    weights = [ones(size(posX,1),1)/size(posX,1); ones(size(negX,1),1)/size(negX,1)] / 2;
    
    if do_LOO
        % Leave-one-out
        N = size(X,1);
        LOO = NaN(2, N);
        accP = 0; accN = 0; numP = 0; numN = 0;
        
        for k = 1 : N
            flag = ones(1, N);
            % flag(k) = 0;
            flag((XInfo(k,1) == XInfo(:,1)) & (XInfo(k,2) == XInfo(:,2))) = 0; % some peaks can be mutliple types (B, C, Y, Z)
            flag = flag > 0;
            if any(strcmp( p_v_n, use_simple_model ) )
                tree = fitctree(X(flag,:),Y(flag));
                [a, b] = tree.predict(X(k,:));
                LOO(1, k) = a;
                LOO(2, k) = b(a+1);
            else
                ClassTreeEns = fitensemble(X(flag,:), Y(flag), 'AdaBoostM1', boostNum, 'Tree', 'Holdout', holdOut, 'Weights', weights(flag));
                if isempty(ClassTreeEns.Trained{1}.predictorImportance)
                    tree = fitctree(X, Y);
                    [a, b] = tree.predict(X(k,:));
                    LOO(1, k) = a;
                    LOO(2, k) = b(a+1);
                else
                    [a, b] = ClassTreeEns.Trained{1}.predict(X(k,:));
                    LOO(1, k) = a;
                    LOO(2, k) = b(a+1);
                end
            end
            if Y(k) == 1
                accP = accP + (Y(k) == LOO(1, k));
                numP = numP + 1;
            else
                accN = accN + (Y(k) == LOO(1, k));
                numN = numN + 1;
            end
            msg = sprintf( '%s vs %s: (%4d,%4d)\t[%d --> %d (%5.2f)]\t%d/%d\t%d/%d\n', posIon, negIon, ...
                XInfo(k,1), XInfo(k,2), Y(k), LOO(1,k), LOO(2,k), accP, numP, accN, numN );
            fprintf( msg );
        end
        disp( [p_v_n, ' LOO finished.'] );
    end
    
    temp = [trainVectors.massFeatures, trainVectors.massFeatures, trainVectors.massFeatures];
    ionClassifier.(p_v_n).X = X;
    ionClassifier.(p_v_n).Y = Y;
    ionClassifier.(p_v_n).XInfo = XInfo;
    ionClassifier.(p_v_n).BoostingNum = boostNum;
    ionClassifier.(p_v_n).BoostingHoldout = holdOut;
    if do_LOO
        ionClassifier.(p_v_n).LOO = LOO;
        ionClassifier.(p_v_n).fail = XInfo( LOO(1,:)' ~= Y, : );
    end
    
    % remove shared samples
    [~, ia] = unique( XInfo, 'rows' );
    
    if any(strcmp( p_v_n, use_simple_model ) )
        ionClassifier.(p_v_n).classifier = fitctree(X(ia,:), Y(ia,:));
        % ionClassifier.(p_v_n).classifier = fitensemble(X(ia,:), Y(ia,:), 'Bag', boostNum, 'Tree', 'Holdout', holdOut, 'Weights', weights(ia));
        ionClassifier.(p_v_n).FeatureImportance = ionClassifier.(p_v_n).classifier.predictorImportance;
    else
        ionClassifier.(p_v_n).classifier = fitensemble(X(ia,:), Y(ia,:), 'AdaBoostM1', boostNum, 'Tree', 'Holdout', holdOut, 'Weights', weights(ia));
        if isempty( ionClassifier.(p_v_n).classifier.Trained{1}.predictorImportance )
            % ionClassifier.(p_v_n).classifier = fitctree(X(ia,:), Y(ia,:));
            ionClassifier.(p_v_n).classifier = fitensemble(X(ia,:), Y(ia,:), 'Bag', boostNum, 'Tree', ...
                               'Holdout', holdOut, 'Weights', weights(ia), 'Type', 'classification');
            % ionClassifier.(p_v_n).FeatureImportance = ionClassifier.(p_v_n).classifier.predictorImportance;
        % else
            % ionClassifier.(p_v_n).FeatureImportance = ionClassifier.(p_v_n).classifier.Trained{1}.predictorImportance;
            % ionClassifier.(p_v_n).SelectedFeatureIdx = find( ionClassifier.(p_v_n).FeatureImportance(1:end-2) > 0 );
            % ionClassifier.(p_v_n).SelectedFeatureMass = temp( ionClassifier.(p_v_n).SelectedFeatureIdx );
        end
        ionClassifier.(p_v_n).FeatureImportance = ionClassifier.(p_v_n).classifier.Trained{1}.predictorImportance;
        ionClassifier.(p_v_n).SelectedFeatureIdx = find( ionClassifier.(p_v_n).FeatureImportance(1:end-2) > 0 );
        ionClassifier.(p_v_n).SelectedFeatureMass = temp( ionClassifier.(p_v_n).SelectedFeatureIdx );
    end
    disp( ionClassifier.(p_v_n) );
    % disp( num2cell([LOOResult.(p_v_n).SelectedFeatureIdx; LOOResult.(p_v_n).SelectedFeatureMass]) );
    
end % round

warning('on', 'stats:classreg:learning:modifier:AdaBoostM1:modify:Terminate');

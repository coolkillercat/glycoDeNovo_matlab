classdef trainMethods < handle
    properties
        methodList = ['cTree', 'randomForest', 'AdaBoostM1'];
    end
    
    methods
    
        function c = train_by_method(X, Y, method, hyperparameters)
            c.classifier = [];
            c.FeatureImportance = [];
            if nargin < 3 
                warning('train_by_method error: insufficient parameters')
                return;
            end
            if nargin == 3
                mBoostNum = 20;
                weights = ones(size(X,1))/size(X,1);
            else
                mBoostNum = hyperparameters{1};
                weights = hyperparameters{2};
            end
            switch(method)
                case 'cTree'
                    c.classifier = fitctree(X, Y);
                    c.FeatureImportance = c.classifier.predictorImportance;
                case 'randomForest'
                    c.classifier = TreeBagger(mBoostNum, X, Y);
                    c.FeatureImportance = c.classifier.predictorImportance;
                case 'AdaBoostM1'
                    c.classifier = fitensemble(X, Y, 'AdaBoostM1', mBoostNum, 'tree', ...
                        'Weights', weights,'type', 'classification');
                    c.FeatureImportance = c.classifier.Trained{1}.predictorImportance;
                otherwise
                    warning('train_by_method error: no such method')
                    return;
            end
        end
    end

end
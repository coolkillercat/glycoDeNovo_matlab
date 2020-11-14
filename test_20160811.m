%% 2016/08
clear;

allpath = 'D:\Projects\Glycomics\data\20160811\new_format\';
alldirs = dir( allpath );
alldirs = alldirs(3:end);
accuracy = 0.02;

for d = 1 : length(alldirs)
    datapath = [allpath, alldirs(d).name, '\'];
    disp( '============================' );
    disp( datapath );
    temp = dir( [datapath, '*.txt'] );
    filenames = cell(1, length(temp));
    for k = 1 : length(temp)
        filenames{k} = temp(k).name;
    end
    gap = zeros(1, length(temp));
    minus2H = zeros(1, length(temp));

    for f = 1 : length( filenames )
        disp( ['Working on ' filenames{f}] );
        specU = CSpectrum.load( [datapath, filenames{f}] );
        specU.protonate( 1, accuracy );
        specU.merge_peaks( accuracy );
        
        reconstructor = CGlycanDeNovo( accuracy, [] );
        reconstructor.mCheckMinusH = 0;
        reconstructor.interpret_peaks( specU );
        reconstructor.reconstruct_formulas();
        
        if isempty(specU.mPeaks(end).mInferredFormulas)
            specU.add_complementary_ions( [], 0);
            specU.merge_peaks( accuracy );
            reconstructor.interpret_peaks( specU );
            reconstructor.reconstruct_formulas();
        end
        
        if isempty(specU.mPeaks(end).mInferredFormulas)
            disp( 'Empty reconstruction' );
        end
        
        name = ['rec.' filenames{f}];
        if ~exist( [datapath, 'results'], 'dir' )
            mkdir( [datapath, 'results'] );
        end
        specU.save_reconstruction( [datapath, 'results', filesep, name], [], reconstructor.mCheckMinus2H*2, reconstructor.mCheckGap );
        specU.filename = strrep( [datapath, 'results', filesep, name], '.txt', '.mat' );
        save( specU.filename, 'specU' );
        fprintf( '\n' );
    end
end

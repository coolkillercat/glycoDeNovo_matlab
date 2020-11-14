function batch_reconstruct( datapath, filenames, optionMinus2H, optionGap, accuracy, forceGap, addComplement, complementOnly )

if nargin < 5, accuracy = 5; end % 5ppm
if nargin < 6 || isempty(forceGap), forceGap = 0; end
if nargin < 7 || isempty(addComplement), addComplement = 1; end
if nargin < 8 || isempty(complementOnly), complementOnly = 0; end

for f = 1 : length( filenames )
    disp( ['Working on ' filenames{f}] );
    specU = CSpectrum.load( [datapath, filenames{f}] );
    specU.protonate( 1 );
    specU.merge_peaks( 0.001 );
    if addComplement
        specU.add_complementary_ions();
    end
    %specU.merge_peaks( 0.001 );
    
    reconstructor = CGlycoDeNovo( accuracy(1), [], optionMinus2H(f), optionGap(f) );
    reconstructor.mCheckMinusH = 0;
    if complementOnly
        reconstructor.mUseComplementOnly = 1;
    end
    reconstructor.interpret_peaks( specU );
    
    reconstructor.reconstruct_formulas();
    
    if isempty(specU.mPeaks(end).mInferredFormulas) && forceGap
        reconstructor.mCheckGap = 2;
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
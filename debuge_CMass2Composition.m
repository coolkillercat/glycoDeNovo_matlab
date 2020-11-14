spec = CSpectrum.load( 'D:\Projects\Glycomics\data\20190116\G2FB.1b28.txt' );
% spec = CSpectrum.load( 'D:\Projects\Glycomics\data\20190116\Man5.1a1.txt' );
% spec = CSpectrum.load( 'D:\Projects\Glycomics\data\20190116\G3S1F.1b34.txt' );
% spec = CSpectrum.load( 'D:\Projects\Glycomics\data\20170615\A2F.Reduced.Na.EED.PM.txt' );

spec.protonate();

m2c.set_reducing_end_modification( spec.mReducingEndModification );
m2c.set_permethylation( spec.mPermethylated );

%%
testS = m2c.correct_spectrum( spec );

%%
reconstructor = CGlycoDeNovo( 5, find(testS.mComposition > 0), 0, 0 );
reconstructor.mCheckMinusH = 0;
reconstructor.mMaxBranchingNum = 2;
reconstructor.mCompositionConstraint = testS.mComposition;
reconstructor.interpret_peaks( testS );

%%
tic;
minus2H = 0;
gap = 0;
for k = 1 : length(testS)    
    disp( ['========== ', num2str(k), ' =========='] );
    reconstructor = CGlycoDeNovo( 5, find(testS(k).mComposition > 0), minus2H, gap );
    reconstructor.mCheckMinusH = 0;
    reconstructor.mMaxBranchingNum = 3;
    reconstructor.mCompositionConstraint = testS(k).mComposition;
    reconstructor.interpret_peaks( testS(k) );
    reconstructor.reconstruct_formulas();
    
    %
    name = ['rec.test' num2str(k) '.txt'];
    testS(k).save_reconstruction( ['D:\Projects\Glycomics\matlab\results', filesep, name], [], reconstructor.mCheckMinus2H*2, reconstructor.mCheckGap );
    debug_reconstruction(testS(k), ['D:\Projects\Glycomics\matlab\results', filesep, strrep(name, '.txt', '.debug.txt')] );
    disp( name );
end
toc;
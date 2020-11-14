function annotate_glycosic_bonds( spectra, result_dir )

if isempty( spectra ), return; end

for k = 1 : length(spectra)
    spec = spectra{k};
    if strcmp(spec.mExperimentMethod, 'EED' ) == 0 || isempty( spec.comment )
        continue;
    end
    
    temp = strsplit( spec.filename, '/' );
    
    result_file = [result_dir, temp{end-1}, '_', strrep(temp{end}, '.txt', '.ann') ];

    g = CGlycan( spec.mPermethylated, spec.mReducingEndModification );
    idx = strfind( spec.comment, '.' );
    if isempty(idx)
        name = spec.comment;
    else
        name = strtrim( spec.comment(idx+1:end) );
    end
    try
        g.parse(name);
    catch
        disp( ['No glycan name. ', num2str(k), ': ', name] );
        continue;
    end
    
    % cleave g, keep B/C/Y/Z ions
    [Ys, Bs, Xs, As, Zs, Cs] = g.cleave();
    
    % Save annotation
    fid = fopen( result_file, 'w' );
    fprintf( fid, '# %s\n', name );
    fprintf( fid, '# Method: %s\n', spec.mExperimentMethod );
    fprintf( fid, '# Metal: Proton\n' );
    if spec.mPermethylated
        fprintf( fid, '# Permethylated\n' );
    end
    fprintf( fid, '# Precursor: H+ %f\n', spec.mPrecursor );
    fprintf( fid, '\nm/z\tz\tintensity\tpeaktype\tion-type\tRE\tNRE\tlinkage\tion-formula\n' );
        
    for m = 1 : length(spec.mPeaks)
        aPeak = spec.mPeaks(m);
        if aPeak.mIsComplement
            peaktype = 'com';
        else
            peaktype = '';
        end
        
        found = internal_save( fid, As, aPeak, peaktype );
        found = found + internal_save( fid, Bs, aPeak, peaktype );
        found = found + internal_save( fid, Cs, aPeak, peaktype );
        found = found + internal_save( fid, Xs, aPeak, peaktype );
        found = found + internal_save( fid, Ys, aPeak, peaktype );
        found = found + internal_save( fid, Zs, aPeak, peaktype );
        if found == 0
            fprintf( fid, '%f\t%d\t%10d\t%s\n', aPeak.mMass, 1, aPeak.mIntensity, peaktype );
        end
    end
    fclose(fid);    
end

    function found = internal_save( fileID, ions, onePeak, pType )
        found = 0;
        masses = [ions.mMass];
        d = abs( masses - onePeak.mMass );
        idxes = find( d < 0.005 );
        if ~isempty(idxes)
            for p = 1 : length(idxes)
                if d(idxes(p)) / masses(idxes(p)) < 0.000005
                    found = 1;
                    ion = ions(idxes(p));
                    fprintf( fileID, '%f\t%d\t%10d\t%s', onePeak.mMass, 1, onePeak.mIntensity, pType );                    
                    fprintf( fileID, '\t%s\t%s\t%s\t%d\t%s\n', ion.mType, ion.mReducingEndMonosaccharide.mClass, ...
                        ion.mNonReducingEndMonosaccharide.mClass, ion.mNonReducingEndMonosaccharide.mLinkedTo.C, ion.mFormula(3:end) );
                end
            end
        end
    end
end
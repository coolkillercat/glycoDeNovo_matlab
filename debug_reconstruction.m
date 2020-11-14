function debug_reconstruction( spec, filename, debugMatchedSpectrum )

if nargin < 3, debugMatchedSpectrum = 0; end
if isempty( spec.comment ), return; end

g = CGlycan(spec.mPermethylated);
g.mReducingEndModification = spec.mReducingEndModification;
g.parse( spec.comment );

if debugMatchedSpectrum && any(g.mComposition - spec.mComposition)
    return;
end
disp( ['Spectrum: ', spec.comment] );
disp( ['Glycan: ', g.mFormula] );
fprintf('\n');

if nargin < 2
    fid = -1;
else
    fid = fopen( filename, 'w' );
end

if fid > -1
    idx = strfind( spec.filename, 'data' );
    fprintf(fid, '%s\n', ['File: ', spec.filename(idx+5:end-4)] );
    fprintf(fid, '%s\n\n', ['Glycan: ', g.mFormula] );
end

if ~spec.mProtonated
    [~, Bs, ~, ~, ~, Cs] = g.cleave( spec.mMetal );
else
    [~, Bs, ~, ~, ~, Cs] = g.cleave();
end

massAccuracyPPM = 0.000005;

% pair B and C ions
temp = [Bs.mMass];
[~, idx] = sort(temp);
Bs = Bs(idx);
temp = [Bs.mMass];
[~, idx] = unique(temp);
Bs = Bs(idx);
ionpairs = cell(1, length(Bs));
for k = 1 : length(Bs)
    for cIon = Cs
        if strcmp( Bs(k).mConciseFormula, cIon.mConciseFormula )
            ionpairs{k} = [Bs(k), cIon];
            break;
        end
    end
end

peakmasses = [spec.mPeaks.mMass];
for k = 1 : length(Bs)
    aIon = ionpairs{k}(1);
    
    msg = ['Ion: ', aIon.mConciseFormula, ', B m/z=', num2str(aIon.mMass), ', C m/z=', num2str(ionpairs{k}(2).mMass)];
    disp( msg );
    if fid ~= -1
        fprintf(fid, '%s\n', msg);
    end
    
    massdiff = peakmasses - aIon.mMass;
    [mv, idx] = min( abs( massdiff ) );
    aPeak = spec.mPeaks(idx);
    if aPeak.mIsComplement
        if aPeak.mComplement == 0
            if isempty( aPeak.mOriPeaks )
                threshold = peakmasses(idx) * massAccuracyPPM;
            else
                threshold = aPeak.mOriPeaks(1).mMass * massAccuracyPPM;
            end
        else
            threshold = peakmasses( -aPeak.mComplement) * massAccuracyPPM;
        end
    else
        threshold = peakmasses(idx) * massAccuracyPPM;
    end
    if ( mv < threshold )
        msg = ['   B Match peak ', num2str(idx), ' m/z=', num2str(peakmasses(idx))];
        disp(msg);
        if fid ~= -1
            fprintf(fid, '%s\n', msg);
        end
    else
        massdiff2 = peakmasses - aIon.mMass + CMass.H*2;
        [mv2, idx2] = min( abs( massdiff2 ) );
        if ( mv2 < threshold )
            msg = ['   B Match-2H peak ', num2str(idx2), ' m/z=', num2str(peakmasses(idx2))];
            disp(msg);
            if fid ~= -1
                fprintf(fid, '%s\n', msg);
            end
        else
            msg = ['   B Miss (peak ', num2str(idx), ' m/z=', num2str(peakmasses(idx))];
            if massdiff(idx) < 0
                msg = [msg, ' by ', num2str(mv-threshold)];
            else
                msg = [msg, ' by -', num2str(mv-threshold)];
            end
            msg = [msg, ', -2H to peak ', num2str(idx2), ' m/z=', num2str(peakmasses(idx2)), ' by '];
            if massdiff2(idx2) < 0
                msg = [msg, num2str(mv2-threshold), '): '];
            else
                msg = [msg, '-', num2str(mv2-threshold), '): '];
            end
            disp(msg);
            if fid ~= -1
                fprintf(fid, '%s\n', msg);
            end
        end
    end
    
    aIon = ionpairs{k}(2);
    massdiff = peakmasses - aIon.mMass;
    [mv, idx] = min( abs( massdiff ) );
    if spec.mPeaks(idx).mIsComplement
        if aPeak.mComplement == 0
            if isempty( aPeak.mOriPeaks )
                threshold = peakmasses(idx) * massAccuracyPPM;
            else
                threshold = aPeak.mOriPeaks(1).mMass * massAccuracyPPM;
            end
        else
            threshold = peakmasses( -aPeak.mComplement) * massAccuracyPPM;
        end
    else
        threshold = peakmasses(idx) * massAccuracyPPM;
    end
    if ( mv < threshold )
        msg = ['   C Match peak ', num2str(idx), ' m/z=', num2str(peakmasses(idx))];
        disp(msg);
        if fid ~= -1
            fprintf(fid, '%s\n', msg);
        end
    else
        massdiff2 = peakmasses - aIon.mMass + CMass.H*2;
        [mv2, idx2] = min( abs( massdiff2 ) );
        if ( mv2 < threshold )
            msg = ['   C Match-2H peak ', num2str(idx2), ' m/z=', num2str(peakmasses(idx))];
            disp(msg);
            if fid ~= -1
                fprintf(fid, '%s\n', msg);
            end
        else
            msg = ['   C Miss (peak ', num2str(idx), ', m/z=', num2str(peakmasses(idx))];
            if massdiff(idx) < 0
                msg = [msg, ' by ', num2str(mv-threshold)];
            else
                msg = [msg, ' by -', num2str(mv-threshold)];
            end
            msg = [msg, ', -2H to peak ', num2str(idx2), ' m/z=', num2str(peakmasses(idx)), ' by '];
            if massdiff2(idx2) < 0
                msg = [msg, num2str(mv2-threshold), '): '];
            else
                msg = [msg, '-', num2str(mv2-threshold), '): '];
            end
            disp(msg);
            if fid ~= -1
                fprintf(fid, '%s\n', msg);
            end
        end
    end
end

if fid > -1
    fclose(fid);
end
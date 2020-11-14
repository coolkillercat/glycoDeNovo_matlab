classdef CGlycan < handle % By Pengyu Hong @ Brandeis University
    properties
        mStem = CMonosaccharide.empty(0.0);
        mBranches = CGlycan.empty(0.0);
        mFormula = '';
        mInferredFormula = '';
        mConciseFormula = '';
        mBranchingCode = '';
        mComposition = zeros(1, 8);

        mMass = 0;
        mNumUnits = 0;

        mPermethylated = 0;
        mIsCleaved = 0;
        
        mReducingEndModification = '';
    end
    
    methods(Static)
%         function massCom = reducing_end_mass_compensation( reducingEndModification, permethylated )
%         % function massCom = reducing_end_mass_compensation( reducingEndModification, permethylated )
%             if nargin < 2 || isempty(permethylated)
%                 permethylated = 0;
%             end
%             massCom = 0;
%             switch reducingEndModification
%                 case CMass.cReducingEndModification_Aminopyridine
%                     massCom = CMass.cMassCompensation_Aminopyridine;
%                 case CMass.cReducingEndModification_Deuterium
%                     massCom = CMass.cMassCompensation_Deuterium;
%                 case CMass.cReducingEndModification_O18
%                     massCom = CMass.cMassCompensation_O18;
%                 case CMass.cReducingEndModification_PRAGS
%                     massCom = CMass.cMassCompensation_PRAGS;
%                 case CMass.cReducingEndModification_Reduced
%                     massCom = CMass.H2 + permethylated * CMass.CH2;
%             end
%         end
        
        function [stem_formula, branch_formulas] = break_formula( formula )
            formula = strtrim( formula );
            
            % find branching points
            branch_pos = []; num = 0;            
            for k = 1 : length( formula )
                if formula(k) == '['
                    if num == 0
                        branch_pos(1,end+1) = k;
                    end
                    num = num + 1;
                elseif formula(k) == ']'
                    num = num - 1;
                    if num == 0
                        branch_pos(2,end) = k;
                    end
                end
            end
            if num ~= 0, throw( MExeption( 'CGlycan.parse', 'Wrong formula: mis-matched [-].') ); end
            
            % pre-assign the stem, will be changed later.
            if ~isempty( branch_pos )
                stem_formula = strtrim( formula( branch_pos(end,end)+1 : end) );
            else
                stem_formula = formula;
            end
            
            % get branch formula
            branch_formulas = cell( 1, size( branch_pos, 2 ) );
            if ~isempty( branch_pos )
                for k = 1 : size( branch_pos, 2 )
                    branch_formulas{k} = strtrim( formula( branch_pos(1, k)+1 : branch_pos(2, k)-1 ) );
                end
            end
        end
       
        function result = checkNLinked( glycan )
            result = 0;
            if isa(glycan, 'char')
                frml = glycan;
            elseif isa(glycan, 'CGlycan')
                frml = glycan.mConciseFormula;
            end
            
            [stem_formula, branch_formulas] = CGlycan.break_formula( frml );
            if length( branch_formulas ) > 2
                return;
            end
            
            if strcmp( stem_formula, 'Hex HexNAc HexNAc' ) && length(branch_formulas) == 2
                if length(branch_formulas{1}) > 4 && strcmp( branch_formulas{1}(end-3:end), ' Hex') && ...
                   length(branch_formulas{2}) > 4 && strcmp( branch_formulas{2}(end-3:end), ' Hex')
                    result = 1;
                end
            elseif strcmp( stem_formula, 'HexNAc HexNAc' )
                if strcmp( branch_formulas{1}, 'Fuc' )
                    [stem_formula, branch_formulas] = CGlycan.break_formula( branch_formulas{2} );
                elseif strcmp( branch_formulas{2}, 'Fuc' )
                    [stem_formula, branch_formulas] = CGlycan.break_formula( branch_formulas{1} );
                else
                    return;
                end
                if strcmp( stem_formula, 'Hex' ) && length(branch_formulas) == 2 && ...
                   length(branch_formulas{1}) > 4 && strcmp( branch_formulas{1}(end-3:end), ' Hex') && ...
                   length(branch_formulas{2}) > 4 && strcmp( branch_formulas{2}(end-3:end), ' Hex')
                    result = 1;
                end
            elseif strcmp( stem_formula, 'HexNAc' )
                if strcmp( branch_formulas{1}, 'Fuc' )
                    [stem_formula, branch_formulas] = CGlycan.break_formula( branch_formulas{2} );
                elseif strcmp( branch_formulas{2}, 'Fuc' )
                    [stem_formula, branch_formulas] = CGlycan.break_formula( branch_formulas{1} );
                else
                    return;
                end
                if strcmp( stem_formula, 'Hex HexNAc' ) && length(branch_formulas) == 2 && ...
                   (strcmp(branch_formulas{1}, 'Hex') || (length(branch_formulas{1}) > 4 && strcmp( branch_formulas{1}(end-3:end), ' Hex'))) && ...
                   (strcmp(branch_formulas{1}, 'Hex') || (length(branch_formulas{2}) > 4 && strcmp( branch_formulas{2}(end-3:end), ' Hex')))
                    result = 1;
                end
            end
        end
    end
    
    methods
        function obj = CGlycan( permethylation, reducingEndModification )
        % function obj = CGlycan( permethylation, reducingEndModification )
            if nargin >= 1
                if ~isempty( permethylation )
                    obj.mPermethylated = permethylation;
                end
                if nargin >= 2 && ~isempty( reducingEndModification )
                    obj.mReducingEndModification = reducingEndModification;
                end
            end
        end
        
        function copy_modifications( obj, sourceObj )
            obj.mPermethylated = sourceObj.mPermethylated;
            obj.mReducingEndModification = sourceObj.mReducingEndModification;
        end
        
        function [monos, childrenMap] = flatten( obj, monos, childrenMap )
        % convert the tree structure of a glycan into a falt structure
        % stored in monos and childrenMap.
            if nargin < 2, monos = CMonosaccharide.empty( obj.mNumUnits, 0 ); end
            if nargin < 3, childrenMap = cell(1, obj.mNumUnits); end
            
            for k = 1 : length( obj.mStem )
                monos( obj.mStem(k).mID ) = obj.mStem(k);
            end
            for k = 1 : length( obj.mStem )-1
                childrenMap{ obj.mStem(k).mID } = obj.mStem(k+1).mID;
            end
            for k = 1 : length( obj.mBranches )
                childrenMap{ obj.mStem(end).mID } = [childrenMap{ obj.mStem(end).mID }, obj.mBranches(k).mStem(1).mID];
                [monos, childrenMap] = obj.mBranches(k).flatten( monos, childrenMap );
            end
        end
        
        function [Ys, Bs, Xs, As, Zs, Cs] = cleave( obj, ionMetal )
            if nargin < 2, ionMetal = 'H'; end
            Ys = []; Bs = []; Xs = []; As = []; Zs = []; Cs = [];
            
            if isa( obj, 'CIon' )
                towardABC = sum( obj.type == 'A' ) > 0 || sum( obj.type == 'B' ) > 0 || sum( obj.type == 'C' ) > 0;
                towardXYZ = sum( obj.type == 'Y' ) > 0 || sum( obj.type == 'X' ) > 0 || sum( obj.type == 'Z' ) > 0;
            else
                towardXYZ = 0;
                towardABC = 0;
            end
            
            % Cleave obj.mStem(1:end-1) into Y, B, X, and A.
            for k = 1 : length( obj.mStem )-1
                aY = CIon( 'Y', ionMetal );
                if towardABC
                    aY.type = [aY.type, obj.type];
                end
                aY.mIsCleaved = 1;
                aY.copy_modifications( obj );
                aY.mStem = obj.mStem(1:k).copy(); 
                aY.mCleaveMono = aY.mStem(k);
                aY.mReducingEndMonosaccharide = aY.mCleaveMono;
                aY.mNonReducingEndMonosaccharide = obj.mStem(k+1);
                aY.update( 1 );
                Ys = [Ys, aY];
                
                aZ = aY.copy(); 
                aZ.mType(1) = 'Z'; 
                aZ.update( 1 );
                Zs = [Zs, aZ];
                
                aB = CIon( 'B', ionMetal );
                if towardXYZ
                    aB.type = [aB.type, obj.type];
                end
                aB.mStem = obj.mStem(k+1:end).copy();
                aB.mCleaveMono = obj.mStem(k);
                if ~isempty( obj.mBranches )
                    aB.set_branches( obj.mBranches.copy() );
                end
                aB.mIsCleaved = 1;
                aB.mReducingEndMonosaccharide = aB.mCleaveMono;
                aB.mNonReducingEndMonosaccharide = obj.mStem(k+1);
                % aY.copy_modifications( obj );
                aB.update( 1 );
                Bs = [Bs, aB];
                
                aC = aB.copy(); 
                aC.mType(1) = 'C'; 
                aC.update( 1 );
                aC.mCleaveMono = obj.mStem(k);
                Cs = [Cs, aC];
                
                if isempty(obj.mStem(k).mCleavage.type) 
                    [tXs, tAs] = obj.XA_cleave_nonend_stem_unit( k, ionMetal );
                    for tion = tXs, tion.mIsCleaved = 1; end
                    for tion = tAs, tion.mIsCleaved = 1; end
                    if towardXYZ || towardABC
                        sbV = [];
                        for m = 3 : 6
                            if ~isempty(obj.mStem(k).mVertices.modification{m})
                                sbV = [sbV, m];
                            end
                        end
                        for tion = tXs
                            if tion.mNumUnits > 1 && ~isempty( intersect(sbV, tion.mStem(end).mCleavage.vertices) )
                                tion.type = [tion.type, obj.type];
                                tion.update( 1 );
                            end
                        end
                        for tion = tAs
                            if tion.mNumUnits > 1
                                tion.type = [tion.type, obj.type];
                                tion.update( 1 );
                            end
                        end
                    end
                    Xs = [Xs, tXs];
                    As = [As, tAs];
                end
            end
            
            % Cleave obj.mStem(end) into X and A
            if isa(obj, 'CIon')
                if isempty( obj.mStem(end).mCleavage.type )
                    sbV = [];
                    for k = 3 : 6
                        if ~isempty(obj.mStem(end).mVertices.modification{k})
                            sbV = [sbV, k];
                        end
                    end
                    [tXs, tAs] = obj.XA_cleave_end_stem_unit( ionMetal );
                    
                    if ~isempty( sbV )
                        flag = zeros(1, length(tXs));
                        for k = 1 : length(tXs)
                            tion = tXs(k);
                            if ~isempty( intersect( tion.mCleaveMono.mCleavage.vertices, sbV ) )
                                tion.type = [tion.type, obj.type];
                                tion.update( 1 );
                                flag(k) = 1;
                            end
                        end
                        tXs = tXs(flag>0);
                    end
                    for tion = tAs
                        onCleavedBranch = 0;
                        for br = tion.mBranches
                            if br.mIsCleaved
                                onCleavedBranch = 1;
                                break;
                            end
                        end
                        if onCleavedBranch
                            tion.type = [tion.type, obj.type];
                            tion.update( 1 );
                        end
                    end
                    Xs = [Xs, tXs];
                    As = [As, tAs];
                end
            else
                [tXs, tAs] = obj.XA_cleave_end_stem_unit( ionMetal );
                Xs = [Xs, tXs];
                As = [As, tAs];
            end
            
            % Cleave whole branches
            [tYs, tBs, tZs, tCs] = obj.cleave_whole_branches( ionMetal );
            Ys = [Ys, tYs]; 
            Bs = [Bs, tBs]; 
            Zs = [Zs, tZs]; 
            Cs = [Cs, tCs];             
            
            % Cleave individual branches
            num = 0;
            for k = 1 : length( obj.mStem(end).mVertices.modification )
                if ~isempty( obj.mStem(end).mVertices.modification{k} ) && ...
                   ( obj.mStem(end).mVertices.modification{k}(1) == '-' )
                    num = num + 1;
                end
            end
            branchCleaved = (num > length(obj.mBranches)+1);
            for k = 1 : length( obj.mBranches )
                % cleave the k-th branch
                [bYs, bBs, bXs, bAs, bZs, bCs] = obj.mBranches(k).cleave( ionMetal );
                for tion = bBs, tion.mIsCleaved = 1; end % this branch is cleaved
                for tion = bAs, tion.mIsCleaved = 1; end % this branch is cleaved
                for tion = bCs, tion.mIsCleaved = 1; end % this branch is cleaved
                
                if obj.mBranches(k).mIsCleaved
                    for tion = bBs
                        tion.type = [tion.typ, obj.type];
                        tion.update( 1 );
                    end
                    for tion = bAs
                        tion.type = [tion.typ, obj.type];
                        tion.update( 1 );
                    end
                    for tion = bCs
                        tion.type = [tion.typ, obj.type];
                        tion.update( 1 );
                    end
                end
                
                Bs = [Bs, bBs]; 
                As = [As, bAs]; 
                Cs = [Cs, bCs];
                
                for m = 1 : length( bXs )
                    aX = CIon( 'X', ionMetal, obj );
                    aX.mIsCleaved = 1;
                    if towardABC || branchCleaved
                        aX.type = [aX.type, obj.type];
                    end
                    aX.mBranches(k).mStem = bXs(m).mStem;
                    aX.mBranches(k).mBranches = bXs(m).mBranches; % don't use set_branches because it will change stem(end)'s vertex modification
                    aX.mBranches(k).update(1);
                    aX.mCleaveMono = bXs(m).mCleaveMono;
                    aX.mReducingEndMonosaccharide = bXs(m).mReducingEndMonosaccharide;
                    aX.mNonReducingEndMonosaccharide = bXs(m).mNonReducingEndMonosaccharide;
                    aX.update(1);
                    aX.mCleaveMono = bXs(m).mCleaveMono;
                    Xs = [Xs, aX];
                end
                for m = 1 : length( bYs )
                    aY = CIon( 'Y', ionMetal, obj );
                    aY.mIsCleaved = 1;
                    if towardABC || branchCleaved
                        aY.type = [aY.type, obj.type];
                    end
                    aY.mBranches(k).mStem = bYs(m).mStem;
                    aY.mBranches(k).mBranches = bYs(m).mBranches; % don't use set_branches because it will change stem(end)'s vertex modification
                    aY.mBranches(k).update(1);
                    aY.mCleaveMono = bYs(m).mCleaveMono;
                    aY.mReducingEndMonosaccharide = bYs(m).mReducingEndMonosaccharide;
                    aY.mNonReducingEndMonosaccharide = bYs(m).mNonReducingEndMonosaccharide;
                    aY.update();
                    Ys = [Ys, aY];
                end
                for m = 1 : length( bZs )
                    aZ = CIon( 'Z', ionMetal, obj );
                    aZ.mIsCleaved = 1;
                    if towardABC || branchCleaved
                        aZ.type = [aZ.type, obj.type];
                    end
                    aZ.mBranches(k).mStem = bZs(m).mStem;
                    aZ.mBranches(k).mBranches = bZs(m).mBranches; % don't use set_branches because it will change stem(end)'s vertex modification
                    aZ.mBranches(k).update(1);
                    aZ.mCleaveMono = bZs(m).mCleaveMono;
                    aZ.mReducingEndMonosaccharide = bZs(m).mReducingEndMonosaccharide;
                    aZ.mNonReducingEndMonosaccharide = bZs(m).mNonReducingEndMonosaccharide;
                    aZ.update( 1 );
                    Zs = [Zs, aZ];
                end
            end
            
            % Remove Redundant Ions
            if ( obj.mStem(1).mID == 1 )
                As = CIon.unique( As );
                Bs = CIon.unique( Bs );
                Cs = CIon.unique( Cs );
                Xs = CIon.unique( Xs );
                Ys = CIon.unique( Ys );
                Zs = CIon.unique( Zs );
            end
        end

        function [As, Bs, Cs, Xs, Ys, Zs] = cleave_unit( obj, unitID, ionMetal ) % Cleave the unit indicated by unitID.
            if nargin < 3, ionMetal = ''; end
            Xs = []; Ys = []; Zs = []; As = []; Bs = []; Cs = [];
            
            % Cleave obj.mStem(1:end-1) into X, and A.
            for k = 1 : length( obj.mStem )-1
                if obj.mStem(k).mID == unitID
                    Ys = CIon( 'Y', ionMetal );
                    Ys.copy_modifications( obj );
                    Ys.mStem = obj.mStem(1:k).copy(); 
                    Ys.update(1);
                    Zs = Ys.copy(); 
                    Zs.type = 'Z'; 
                    Zs.update(1);
                    
                    Bs = CIon( 'B', ionMetal );
                    Bs.mStem = obj.mStem(k+1:end).copy();
                    Bs.set_branches( obj.mBranches.copy() );
                    Bs.copy_modifications( obj );
                    Bs.update(1);
                    Cs = Bs.copy(); 
                    Cs.type = 'C'; 
                    Cs.update(1);
                    
                    [Xs, As] = obj.XA_cleave_nonend_stem_unit( k, ionMetal );
                    return;
                end
            end
            
            % Cleave obj.mStem(end) into X and A
            if obj.mStem(end).mID == unitID
                [Ys, Bs, Zs, Cs] = obj.cleave_whole_branches( ionMetal );
                [Xs, As] = obj.XA_cleave_end_stem_unit( ionMetal );
                return;
            end
            
            % Cleave obj.mBranches
            for k = 1 : length( obj.mBranches )
                [As, Bs, Cs, Xs, Ys, Zs] = obj.mBranches(k).cleave_unit( unitID, ionMetal );
                if ~isempty( Xs )
                    for m = 1 : length( Ys )
                        aY = CIon( 'Y', ionMetal, obj );
                        aY.mBranches(k).mStem = Ys(m).mStem;
                        aY.mBranches(k).mBranches = Ys(m).mBranches;
                        aY.mBranches(k).update();
                        aY.mCleaveMono = Ys(m).mCleaveMono;
                        aY.update(1);
                        Ys(m) = aY;
                        aZ = aY.copy(); aZ.type = 'Z'; aZ.update(1);
                        Zs(m) = aZ;
                    end
                    for m = 1 : length( Xs )
                        aX = CIon( 'X', ionMetal, obj );
                        aX.mBranches(k).mStem = Xs(m).mStem;
                        aX.mBranches(k).mBranches = Xs(m).mBranches;
                        aX.mBranches(k).update();
                        aX.mCleaveMono = Xs(m).mCleaveMono;
                        aX.update(1);
                        Xs(m) = aX;
                    end
                    return;
                end
            end
        end
               
        function [Xs, As] = XA_cleave_nonend_stem_unit( obj, dth_stem, ionMetal )
            if nargin < 3, ionMetal = ''; end
            if dth_stem == length( obj.mStem ) % obj.mStem(dth_stem) should not be a branching unit.
                throw( MException('CGlycan:XA_cleave_nonend_stem_unit', 'Tried to cleave a branching unit. Use XA_cleave_stem_end() instead.') );
            end
            if isempty( obj.mStem(dth_stem).mLinkedIn.V )
                throw( MException('CGlycan:XA_cleave_nonend_stem_unit', 'No linkedIn information.') );
            end

            Xs = []; As = [];
            if ~isempty( obj.mStem(dth_stem).mCleavage.type )  % has been cleaved
                 return;
            end
            
            for m = 1 : length( CCrossRingCleavage.cCRC )
                [x, a] = obj.mStem(dth_stem).cleave( CCrossRingCleavage.cCRC{m} );
                
                aX = CIon('X', ionMetal);
                aX.copy_modifications( obj );
                aX.mReducingEndMonosaccharide = obj.mStem(dth_stem);
                aX.mNonReducingEndMonosaccharide = obj.mStem(dth_stem);
                if dth_stem == 1
                    aX.mStem = x;
                else
                    % aX.mStem = obj.mStem(1:dth_stem-1).copy;
                    aX.mStem = obj.mStem(1:dth_stem-1);
                    aX.mStem(end+1) = x;
                end
                aX.mCleaveMono = x;
                if sum( aX.mStem(end).mCleavage.vertices == obj.mStem(dth_stem+1).mLinkedTo.V ) > 0
                    % aX.mStem = [aX.mStem, obj.mStem(dth_stem+1:end).copy];
                    aX.mStem = [aX.mStem, obj.mStem(dth_stem+1:end)];
                    if ~isempty( obj.mBranches )
                        % aX.mBranches = obj.mBranches.copy;
                        aX.set_branches( obj.mBranches );
                    end
                else
                    x.set_linkedIn([], []);
                end
                aX.update(1); 
                Xs = [Xs, aX];
                
                aA = CIon( 'A', ionMetal );
                aA.copy_modifications( obj );
                aA.mReducingEndMonosaccharide = obj.mStem(dth_stem);
                aA.mNonReducingEndMonosaccharide = obj.mStem(dth_stem);
                aA.mStem(1) = a;
                aA.mCleaveMono = a;
                if sum( aA.mStem(1).mCleavage.vertices == obj.mStem(dth_stem+1).mLinkedTo.V ) > 0
                    % aA.mStem = [aA.mStem, obj.mStem(dth_stem+1:end).copy];
                    aA.mStem = [aA.mStem, obj.mStem(dth_stem+1:end)];
                    if ~isempty( obj.mBranches )
                        % aA.mBranches = obj.mBranches.copy;
                        aA.set_branches( obj.mBranches );
                    end
                end
                aA.update(1); As = [As, aA];
            end
        end

        function [Xs, As] = XA_cleave_end_stem_unit( obj, ionMetal )
            if nargin < 2, ionMetal = ''; end
            Xs = []; As = [];
            if ~isempty( obj.mStem(end).mCleavage.type )  % has been cleaved
                return;
            end
            
            for m = 1 : length( CCrossRingCleavage.cCRC )
                [x, a] = obj.mStem(end).cleave( CCrossRingCleavage.cCRC{m} );
                % x.set_linkedIn( [], [] );
                a.mLinkedTo.V = []; 
                a.mLinkedTo.type = [];
                
                aX = CIon('X', ionMetal);
                aX.copy_modifications( obj );
                % aX.mStem = obj.mStem(1:end-1).copy;
                aX.mStem = [obj.mStem(1:end-1), x];
                aX.mCleaveMono = x;
                aX.mReducingEndMonosaccharide = obj.mStem(end);
                aX.mNonReducingEndMonosaccharide = obj.mStem(end);
                for b = 1 : length( obj.mBranches )
                    if sum( aX.mStem(end).mCleavage.vertices == obj.mBranches(b).mStem(1).mLinkedTo.V ) > 0
                        aX.mStem(end).remove_linkedIn( obj.mBranches(b).mStem(1).mLinkedTo.V );
                        % aX.add_branches( obj.mBranches(b).copy );
                        aX.add_branches( obj.mBranches(b) );
                    end
                end
                aX.update(1);
                Xs = [Xs, aX];
                
                aA = CIon('A', ionMetal);
                a.mLinkedIn.V = []; a.mLinkedIn.C = []; a.mLinkedIn.type = []; a.mLinkedIn.unitID = [];
                aA.mStem = a;
                aA.mCleaveMono = a;
                aA.copy_modifications( obj );
                aA.mReducingEndMonosaccharide = obj.mStem(end);
                aA.mNonReducingEndMonosaccharide = obj.mStem(end);
                if ~isempty( obj.mBranches ) % remove natural loss
                    for b = 1 : length( obj.mBranches )
                        if sum( aA.mStem(end).mCleavage.vertices == obj.mBranches(b).mStem(1).mLinkedTo.V ) > 0
                            % aA.add_branches( obj.mBranches(b).copy );
                            aA.add_branches( obj.mBranches(b) );
                        end
                    end
                    aA.update(1); 
                    As = [As, aA];
                end
            end
        end
        
        function [Ys, Bs, Zs, Cs] = cleave_whole_branches( obj, ionMetal )
            if nargin < 2, ionMetal = ''; end
            Ys = []; Bs = []; Zs = []; Cs = [];
            
            num = 0;
            for k = 1 : length( obj.mStem(end).mVertices.modification )
                if ~isempty( obj.mStem(end).mVertices.modification{k} ) && ...
                   ( obj.mStem(end).mVertices.modification{k}(1) == '-' )     
                    num = num + 1;
                end
            end
            branchCleaved = (num > length(obj.mBranches)+1);
            abc = isa( obj, 'CIon' ) && ( sum(obj.type == 'B') || sum(obj.type == 'C') || sum(obj.type == 'A') );
            
            for k = 1 : length( obj.mBranches )
                aY = CIon( 'Y', ionMetal );
                if branchCleaved || abc
                    aY.type = [aY.type, obj.type];
                end
                aY.copy_modifications( obj );
                aY.mStem = [obj.mStem(1:end-1), obj.mStem(end).copy()]; % again Replicate the monosaccharide that will be changed.
                aY.mCleaveMono = aY.mStem(end);
                % aY.mBranches = obj.mBranches( [1:k-1, k+1:end] ).copy();
                % aY.mBranches = obj.mBranches( [1:k-1, k+1:end] );
                aY.set_branches( obj.mBranches( [1:k-1, k+1:end] ) );
                if obj.mStem(end).mPermethylated % left CH2 scars
                    aY.mStem(end).mVertices.modification{ obj.mBranches(k).mStem(1).mLinkedTo.V } = '-CH2';
                    aY.mStem(end).update();
                end
                aY.mReducingEndMonosaccharide = obj.mStem(end);
                aY.mNonReducingEndMonosaccharide = obj.mBranches(k).mStem(1);
                aY.update(1);
                Ys = [Ys, aY];
                aZ = aY.copy;
                aZ.mType(1) = 'Z';
                aZ.update(1);
                Zs = [Zs, aZ];
                
                aB = CIon( 'B', ionMetal, obj.mBranches(k) );
                if obj.mBranches(k).mIsCleaved
                    aB.type = ['B', obj.type];
                end
                aB.mIsCleaved = 1;
                aB.copy_modifications( obj );
                aB.mStem(1) = aB.mStem(1).copy; % again Replicate the monosaccharide that will be changed.
                aB.mStem(1).mLinkedTo.V = [];
                aB.mStem(1).mLinkedTo.C = [];
                aB.mStem(1).mLinkedTo.type = [];
                aB.mStem(1).mLinkedTo.toUnitID = [];
                aB.mReducingEndMonosaccharide = obj.mStem(end);
                aB.mNonReducingEndMonosaccharide = obj.mBranches(k).mStem(1);
                aB.mCleaveMono = aY.mCleaveMono;
                aB.update(1);
                Bs = [Bs, aB];
                aC = aB.copy; 
                aC.mType(1) = 'C'; 
                aC.update(1);
                Cs = [Cs, aC];
            end
        end

        function order_branch_by_name( this )
            num = length(this.mBranches);
            for k = 1 : num
                this.mBranches(k).order_branch_by_name();
            end
            this.update_formula(1);
        end
        
        function results = simulate_spectrum( obj, ionMetal, ionTypes, fragmentTwice, noise )
            if nargin < 2, ionMetal = ''; end
            if nargin < 3 || isempty( ionTypes ), ionTypes = { 'Y', 'B', 'X', 'A', 'Z', 'C' }; end
            if nargin < 4 || isempty( fragmentTwice ), fragmentTwice = 0; end
            if nargin < 5 || isempty( noise ), noise = 0; end
            
            results = CSpectrumPeak.empty(1,0);

            [y_ions, b_ions, x_ions, a_ions, z_ions, c_ions] = obj.cleave( ionMetal );
            
            if sum( strcmp( ionTypes, 'Y' ) ) > 0
                for k = 1 : length( y_ions )
                    results(end+1) = CSpectrumPeak();
                    results(end).ions = y_ions(k);
                    results(end).types = {'Y'};
                    results(end).mMass = y_ions(k).mMass;
                end
            end
            
            if sum( strcmp( ionTypes, 'B' ) ) > 0
                for k = 1 : length( b_ions )
                    results(end+1) = CSpectrumPeak();
                    results(end).ions = b_ions(k);
                    results(end).types = {'B'};
                    results(end).mMass = b_ions(k).mMass;
                end
            end
            
            if sum( strcmp( ionTypes, 'X' ) ) > 0
                for k = 1 : length( x_ions )
                    results(end+1) = CSpectrumPeak();
                    results(end).ions = x_ions(k);
                    results(end).types = {'X'};
                    results(end).mMass = x_ions(k).mMass;
                end
            end
            
            if sum( strcmp( ionTypes, 'A' ) ) > 0
                for k = 1 : length( a_ions )
                    results(end+1) = CSpectrumPeak();
                    results(end).ions = a_ions(k);
                    results(end).types = {'A'};
                    results(end).mMass = a_ions(k).mMass;
                end
            end
            
            if sum( strcmp( ionTypes, 'Z' ) ) > 0
                for k = 1 : length( z_ions )
                    results(end+1) = CSpectrumPeak();
                    results(end).ions = z_ions(k);
                    results(end).types = {'Z'};
                    results(end).mMass = z_ions(k).mMass;
                end
            end
            
            if sum( strcmp( ionTypes, 'C' ) ) > 0
                for k = 1 : length( c_ions )
                    results(end+1) = CSpectrumPeak();
                    results(end).ions = c_ions(k);
                    results(end).types = {'C'};
                    results(end).mMass = c_ions(k).mMass;
                end
            end
            
            total = CIon('T', ionMetal, obj);
            results(end+1) = CSpectrumPeak();
            results(end).ions = total;
            results(end).types = {'T'};
            results(end).mMass = total.mMass;
            
            % merge ions with the same mass and remove unstable ions (less than a unit)
            flag = ones(1, length( results ));
            for k = 1 : length( results ) - 2
                if flag(k) == 0, continue; end
                for m = k + 1 : length( results ) - 1
                    if flag(m) == 0, continue; end
                    if results(m).ions(1).mNumUnits == 1 && ( results(m).ions(1).type == 'X' || results(m).ions(1).type == 'A' )
                        flag(m) = 0;
                        continue;
                    end
                    if abs( results(k).mMass - results(m).mMass ) < 0.0001
                        results(k).ions = [results(k).ions, results(m).ions];
                        results(k).types = [results(k).types, results(m).types];
                        flag(m) = 0;
                    end
                end
            end
            results = results(flag > 0);
            results = results.sort();
            
            % remove identical ions
            for k = 1 : length( results ) - 1
                num = length( results(k).ions );
                flag = ones(1, num);
                for m = 1 : num-1
                    if flag(m) == 0, continue; end
                    for n = m+1 : num
                        if strcmp( results(k).ions(m).mFormula, results(k).ions(n).mFormula )
                            flag(n) = 0;
                        end
                    end
                end
                results(k).ions = results(k).ions(flag>0);
                results(k).types = results(k).types(flag>0);
            end
            
            if fragmentTwice > 0
                internalIons = [];
                
                for k = 1 : length(results)-1 % don't touch the last one, it is the precursor.
                    for ion = results(k).ions
                        if ion.mNumUnits == 1, continue; end
                        [iys, ibs, ixs, ias, izs, ics] = ion.cleave( ionMetal );
                        temp = [iys, ibs, ixs, ias, izs, ics];
                        
                        % remove 'Y', 'X', 'Z', 'A', 'B', 'C'
                        flag = zeros(1, length(temp));
                        for m = 1 : length(temp)
                            if length(temp(m).type) > 1
                                if temp(m).mNumUnits == 1
                                    if sum( temp(m).type == 'A' ) == 0 && sum( temp(m).type == 'X' ) == 0
                                        flag(m) = 1;
                                    end
                                elseif temp(m).mNumUnits == 2
                                    if ~( strcmp(temp(m).type, 'AX') || strcmp(temp(m).type, 'XA') || ...
                                          strcmp(temp(m).type, 'AA')  || strcmp(temp(m).type, 'XX') )
                                        flag(m) = 1;
                                    end
                                else
                                    flag(m) = 1;
                                end
                            end
                        end
                        internalIons = [internalIons, temp(flag>0)];
                    end
                end
                
                % remove identical internal ions
                num = length( internalIons );
                flag = ones(1, num);
                for k = 1 : num-1
                    if flag(k)==0, continue; end
                    for m = k+1 : num
                        if flag(m) == 0, continue; end
                        if internalIons(k).mMass ~= internalIons(m).mMass
                            continue;
                        end
                        if strcmp(internalIons(k).type, internalIons(m).type) == 0
                            continue;
                        end
                        if strcmp( internalIons(k).mFormula, internalIons(m).mFormula )
                            flag(m) = 0;
                        end
                    end
                end
                internalIons = internalIons(flag>0);
                
                for k = 1 : length( internalIons )
                    idx = find( abs( internalIons(k).mMass - [results.mMass]) < 0.0001 );
                    if isempty( idx )
                        results(end+1) = CSpectrumPeak();
                        results(end).ions = internalIons(k);
                        results(end).types = {internalIons(k).type};
                        results(end).mMass = internalIons(k).mMass;
                    else
                        results(idx(1)).ions = [results(idx(1)).ions, internalIons(k)];
                        results(idx(1)).types = [results(idx(1)).types, {internalIons(k).type}];
                    end
                end
                
                % order results
                results = results.sort();
            end
            
            % add uncertainty
            if noise > 0
                for k = 1 : length( results )
                    results(k).mMass = results(k).mMass + noise * (rand - 0.5);
                end
            end
        end
        
        function update( obj, fast_mode )
        % If fast_mode == 1, do not update its branches.
            if isempty( obj ), return; end
            if nargin < 2, fast_mode = 0; end
            
            obj.mNumUnits = length( obj.mStem );
            obj.mMass = 0;
            
            obj.update_formula( fast_mode );
            
            if ~isempty(obj.mStem) && obj.mStem(1).mID == 1                
                obj.mMass = obj.mMass + CMass.get_mass_compensation( obj.mReducingEndModification, obj.mPermethylated );
            end
        end
           
        function update_formula( obj, fast_mode )
        % If fast_mode == 1, do not update its branches.
            obj.mComposition = zeros(1, 8);
            newFormula = '';
            newInferredFormula = '';
            newConciseFormula  = '';
            for k = length( obj.mStem ) : -1 : 1
                obj.mComposition( obj.mStem(k).mClassID ) = obj.mComposition( obj.mStem(k).mClassID ) + 1;
                if fast_mode == 0
                    obj.mStem(k).update();
                end
                if k < length( obj.mStem )
                    newFormula = [newFormula, ' '];
                    newInferredFormula = [newInferredFormula, ' '];
                    newConciseFormula = [newConciseFormula, ' '];
                end
                newFormula =[newFormula, obj.mStem(k).mFormula];
                newInferredFormula = [newInferredFormula, obj.mStem(k).mInferredFormula];
                newConciseFormula = [newConciseFormula, obj.mStem(k).mConciseFormula];
                obj.mMass = obj.mMass + obj.mStem(k).mMass;
            end
            obj.mMass = obj.mMass - CMass.H2O * (length( obj.mStem ) - 1); %lose a H2O when a glycosidic bond is formed.
            
            for k = 1 : length( obj.mBranches )
                if isempty( obj.mBranches(k) )
                    newFormula = ['[] ', newFormula];
                    newInferredFormula = ['[] ', newInferredFormula];
                    newConciseFormula = ['[] ', newConciseFormula];
                else
                    if fast_mode == 0
                        obj.mBranches(k).update();
                    end
                    newFormula = ['[', obj.mBranches(k).mFormula, '] ', newFormula];
                    newInferredFormula = ['[', obj.mBranches(k).mInferredFormula, '] ', newInferredFormula];
                    newConciseFormula = ['[', obj.mBranches(k).mConciseFormula, '] ', newConciseFormula];
                    obj.mMass = obj.mMass + obj.mBranches(k).mMass - CMass.H2O; %lose a H2O when a glycosidic bond is formed.
                    obj.mComposition = obj.mComposition + obj.mBranches(k).mComposition;
                end
                obj.mNumUnits = obj.mNumUnits + obj.mBranches(k).mNumUnits;
            end
            obj.mFormula = newFormula;
            obj.mInferredFormula = newInferredFormula;
            obj.mConciseFormula = newConciseFormula;
            if ~isempty(obj.mStem) && obj.mStem(1).mID == 1
                obj.updateUnitID();
            end
        end

        function nextID = updateUnitID( obj, baseID )
            if nargin < 2 || isempty( baseID )
                baseID = 1;
            end
            for k = 1 : length( obj.mStem )
                obj.mStem(k).mID = baseID;
                baseID = baseID + 1;
            end
            for k = 1 : length( obj.mBranches )
                baseID = obj.mBranches(k).updateUnitID( baseID );
            end
            nextID = baseID;
        end

        function result = copy( obj )
            if isempty( obj )
                result = [];
            else
                fields = fieldnames( obj(1) );
                num = length(obj);
                result = CGlycan.empty( 0, num );
                for k = 1 : num
                    for f = 1 : length( fields )
                        if ~isempty( obj(k).(fields{f}) )
                            if strcmp( fields{f}, 'stem' ) || strcmp( fields{f}, 'branches' )
                                result(k).(fields{f}) = obj(k).(fields{f}).copy;
                            elseif ~strcmp( fields{f}(1:2), 'cR' )
                                result(k).(fields{f}) = obj(k).(fields{f});
                            end
                        end
                    end
                end
            end
        end
        
        function mono = find_mono_by_id(obj, id)
            mono = [];
            if id < 0, return; end
            
            for k = 1 : length( obj.mStem )
                if obj.mStem(k).mID == id
                    mono = obj.mStem(k); return;
                end
            end
            
            for k = 1 : length( obj.mBranches )
                mono = obj.mBranches(k).find_mono_by_id( id );
                if ~isempty( mono ) return; end
            end
        end
        
        function startNodeID = parse( obj, formula, startNodeID, dist_reducing_end )
            if nargin < 3 || isempty(startNodeID), startNodeID = 1; end
            if nargin < 4 || isempty(dist_reducing_end), dist_reducing_end = 0; end
            
            if startNodeID == 1
                formula = strtrim( formula );
                formula = strrep( formula, 'NeuAc', 'Neu5Ac' );
                idxes = strfind( formula, ')' );
                for k = length(idxes) : -1 : 1
                    if idxes(k) < length(formula)
                       if ~strcmp( formula(idxes(k)+1), ']' )
                           formula = [formula(1:idxes(k)), ' ', formula(idxes(k)+1:end)];
                       end
                    end
                end
            end
            
            obj.mFormula = formula;
            obj.mStem = []; obj.mBranches = [];
            
            % find branching points
            branch_pos = []; num = 0;            
            for k = 1 : length( formula )
                if formula(k) == '['
                    if num == 0
                        branch_pos(1,end+1) = k;
                    end
                    num = num + 1;
                elseif formula(k) == ']'
                    num = num - 1;
                    if num == 0
                        branch_pos(2,end) = k;
                    end
                end
            end
            if num ~= 0, throw( MExeption( 'CGlycan.parse', 'Wrong formula: mis-matched [-].') ); end
            
            % pre-assign the stem, will be changed later.
            if ~isempty( branch_pos )
                stem_formula = strtrim( formula( branch_pos(end,end)+1 : end) );
            else
                stem_formula = formula;
            end
            
            % get branch formula
            branch_formula = cell( 1, size( branch_pos, 2 ) );
            if ~isempty( branch_pos )
                for k = 1 : size( branch_pos, 2 )
                    branch_formula{k} = strtrim( formula( branch_pos(1, k)+1 : branch_pos(2, k)-1 ) );
                end
            end
            
            % parse the stem
            stem_formula = strtrim( stem_formula );
            units = strsplit( stem_formula, ' ' );
            units = fliplr( units );
            obj.mStem = CMonosaccharide.empty(0, length(units));
            for k = 1 : length( units )
                ind = strfind( units{k}, '(' );
                if isempty( ind )
                    ind = strfind( units{k}, '-' );
                    linked_to_carbon = [];
                    if ~isempty( ind )
                        linked_to_carbon = str2double( units{k}(ind+1:end) );
                        units{k} = units{k}(1:ind-1);
                    end
                else
                    temp = units{k}(ind+1:end-1);
                    units{k} = units{k}(1:ind-1);
                    ind = strfind( temp, '-' );
                    linked_to_carbon = [];
                    if ~isempty( ind )
                        linked_to_carbon = str2double( temp(ind+1:end) );
                    end
                end
                
                obj.mStem(k) = CMonosaccharide( units{k}, obj.mPermethylated );
                obj.mStem(k).mID = startNodeID;
                obj.mStem(k).mDistance2Root = dist_reducing_end + k;
                obj.mComposition( obj.mStem(k).mClassID ) = obj.mComposition( obj.mStem(k).mClassID ) + 1;
                if ~isempty( linked_to_carbon )
                    if k > 1 % k == 1 is the root of a glycan or the roor of a branch
                        v = obj.mStem(k-1).mCarbon2Vertex( linked_to_carbon );
                        CMonosaccharide.link_monosaccharides( obj.mStem(k-1), v, obj.mStem(k), 1 );
                    else
                        obj.mStem(k).mLinkedTo.C = linked_to_carbon;
                        obj.mStem(k).mLinkedTo.type = 1;
                    end
                elseif k > 1
                    CMonosaccharide.link_monosaccharides( obj.mStem(k-1), 4, obj.mStem(k), -1 );
                end
                startNodeID = startNodeID + 1;
            end
            obj.mNumUnits = length( obj.mStem );

            if (obj.mStem(1).mID == 1) && strcmp( obj.mReducingEndModification, CMass.cReducingEndModification_Reduced )
                obj.mStem(1).set_reduced( 1 );
            end
            
            % parse the branches
            dist_2_leaf = 0;
            if isempty( branch_formula )
                obj.mBranches = CGlycan.empty(0,0);
            else
                newBranches = CGlycan.empty(0, length(branch_formula));
                for k = 1 : length( branch_formula ) 
                    % parse a branch
                    newBranches(k) = CGlycan;
                    newBranches(k).mPermethylated = obj.mPermethylated;
                    startNodeID = newBranches(k).parse( branch_formula{k}, startNodeID, obj.mStem(end).mDistance2Root );
                    if isempty( newBranches(k).mStem(1).mLinkedTo.V ) && ~isempty( newBranches(k).mStem(1).mLinkedTo.C )
                        newBranches(k).mStem(1).mLinkedTo.V = obj.mStem(end).mCarbon2Vertex( newBranches(k).mStem(1).mLinkedTo.C );
                    end
                    obj.mNumUnits = obj.mNumUnits + newBranches(k).mNumUnits;
                    if dist_2_leaf < newBranches(k).mStem(1).mDistance2Leaf
                        dist_2_leaf = newBranches(k).mStem(1).mDistance2Leaf;
                    end
                end
                obj.set_branches( newBranches ); % establish the links
            end
            for k = 1 : length( obj.mStem )
                obj.mStem(end-k+1).mDistance2Leaf = dist_2_leaf + k;
            end

            obj.update() % update formula and mass
        end
        
        function set_stem( obj, newStem )
        % Note that the newStem may not be a complete one.
            obj.mStem = newStem;
            if ~isempty( obj.mStem )
                nB = length( obj.mBranches );
                if nB
                    vertices = zeros(1,nB);
                    types = zeros(1,nB);
                    for k = 1 : nB
                        vertices(k) = obj.mBranches(k).mStem(1).mLinkedTo.V;
                        types(k) = obj.mBranches(k).mStem(1).mLinkedTo.type;
                    end
                    obj.mStem(end).set_linkage( vertices, types );
                    obj.mStem(end).update();
                end
            end
            obj.update(1);
        end
        
        function add_branches( obj, branches )
            if isempty( branches ), return; end
            if length(branches) + length(obj.mBranches) > length( obj.mStem(end).mLegalLinkedInVs )
                throw( MException('CGlycan:add_branches1', 'Tried to add too many branches') );
            end
            
            available = setdiff( obj.mStem(end).mLegalLinkedInVs, obj.mStem(end).mLinkedIn.V );
            for k = 1 : length( branches )
                if isempty( branches(k).mStem(1).mLinkedTo.V )
                    obj.mStem(end).linkage.inV = [obj.mStem(end).linkage.inV, available(1)];
                    branches(k).mStem(1).set_outLinkage( available(1), -1, obj.mStem(end).mID );
                    branches(k).update(1);
                    obj.mBranches = [obj.mBranches, branches(k)];
                    available = available(2:end);
                else
                    if branches(k).mStem(1).mLinkedTo.type > 0
                        if sum( available == branches(k).mStem(1).mLinkedTo.V ) == 0
                            throw( MException('CGlycan:add_branches2', 'Branch position occupied.') );
                        else
                            obj.mStem(end).add_linkedIn( branches(k).mStem(1).mLinkedTo.V, branches(k).mStem(1).mLinkedTo.type, branches(k).mStem(1).mID );
                            available = setdiff( available, branches(k).mStem(1).mLinkedTo.V );
                            obj.mBranches = [obj.mBranches, branches(k)];
                        end
                    else
                        if sum( available == branches(k).mStem(1).mLinkedTo.V ) > 0
                            obj.mStem(end).add_linkedIn( branches(k).mStem(1).mLinkedTo.V, branches(k).mStem(1).mLinkedTo.type, branches(k).mStem(1).mID );
                            available = setdiff( available, branches(k).mStem(1).mLinkedTo.V );
                        else
                            branches(k).mStem(1).mLinkedTo.V = available(1);
                            branches(k).mStem(1).mLinkedTo.type = -1;
                            branches(k).mStem(1).update();
                            obj.mStem(end).add_linkedIn( branches(k).mStem(1).mLinkedTo.V, branches(k).mStem(1).mLinkedTo.type, branches(k).mStem(1).mID );
                            available = available(2:end);
                        end
                        obj.mBranches = [obj.mBranches, branches(k)];
                    end
                end
                if obj.mStem(end).mPermethylated
                    obj.mStem(end).mVertices.modification{branches(k).mStem(1).mLinkedTo.V} = '-CH2';
                end
            end
            
            % sort branches
            obj.mStem(end).mLinkedIn.V = sort( obj.mStem(end).mLinkedIn.V );
            for k = 1 : length( obj.mBranches ) - 1
                for m = k+1 : length( obj.mBranches )
                    if obj.mBranches(k).mStem(1).mLinkedTo.V > obj.mBranches(m).mStem(1).mLinkedTo.V
                        temp = obj.mBranches(k); obj.mBranches(k) = obj.mBranches(m); obj.mBranches(m) = temp;
                    end
                end
            end
            
            obj.update_branching_info();
            obj.update(1);
        end
        
        function set_branches( obj, newBranches )
            bs = CGlycan.empty(0, length( newBranches )); % in case some of new branches are ions.
            for k = 1 : length( newBranches )
                if isa( class(newBranches(k)), 'CIon' )
                    bs(k) = CGlycan();
                    bs(k).mStem = newBranches(k).mStem;
                    bs(k).mBranches = newBranches(k).mBranches;
                    bs(k).update(1);
                else
                    bs(k) = newBranches(k);
                end
            end
            
            % order branches according to their linkedTo.V
            for k = 1 : length( bs ) - 1
                for m = k+1 : length( bs )
                    if isempty( bs(m).mStem(1).mLinkedTo.V )
                        continue;
                    end
                    if isempty( bs(k).mStem(1).mLinkedTo.V ) || ...
                       ( bs(k).mStem(1).mLinkedTo.V > bs(m).mStem(1).mLinkedTo.V )
                        temp = bs(k); bs(k) = bs(m); bs(m) = temp;
                    end
                end
            end
            obj.mBranches = bs;
            
            obj.mStem(end).clear_linkedIn( );
            flag = ones(1, length(bs));
            for b = 1 : length( bs )
                v = bs(b).mStem(1).mLinkedTo.V;
                if ~isempty( v )
                    CMonosaccharide.link_monosaccharides( obj.mStem(end), v, bs(b).mStem(1), 1 );
                    flag(b) = 0;
                end
            end
            
            % Assign default values if some branches do not have the linked_to information
            temp = bs(flag > 0);
            possible = intersect( obj.mStem(end).mLegalLinkedInVs, obj.mStem(end).mCleavage.vertices );
            possible = setdiff( possible, obj.mStem(end).mLinkedIn.V );
            for b = 1 : length( temp )
                CMonosaccharide.link_monosaccharides( obj.mStem(end), possible(b), temp(b).mStem(1), -1 );
            end
            
            obj.update();
        end
        
        function set_inferred_branches( obj, newBranches )
            num = length( newBranches );
            availableVs = intersect( obj.mStem(end).mLegalLinkedInVs, obj.mStem(end).mCleavage.vertices );
            if length( availableVs ) < num
                throw( MException('CGlycan:set_inferred_branches', 'Two many branches.') );
            end

            bs = CGlycan.empty( 0, num ); % in case some of new branches are Ions. copy them into Glycans
            for k = 1 : num
                if isa( newBranches(k), 'CIon' )
                    bs(k) = CGlycan();
                    bs(k).mStem = newBranches(k).mStem;
                    bs(k).mBranches = newBranches(k).mBranches;
                    bs(k).update(1);
                else
                    bs(k) = newBranches(k);
                end
            end
            
            for k = 1 : num-1
                for m = k+1 : num
                    if bs(k).mNumUnits > bs(m).mNumUnits
                        temp = bs(m); bs(m) = bs(k); bs(k) = temp;
                    elseif bs(k).mMass > bs(m).mMass
                        temp = bs(m); bs(m) = bs(k); bs(k) = temp;
                    end
                end
            end
            obj.mBranches = bs;
            
            for b = 1 : num
                CMonosaccharide.link_monosaccharides( obj.mStem(end), availableVs(b), bs(b).mStem(1), -1 );
            end

            for b = 1 : num
                bs(b).mStem(1).update();
                bs(b).update(1);
            end           
            obj.update(1);
        end

        function set_linkage( obj, linkageInfo )
        % linkageInfo(k) is the linkedTo.V of the k-th unit
            obj.mStem(1).clear_linkedIn();
            for k = 2 : length( obj.mStem )
                CMonosaccharide.link_monosaccharides( obj.mStem(k-1), linkageInfo(obj.mStem(k).mID), obj.mStem(k), 1 );
                obj.mStem(k).clear_linkedIn();
            end
            for b = 1 : length( obj.mBranches )
                obj.mBranches(b).set_linkage( linkageInfo );
                id = obj.mBranches(b).mStem(1).mID;
                CMonosaccharide.link_monosaccharides( obj.mStem(end), linkageInfo(id), obj.mBranches(b).mStem(1), 1 );
                obj.mBranches(b).update(1);
            end
            obj.update(1);
        end
        
        function remove_branches( obj, branches )
            changed = 0;
            for k = 1 : length( branches )
                for m = length( obj.mBranches ) : -1 : 1
                    if obj.mBranches(m) == branches(k)
                        changed = 1;
                        if m == 1
                            obj.mStem(end).linkage.inV = obj.mStem(end).linkage.inV(2:end);
                            obj.mStem(end).linkage.inV_type = obj.mStem(end).linkage.inV_type(2:end);
                        elseif m == length( obj.mBranches )
                            obj.mStem(end).linkage.inV = obj.mStem(end).linkage.inV(1:end-1);
                            obj.mStem(end).linkage.inV_type = obj.mStem(end).linkage.inV_type(1:end-1);
                        else
                            obj.mStem(end).linkage.inV = [obj.mStem(end).linkage.inV(1:m-1), obj.mStem(end).linkage.inV(m+1:end)];
                            obj.mStem(end).linkage.inV_type = [obj.mStem(end).linkage.inV_type(1:m-1), obj.mStem(end).linkage.inV_type(m+1:end)];
                        end
                        break;
                    end
                end
            end
            if changed
                obj.update(1);
            end
        end
        
        function update_branching_info( obj )
            bo = { 'r', 'b', 'a' }; 
            bo = bo(3-length(obj.mBranches)+1:end);
            for k = 1 : length( obj.mBranches )
                obj.mBranches(k).mBranchingCode = [obj.mBranchingCode, bo{k}];
                obj.mBranches(k).update_branching_info();
            end
        end
        
        function [mono, children] = get_local_structure( obj, id )
            mono = []; children = [];
            for k = 1 : length( obj.mStem ) - 1
                if obj.mStem(k).mID == id
                    mono = obj.mStem(k);
                    children = obj.mStem(k+1);
                    return
                end
            end
            
            if obj.mStem(end).mID == id
                mono = obj.mStem(end);
                for k = 1 : length(obj.mBranches)
                    children = [children, obj.mBranches(k).mStem(1)];
                end
                return;
            end
            
            for k = 1 : length(obj.mBranches)
                [mono, children] = obj.mBranches(k).get_local_structure(id);
                if ~isempty( mono )
                    return;
                end
            end
        end
        
        function disp_inferred_linkage( obj, header )
            if nargin < 2, header = '+'; end
            
            if isempty( obj.mBranches ) && length( obj.mStem ) <= 1
                return;
            end
            
            disp( [header, ' ', obj.mConciseFormula] );
            for k = 1 : length( obj.mStem )-1
                s = [header, ' ', num2str(obj.mStem(k).mID), ':', obj.mStem(k).class, ' <- ', num2str(obj.mStem(k+1).mID), ':', obj.mStem(k+1).class];
                disp(s);
                s = blanks( length(s) );
                for m = 1 : length( obj.mStem(k).inf_inLinkage )
                    temp = sprintf( '%d', obj.mStem(k).inf_inLinkage(m).C ); 
                    temp = [temp, sprintf( ' <- [%d:', length(obj.mStem(k).inf_inLinkage(m).peaks) ) ];
                    temp = [temp, ' ', sprintf( '%d,', obj.mStem(k).inf_inLinkage(m).peaks )];
                    temp(end) = ']';
                    disp( [s, temp] );
                end
            end
            
            if ~isempty(obj.mBranches)
                s = [header, ' ', num2str(obj.mStem(end).mID), ':', obj.mStem(end).class, ' <<=='];
                disp(s);
                s = blanks( length(s) );
                for m = 1 : length( obj.mStem(end).inf_inLinkage )
                    for v = 1 : length( obj.mStem(end).inf_inLinkage(m).V )
                        cs = obj.mStem(end).get_carbons( obj.mStem(end).inf_inLinkage(m).V(v) );
                        if v == 1
                            temp = num2str( obj.mStem(end).inf_inLinkage(m).C(v) );
                        else
                            temp = [temp, ';', num2str( obj.mStem(end).inf_inLinkage(m).C(v) )];
                        end
                    end
                    temp = [temp, sprintf( ' <- [%d:', length(obj.mStem(end).inf_inLinkage(m).peaks) )];
                    temp = [temp, ' ', sprintf( '%d,', obj.mStem(end).inf_inLinkage(m).peaks )];
                    temp(end) = ']';
                    disp( [s, temp] );
                end
            end

            for b = 1 : length(obj.mBranches)
                obj.mBranches(b).disp_inferred_linkage( [header, '+'] );
            end
        end
        
        function result = isNLinked( obj )
        % test if a glycan is an N-linked glycan
            result = 1;
            if obj.mNumUnits < 5
                result = 0;
            elseif ~strcmp( obj.mStem(1).class, 'HexNAc' )
                result = 0;
            else
                if length(obj.mStem) == 1
                    if length(obj.mBranches) > 2
                        result = 0;
                    elseif strcmp( obj.mBranches(1).mStem(1).class, 'HexNAc' ) && ...
                          (obj.mBranches(2).mNumUnits == 1) && strcmp( obj.mBranches(2).mStem(1).class, 'Fuc' ) && ...
                          isempty( obj.mBranches(1).mBranches ) && strcmp( obj.mBranches(1).mStem(2).class, 'Hex' )
                        return;
                    elseif strcmp( obj.mBranches(2).mStem(1).class, 'HexNAc' ) && ...
                          (obj.mBranches(1).mNumUnits == 1) && strcmp( obj.mBranches(1).mStem(1).class, 'Fuc' ) && ...
                          isempty( obj.mBranches(2).mBranches ) && strcmp( obj.mBranches(2).mStem(2).class, 'Hex' )
                        return;
                    else
                        result = 0;
                    end
                elseif strcmp( obj.mStem(2).class, 'HexNAc' );
                    if length(obj.mStem) == 2
                        if length(obj.mBranches) ~= 2
                            result = 0;
                        elseif strcmp( obj.mBranches(1).mStem(1).class, 'Hex' ) && ...
                                (obj.mBranches(2).mNumUnits == 1) && strcmp( obj.mBranches(2).mStem(1).class, 'Fuc' )
                            return;
                        elseif strcmp( obj.mBranches(2).mStem(1).class, 'Hex' ) && ...
                                (obj.mBranches(1).mNumUnits == 1) && strcmp( obj.mBranches(1).mStem(1).class, 'Fuc' )
                            return;
                        else
                            result = 0;
                        end
                    elseif length(obj.mStem) == 3
                        result = strcmp( obj.mStem(3).class, 'Hex' );
                    else
                        result = 0;
                    end
                else
                    result = 0;
                end
            end
        end
    end
end
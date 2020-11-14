classdef CIon < CGlycan  % By Pengyu Hong @ Brandeis University
    
    properties
        mType
        mSource = '';
        mMetal = '';
        mMetalNum = 1;
        mCleaveMono = CMonosaccharide.empty(0,0);
        mScore = 0;
        note = [];
        
        mReducingEndMonosaccharide = [];
        mNonReducingEndMonosaccharide = [];
    end
    
    methods
        function obj = CIon( t, ionMetal, glycan )
            if nargin < 1 || isempty( t )
                t = 'U';
            end
            obj.mType = t;
            
            if nargin < 2 || isempty( ionMetal )
                ionMetal = ''; % neutral mass
            end
            
            switch ionMetal
                case {'Proton', 'H'}
                    ionMetal = 'H';
                case {'Na', 'Sodium'}
                    ionMetal = 'Na';
                case {'Ce', 'Cesium'}
                    ionMetal = 'Ce';
            end
            
            obj.mMetal = ionMetal;
            
            if nargin == 3 && ~isempty( glycan )
                obj.mReducingEndModification = glycan.mReducingEndModification;
                obj.mPermethylated = glycan.mPermethylated;
                obj.mStem = glycan.mStem.copy();
                for k = 2 : length( obj.mStem )
                    obj.mStem(k-1).clear_linkedIn();
                    CMonosaccharide.link_monosaccharides( obj.mStem(k-1), obj.mStem(k).mLinkedTo.V, obj.mStem(k), obj.mStem(k).mLinkedTo.type );
                end
                obj.mBranches = glycan.mBranches.copy();
            else
                obj.mStem = CMonosaccharide.empty(0.0);
                obj.mBranches = CGlycan.empty(0.0);
            end
            obj.update(1);
        end
        
        function result = copy( obj )
            if isempty( obj )
                result = [];
            else
                fields = fieldnames( obj(1) );
                num = length(obj);
                result = CIon.empty( 0, num );
                for k = 1 : num
                    for f = 1 : length( fields )
                        if ~isempty( obj(k).(fields{f}) )
                            result(k).(fields{f}) = obj(k).(fields{f});
                            % if ( strcmp( fields{f}, 'stem' ) || strcmp( fields{f}, 'branches' ) )
                            %    result(k).(fields{f}) = obj(k).(fields{f}).copy;
                            % elseif ~strcmp( fields{f}, 'cleaveMono' ) && ~strcmp( fields{f}(1:2), 'cReducingEnd' )
                            %    result(k).(fields{f}) = obj(k).(fields{f});                            
                            % else
                            %    stop = 1;
                            % end
                        end
                    end
                    % if ~isempty( obj(k).mCleaveMono )
                    %    result(k).mCleaveMono = result(k).find_mono_by_id( obj(k).mCleaveMono.mID );
                    % end
                end
            end
        end
        
        function result = internal_cleave( obj )
            result = [];
            for k = 1 : length( obj.mStem )-1
                aIon = [];
                if obj.mType == 'B' || obj.mType == 'A' % cleave in the Y direction
                    aIon = CIon( [obj.mType, 'Y'] );
                    aIon.mMetal = obj.mMetal;
                    aIon.set_stem( obj.mStem(1:k).copy );
                elseif obj.mType == 'Y' || obj.mType == 'X' % cleave in the B direction
                    aIon = CIon( [obj.mType, 'B'] );
                    aIon.mMetal = obj.mMetal;
                    aIon.mStem = obj.mStem(k+1:end).copy();
                    if ~isempty( obj.mBranches )
                        aIon.set_branches( obj.mBranches.copy );
                    else
                        aIon.update();
                    end
                end
                if ~isempty( aIon )
                    result = [result, aIon];
                end
            end
            
            % Cleave g.mBranches
            for k = 1 : length( obj.mBranches )
                if obj.mBranches(k).cleaved == 0
                    [Ys, ~, ~, ~] = obj.mBranches(k).cleave();
                    for m = 1 : length( Ys )
                        aIon = CIon( [obj.mType, 'Y'] );
                        aIon.mMetal = obj.mMetal;
                        aIon.mStem = obj.mStem.copy();
                        b = CGlycan; b.mStem = Ys(m).mStem; b.mBranches = Ys(m).mBranches;
                        temp = [obj.mBranches(1:k-1), b, obj.mBranches(k+1:end)];
                        aIon.set_branches( temp.copy );
                        result = [result, aIon];
                    end
                end
            end
        end
        
        function update( obj, fast_mode )
            if nargin < 2
                fast_mode = 0;
            end
            update@CGlycan( obj, fast_mode );
            
            metalMass = 0;
            switch obj.mMetal
                case 'H'
                    metalMass = CMass.Proton;
                case 'Proton'
                    metalMass = CMass.Proton;
                case 'Ce'
                    metalMass = CMass.Cesium;
                case 'Li'
                    metalMass = CMass.Lithium;
                case 'Na'
                    metalMass = CMass.Sodium;
            end
            obj.mMass = obj.mMass + metalMass * obj.mMetalNum;
            
            if sum( obj.mType == 'B' ) > 0
                obj.mMass = obj.mMass - CMass.H2O;
            elseif sum( obj.mType == 'Z' ) > 0
                obj.mMass = obj.mMass - CMass.H2O;
            end
            
            %if ~isempty(obj.mReducingEndModification) && ~isempty( obj.mStem ) && obj.mStem(1).id == 1
            %    obj.mMass = obj.mMass + CGlycan.reducing_end_mass_compensation( obj.mReducingEndModification, obj.mPermethylated );
            %end
            
            obj.mFormula = [obj.mType, ': ', obj.mFormula];
            obj.mInferredFormula = [obj.mType, ': ', obj.mInferredFormula];
            if obj.mMetalNum == 1
                if strcmp(obj.mMetal, 'H') == 0 && strcmp( obj.mMetal, 'Proton' ) == 0
                    obj.mFormula = [obj.mFormula, ' +', obj.mMetal];
                end
            else
                obj.mFormula = [obj.mFormula, ' +', obj.mMetal, 'x', num2str(obj.mMetalNum)];
            end
        end
        
%         function disp( obj )
%             for k = 1 : length( obj )
%                 disp( [obj(k).mFormula, '. [', num2str(obj(k).mMass), ']'] );
%             end
%         end
    end

    methods(Static)
        function ions = unique( ions )
            len = length(ions);
            if len > 1
                flag = ones(1, length(ions));
                masses = [ions.mMass];
                umass = unique( masses );
                for k = 1 : length(umass)
                    idxes = find( masses == umass(k) );
                    if length(idxes) > 1
                        for ka = 1 : length(idxes)-1
                            if flag(idxes(ka)) == 0, continue; end
                            for kb = ka+1 : length(idxes)
                                if flag(idxes(kb))
                                    if strcmp( ions(idxes(ka)).mFormula, ions(idxes(kb)).mFormula )
                                        flag(idxes(kb)) = 0;
                                    end
                                end
                            end
                        end
                    end
                end
                ions = ions(flag>0);
            end
        end
    end
end
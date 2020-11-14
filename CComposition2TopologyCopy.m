classdef CComposition2TopologyCopy < handle % By Pengyu Hong @ Brandeis University
    properties
        mMaxNumberBranches;
    end
    
    properties % (Access = protected)
        mMap_Composition2TopologySuperSet;
        mMap_BisectCompositions;
        mMonosaccharides;
        mPromptIdx = 0;
    end
    
    methods
        function this = CComposition2TopologyCopy()
            this.mMaxNumberBranches = 2;
            
            this.mMap_Composition2TopologySuperSet = containers.Map({1}, {[]}); % Set key type = integer. Matlab default key type is string
            this.mMap_Composition2TopologySuperSet.remove(1);
            
            this.mMap_BisectCompositions = containers.Map({1}, {[]});
            this.mMap_BisectCompositions.remove(1);
            
            this.mMonosaccharides = CMonosaccharide.empty(0, CMonosaccharideSet.cNumberMonosaccharideClasses);
            for k = 1 : CMonosaccharideSet.cNumberMonosaccharideClasses
                this.mMonosaccharides(k) = CMonosaccharide( CMonosaccharideSet.cMonoClasses{k}, 1 );
            end
        end

        function topoSuperSet = create( this, composition ) 
        % This is the main entry of the recursive algorithm, which create a
        % sub-topology rooted at a monosaccharide with ID = rootMonoClassID.
        % -- composition: the composition of the topologies to be created.
        % -- preMonoClassID: what the created topologies will be attached to.
            if isempty(composition) || sum(composition) == 0
                return;
            end
            compositionKey = composition * [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000]';
            if this.mMap_Composition2TopologySuperSet.isKey( compositionKey )
               topoSuperSet = this.mMap_Composition2TopologySuperSet( compositionKey );
               %fprintf('\n');
               %topoSuperSet.printtopology();
               return
            end
            
            topoSuperSet = CTopologySuperSet;
            this.mPromptIdx = 0;
            
            % create rooted sub-topologies
            for root = 1 : CMonosaccharideSet.cNumberMonosaccharideClasses
                if composition(root) == 0, continue; end
                newComposition = composition;
                newComposition(root) = newComposition(root) - 1;
                
                newKey = newComposition * [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000]'; % Doesn't work for more than 9.
                if newKey == 0
                    newTS = CTopologySet;
                    newTS.mRootMono = this.mMonosaccharides(root);
                    newTS.mRootMonoClassID = root;
                    topoSuperSet.mTopologySets(end+1) = newTS;
                    %fprintf('\n');
                    %newTS.printtopology();
                    break;
                else
                    newTSS = this.create_rooted( newComposition, root );
                    if ~isempty( newTSS )
                        %fprintf('\n');
                        %newTSS.printtopology();
                        topoSuperSet.mTopologySets = [topoSuperSet.mTopologySets, newTSS.mTopologySets];
                    end
                end
            end
            
            if ~isempty( topoSuperSet.mTopologySets )                
                this.mMap_Composition2TopologySuperSet( compositionKey ) = topoSuperSet;
            end
        end
        
        function topoSuperSet = create_rooted( this, composition, rootMonoClassID ) 
        % Create a set of bi-branch sub-topologies rooted at rootMonoClassID.
        % -- composition: the composition of the topologies to be created.
        % -- preMonoClassID: what the created topologies will be attached to.
            if isempty(composition) || sum(composition) == 0
                topoSuperSet = [];
                return;
            end

            topoSuperSet = CTopologySuperSet;
            
            composition1 = zeros(1, CMonosaccharideSet.cNumberMonosaccharideClasses);
            for n1 = 0 : composition(1)
                composition1(1) = n1;
                for n2 = 0 : composition(2)
                    composition1(2) = n2;
                    for n3 = 0 : composition(3)
                        composition1(3) = n3;
                        for n4 = 0 : composition(4)
                            composition1(4) = n4;
                            for n5 = 0 : composition(5)
                                composition1(5) = n5;
                                for n6 = 0 : composition(6)
                                    composition1(6) = n6;
                                    for n7 = 0 : composition(7)
                                        composition1(7) = n7;
                                        for n8 = 0 : composition(8)
                                            composition1(8) = n8;
                                            if sum(composition1) == 0
                                                continue;
                                            end
                                            composition2 = composition - composition1;
                                            
                                            compositionKey1 = composition1 * [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000]'; % Doesn't work for more than 9.
                                            compositionKey2 = composition2 * [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000]'; % Doesn't work for more than 9.
                                            compositionRootKey1 = compositionKey1 + rootMonoClassID * 100000000;
                                            compositionRootKey2 = compositionKey2 + rootMonoClassID * 100000000;

                                            % avoid redundance as we are dealing with topology
                                            if ( compositionKey1 <= compositionKey2 )
                                                key1 = compositionKey1; 
                                                key2 = compositionKey2;
                                                rootKey1 = compositionRootKey1;
                                                rootKey2 = compositionRootKey2;
                                            else
                                                key1 = compositionKey2; 
                                                key2 = compositionKey1;
                                                rootKey1 = compositionRootKey2;
                                                rootKey2 = compositionRootKey1;
                                            end
                                            
                                            if this.mMap_BisectCompositions.isKey( rootKey1 )
                                                branch2key = this.mMap_BisectCompositions( rootKey1 );
                                                if any(branch2key == rootKey2)
                                                    continue;
                                                else
                                                    this.mMap_BisectCompositions( rootKey1 ) = [branch2key, rootKey2]; % Record it.
                                                end
                                            else
                                                this.mMap_BisectCompositions( rootKey1 ) = rootKey2; % Record it.
                                            end
                                            
                                            if this.mMap_Composition2TopologySuperSet.isKey( compositionKey1 )
                                                branchTSS1 = this.mMap_Composition2TopologySuperSet( compositionKey1 );                                                
                                            else
                                                branchTSS1 = this.create( composition1 ); % this.create() saves the topologies of composition1 
                                                                                          % in this.mMap_Composition2TopologySuperSet
                                            end
                                            if compositionKey2 ~= 0
                                                if this.mMap_Composition2TopologySuperSet.isKey( compositionKey2 )
                                                    branchTSS2 = this.mMap_Composition2TopologySuperSet( compositionKey2 );
                                                else
                                                    branchTSS2 = this.create( composition2 ); % this.create() saves topologies of composition2
                                                                                              % in this.mMap_Composition2TopologySuperSet
                                                end
                                            else
                                                branchTSS2 = [];
                                            end
                                            
                                            aTS = CTopologySet;
                                            aTS.mRootMono = this.mMonosaccharides(rootMonoClassID);
                                            aTS.mRootMonoClassID = rootMonoClassID;
                                            aTS.mSources{1} = branchTSS1;
                                            aTS.mSources{2} = branchTSS2;
                                            %fprintf('\n');
                                            %aTS.printtopology();
                                            
                                            topoSuperSet.mTopologySets(end+1) = aTS;
                                            % show the process is still runing.
                                            this.mPromptIdx = this.mPromptIdx + 1;
                                            if mod(this.mPromptIdx, 1000) == 0
                                                fprintf('.');
                                                if mod(this.mPromptIdx, 60000) == 0
                                                    fprintf(' %d\n', this.mPromptIdx );
                                                end
                                            end
                                        end % for n8
                                    end % for n7
                                end % for n6
                            end % for n5
                        end % for n4
                    end % for n3
                end % for n2
            end % for n1
        end % function create_rooted()
    end
end

classdef GWAWrapper < handle %By Hui Sun
    %GWAWrapper: Wrapper for XML-based GWA file generation
    %Every time you want to output some formulas to GWB, make a new instance of that class and add peak annotations to that object by providing peak intensity, mz_ratio and glycan formula, mass, type, score, then use writeGWA to write it to GWA file. 
    properties (Access = private)
        gwaNode;
        root;
        glycan
        pac;
    end
    
    methods
        function obj = GWAWrapper()
            obj.gwaNode = com.mathworks.xml.XMLUtils.createDocument('AnnotatedPeakList');
            obj.root = obj.gwaNode.getDocumentElement;
            annos = obj.gwaNode.createElement('Annotations');
            obj.root.appendChild(annos);
            obj.glycan = obj.gwaNode.createElement('Glycan');
            obj.glycan.setAttribute('structure','');
            annos.appendChild(obj.glycan);
            obj.pac = obj.gwaNode.createElement('PeakAnnotationCollection');
            annos.appendChild(obj.pac);
        end
        
        function add_comment( obj, comment )
            obj.root.appendChild(obj.gwaNode.createComment(comment));
        end
        
        function set_structure ( obj, formula )
            obj.glycan.setAttribute('structure', GWAWrapper.generate_gws(formula));
        end
        
        function add_peakanno(obj, intensity, mz_ratio, formula, mass, type, score)
            % First put PA, P, A, then put fragment into A
            currPa = obj.gwaNode.createElement('PeakAnnotation');
            obj.pac.appendChild(currPa);
            
            currP = obj.gwaNode.createElement('Peak');
            currPa.appendChild(currP);
            currP.setAttribute('mz_ratio',num2str(mz_ratio));
            currP.setAttribute('intensity',num2str(intensity));
            
            currA = obj.gwaNode.createElement('Annotation');
            currPa.appendChild(currA);
            
            currFrag = obj.gwaNode.createElement('FragmentEntry');
            currA.appendChild(currFrag);
            
            currFrag.setAttribute('fragment', GWAWrapper.generate_gws(formula));
            currFrag.setAttribute('mass',num2str(mass));
            currFrag.setAttribute('mz_ratio', num2str(mz_ratio));
            currFrag.setAttribute('name', type);
            currFrag.setAttribute('score', num2str(score));
        end
        
        function writeGWA(obj,filename)
            xmlwrite(filename, obj.gwaNode);
        end
    end
    
    methods(Static, Access=private)
        function preGWS = convert( formula )
            % the CGlycan parser can only parse one layer of branch so recursive call is necessary
            % the parser just get the stem instead of the root hence i need to extract the root one by one
            [stem, branches] = CGlycan.break_formula(formula);
            roots = strsplit(stem);
            root = [];
            for i = 1:length(roots)
                root = ['--4b1D-', roots{i}, ',p', root];
            end
            branch = [];
            if ~isempty(branches)
                for idx = 1:length(branches)
                    if idx == 1
                        branch = GWAWrapper.convert(branches{idx});
                    else
                        branch = ['(',branch,')',GWAWrapper.convert(branches{idx})];
                    end
                end
            end
            preGWS = [root, branch];
        end
        
        function GWS = generate_gws(formula)
            if isempty(formula)
                GWS='';
            else
                preGWS = GWAWrapper.convert(formula);
                GWS = ['freeEnd--?', preGWS(4:end)];
            end
        end
    end
end


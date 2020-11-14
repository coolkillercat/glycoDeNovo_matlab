a=CComposition2TopologyCopy2(); %Don't remove any redundency
% b=CComposition2Topology(); %Record composition1 & composition2
c=CComposition2TopologyCopy(); %Record b&root
% d=CComposition2TopologyCopy(); %Record b&root
%composition=[3,1,1,1,1,1,0,0]
% 
tic;test_raw=a.create(composition);toc
tic;test_root=c.create(composition);toc
% %test_root2=d.create(composition);
% 
tic;test_raw.reconstruct_formulas;toc
tic;test_root.reconstruct_formulas;toc
% 
%tic;test_root2.reconstruct_formulas2;toc
% 
isequal(test_raw.mFormulas,test_root.mFormulas)
% test_root.reconstruct_formulas();
% %isequal(test_raw.mFormulas,test_root2.mFormulas)
%d=CComposition2TopologyCopy();
%test_root2=d.create([1,1,1,1,0,0,0,0])
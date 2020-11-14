tic
c2t = CComposition2Topology;
tss = c2t.create( [2 3 4 4 0 2 0 0] );
tss.reconstruct_formulas;
toc
%%
tss2 = c2t.create( [0 0 3 4 0 0 0 0] );
tss2.reconstruct_formulas;

%%
g = CGlycan(1); % permethylated glycans
% g = CGlycan; for native glycans.

g.parse( tss.mFormulas{2} );

[~, Bs, ~, ~, ~, Cs] = g.cleave( 'Proton' );
spectrumMass = [ [Bs.mMass], [Cs.mMass] ];

spectrumMass = floor( spectrumMass * 100000 ) / 100000;
spectrumMass = unique( spectrumMass );

%%

% for allpossible compositions
%   produce all possible topologies
%   for each topology
%       produce its spectrumMass
%       save the topology and its spectrumMass


% Second way
% cm2c = CMass2Composition;
% cm2c.build_mass_2_composition( mapfile ) % create all possible % compositions that satisfy CMass2Composition.mMonoMaxNums
% uniqueMass = unique( floor(cm2c.mMasses * 100000) / 100000 );
% for k = 1 : length(uniqueMass)
%     currentMass = uniqueMass(k);
%     idxes = find( abs(cm2c.mMasses - currentMass) < currentMass * 0.000005 )
%     compositions = cm2c.mCompositions(idxes,:);
% end

%%
% Given a massrange [massLow, massHigh]
% retrieve all topologies, whose mass is in [massLow, massHigh], and their
% spectrumMass



par.dispT = 1;

par.write_con   = 'on';
par.plotall     = 'off';
par.savegraph   = 'off';
par.savegraph   = 'on';

% STATIC LOOP
par.minI = 2;    % necessary for the constraints to have time to adjust
par.maxI = 10000;
par.tolI = 1e-6; % make sure that par.tolI is not higher than par.tolK

% OUTER STATIC LOOP: to update the derivatives using past iterations
par.minK = 2;
par.maxK = 100;
par.tolK = 1e-4;

par.damp1 = 0.1;
par.damp1_ = 0.0;
par.dampm1 = 0.0;

par.damp2 = 0.1;
par.damp2_ = 0.0;
par.dampm2 = 0.0;

% HJB LOOP: to update the derivatives using past iterations
par.minT1 = 10;
par.minT2 = 10;

par.maxT = 1000;
par.tolT = 0.01;

% with implicit method, dt can be large; however, it cannot be too big since we use guesses in the static
% loop; if there are convergence issues, reduce dt and try again; dt2 should be small
par.dt  = 1; 
par.dt2 = 1;

% STENCIL
par.maxD = 1000;
par.tolD = 1e-4;

% maxP refers to how far we go in the direction in stencil; refer to the algorithm in paper for details
par.maxP  = 2;
par.nextr = par.maxP; % taille de l'extrapolation

% STATE VARIABLES GRID
vec1a = 0.01;
vec1b = 0.99;
vec2a = 0.01;
vec2b = 0.99;

par.n1 = 30;
par.n2 = 30; 

vec1 = linspace(vec1a,vec1b,par.n1); % linear spacing is generally more stable
vec2 = linspace(vec2a,vec2b,par.n2);



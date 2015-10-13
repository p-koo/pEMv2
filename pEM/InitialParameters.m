function [D0 P0 S0] = InitialParameters(numStates)
%-------------------------------------------------------------------------- 
% Summary: This function returns the initial diffusion coefficient, static
% localization noise, and population fraction values for a given number of
% states.  The user should adjust these numbers according to particular 
% system under investigation for better convergence.
% 
% Input:
%       numStates = number of diffusive states to analyze
%
% Output:
%       D0 = vector of diffusion coefficients
%       P0 = vector of population fractions 
%       S0 = vector of static localization noises
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

% 1 state model
D{1} = .3;
P{1} = 1;
S{1} = .05;

% 2 state model
D{2} = [.1 .5];
P{2} = [.5 .5];
S{2} = [.05 .05];

% 3 state model
D{3} = [.01 .3 .6];
P{3} = [.3 .4 .3];
S{3} = [.05 .05 .05];

% 4 state model
D{4} = [.01 .2 .4 .6];
P{4} = [.25 .25 .25 .25];
S{4} = [.05 .05 .05 .05];

% 5 state model
D{5} = [.01 .1 .2 .4 .6];
P{5} = [.2 .2 .2 .2 .2];
S{5} = [.05 .05 .05 .05 .05];

% 6 state model
D{6} = [.01 .05 .15 .3 .5 .8];
P{6} = [.2 .1 .1 .2 .2 .2];
S{6} = [.05 .05 .05 .05 .05 .05];

% 7 state model
D{7} = [.001 .05 .15 .3 .5 .7 1];
P{7} = [.2 .1 .1 .2 .2 .1 .1];
S{7} = [.05 .05 .05 .05 .05 .05 .05];

% 8 state model
D{8} = [.0001 .03 .1 .2 .4 .6 .8 1];
P{8} = [.1 .1 .2 .1 .2 .1 .1 .1];
S{8} = [.05 .05 .05 .05 .05 .05 .05 .05];

% 9 state model
D{9} = [.0001 .03 .1 .2 .4 .6 .8 1 1.2];
P{9} = [.1 .1 .1 .1 .1 .1 .2 .1 .1 ];
S{9} = [.05 .05 .05 .05 .05 .05 .05 .05 .05];

% 10 state model
D{10} = [.0001 .03 .1 .2 .35 .5  .75 1 1.2 2];
P{10} = [.1 .1 .1 .1 .1 .1 .1 .1 .1 .1];
S{10} = [.05 .05 .05 .05 .05 .05 .05 .05 .05 .05];


D0 = D{numStates};
P0 = P{numStates};
S0 = S{numStates};

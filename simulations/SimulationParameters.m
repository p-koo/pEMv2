function  simParams = SimulationParameters(set,uniformA)

% select diffusive state data set
switch set 
    case 1
		Dindex = [.05 .15 .25 .4];
		Sindex = [.04 .04 .04 .04];
		Lindex = [.13 0 0 0];
		LnoiseIndex = [.01 0 0 0 ];
		Aindex = [1  1  .9 .6];
		Vindex = [0  0  0 0];
        Pindex = ones(1,4)/4;
    case 2
        Dindex = [.005 .1  .3];
        Sindex = [.01 .01  .01];
        Lindex = [.05  .2  0];
        LnoiseIndex = [.01 .01  0 ];
        Aindex = [1 1 1  ];
        Vindex = [0 0 0  ];
        Pindex = ones(1,3)/3;
    case 3
        Dindex = [.001  .03   .2     .45];
        Sindex = [.04  .04    .04    .04];
        Lindex = [0 0 0 0];
        LnoiseIndex = [0 0 0 0];
        Aindex = [1 .7 1 .9];
        Vindex = [0 0 0 0 ];
        Pindex = ones(1,4)/4;
    case 4
        Dindex = [0.0001 0.01 0.07 0.17];
        Sindex = [.05  .05    .05    .05];
        Lindex = [0 0 0 0];
        LnoiseIndex = [0 0 0 0];
        Aindex = [1 1 1 1];
        Vindex = [0 0 0 0 ];
        Pindex = [0.3 0.1 0.2 0.4];
    
    case 5
        Dindex = [0.0001 0.05 0.17 0.17];
        Sindex = [.03  .03    .03    .03];
        Lindex = [0 0.15  0  .2];
        LnoiseIndex = [0 .01 0 0.01];
        Aindex = [1 1 1 1];
        Vindex = [0 0 0 0 ];
        Pindex = [0.25 0.25 0.25 0.25];    
    case 6
        Dindex = [.05 .15];
        Sindex = [.03  .03];
        Lindex = [.1  0];
        Dindex = [.05];
        Lindex = [.15];
        LnoiseIndex = [.01 0 ];
        Aindex = [1 1];
        Vindex = [0 0];
        Pindex = [.6 .4];
        Pindex = 1;
    case 7
        Dindex = [.05];
        Sindex = [.03];
        Lindex = [0 0 0];
        LnoiseIndex = [0 0 0];
        Aindex = [.6];
        Vindex = [0 0 0 0 ];
        Pindex = [1];
    case 8
        Dindex = [.05 .15 .1];
        Sindex = [.03  .03 .03];
        Lindex = [.1  0 0];
        LnoiseIndex = [.01 0 0];
        Aindex = [1 1 .6];
        Vindex = [0 0 0];
        Pindex = [.4 .3 .3];
    case 9
        Dindex = [.05 .1];
        Sindex = [.03 .03];
        Lindex = [.1  0];
        LnoiseIndex = [.01 0];
        Aindex = [1 1];
        Vindex = [0 0];
        Pindex = [.7 .3];
                
    otherwise 
        disp('no cases found');
end
simParams.numStates = length(Dindex);

% get diffusion mode
dmode = ones(1,simParams.numStates);
index = find(Lindex ~= 0);
dmode(index) = 2;
index = find(Vindex ~= 0);
dmode(index) = 3;
index = find(Aindex ~= 1);
dmode(index) = 4;

% setup state structure 
state = struct;
for i = 1:simParams.numStates
    state(i).mode = dmode(i);
    state(i).D = Dindex(i);
    state(i).locNoise = Sindex(i);
    state(i).L = Lindex(i);
    state(i).Lnoise = LnoiseIndex(i);
    state(i).v = Vindex(i);
    state(i).A = Aindex(i);
end
simParams.state = state;
simParams.Pindex = Pindex;

% generate transition matrix
if uniformA ~= 0
    A = InitializeTransitionMatrix(simParams.numStates,uniformA);
else
    switch simParams.numStates
        case 2
            A = [.992 .008; .1 .9];
        case 3
            A = [.9950 .0010 .004; .001 .9950 .004; .015 .015 .97];
        case 4
            A = [1 .004 .004 .001; .005 .97 .001 .005; .008 .001 .98 .008; .0001 .005 .005 1];
        case 5
            A = [1 .004 .004 .001 .0001; .005 .97 .001 .005  .005; .008 .001 .98 .008 .005; .0001 .005 .005 1 .0001; .0001 .005 .005 .0001 .98];
        otherwise
    end
end
A = A./repmat(sum(A,2),1,size(A,2));
simParams.A = A;



function  A = TransitionMatrix(simParams,dt,uniformA)

% generate transition matrix
if uniformA == 1
    A = InitializeTransitionMatrix(simParams.numStates,.997);
else
    switch simParams.numStates
        case 2
            k = [0.1460   87.7875; 
                41.8652    1.9156];            
        case 3
            k = [ 0.1629   87.9504  125.7585; 
                96.1503    0.1846   96.1503; 
                125.7585   87.9504    0.1629];
        case 4
            k = [0.1629  100.5531  100.5531  125.7585;
               95.9843    0.2050  125.2468   95.9843;
               87.7329  125.5409    0.3127   87.7329;
              167.6434   96.5158   96.5158    0.1827];
        case 5
            k = [  0.1647  100.5549  100.5549  125.7603  167.6254;
               96.0767    0.2975  125.3392   96.0767   96.0767;
               87.8238  125.6319    0.4036   87.8238   96.3694;
              167.6452   96.5176   96.5176    0.1845  167.6452;
              167.2817   96.1540   96.1540  167.2817    0.1883];
        otherwise
    end
    A = exp(k*dt);
end


 
% Assume a target positioned at x = 1, travelling with speed v = 0.1
state = [ 1000+3.8676; -9.2; 1500-11.7457; -9.2; (2*pi/180)/8 ];

% Create an instance of a 1D Constant Velocity model with observation noise
% variance equal to 0.1 and specifying a mapping between the 1st index of
% the measurement vector to the 1st index of the state vector and from the
% 2nd index of the measurement vector to the 3rd index of the state vector.
ct = CoordinatedTurnModelX('VelocityErrVariance',5,'TurnRateErrVariance',(pi/180));

% View the transition matrix and process covariance matrices
Q = ct.covariance();

% Predict the target's position and velocity after the interval has passed
newState  = ct.feval(state);

% Do the same as above, but this time add process noise to the prediction
newState2 = ct.feval(state,true);

% Generate 50 random noise samples from the dynamic model
noise = ct.random(50);

% Check how likely the predictions we made are
lik = ct.pdf(newState,state);
lik2 = ct.pdf(newState2,state); % HINT: newState2 should be less likely

trajectory = state;
for(i=2:10000)
    trajectory(:,end+1) = ct.feval(trajectory(:,end));
    plot(trajectory(1,:),trajectory(3,:))
    pause(0.1)
end

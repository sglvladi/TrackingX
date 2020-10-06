% Define an intial state vector
state = [1;0.1];

% Create an instance of a 2D Random Walk model, with i.i.d. noise
rw = RandomWalkX('NumDims',2,'ErrVariance',0.1);

% View the transition matrix and process covariance matrices
F = rw.feval();
Q = rw.covar();

% Predict the target's position and velocity after the interval has passed
newState  = rw.feval(state);

% Do the same as above, but this time add process noise to the prediction
newState2 = rw.feval(state,true);

% Generate 50 random noise samples from the dynamic model
noise = rw.random(50);

% Check how likely the predictions we made are
lik = rw.pdf(newState,state);
lik2 = rw.pdf(newState2,state); % HINT: newState2 should be less likely
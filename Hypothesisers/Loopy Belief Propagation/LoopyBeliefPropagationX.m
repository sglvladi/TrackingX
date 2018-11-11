classdef LoopyBeliefPropagationX < HypothesiserX
% LoopyBeliefPropagationX Abstract class
%
% Summary of LoopyBeliefPropagationX:
% This is an implementation of the Loopy Belief Propagation
% algorithm described in [1], based on a modified version of the code
% provided on [2]
%
% LoopyBeliefPropagationX Properties:
%   None
%
% LoopyBeliefPropagationX Methods:
%   + LoopyBeliefPropagationX - Constructor method
%
% (+) denotes puplic properties/methods
%
% [1] J. L. Williams, "Marginal multi-bernoulli filters: RFS derivation of MHT, JIPDA, 
%     and association-based member," in IEEE Transactions on Aerospace and Electronic Systems, 
%     vol. 51, no. 3, pp. 1664-1687, July 2015.
% [2] J. L. Williams and R. A. Lau, "Convergence of loopy belief propagation for data association," 
%     2010 Sixth International Conference on Intelligent Sensors, Sensor Networks and Information Processing, 
%     Brisbane, QLD, 2010, pp. 175-180.
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        NetObj
        LikelihoodMatrix
        AssocWeightsMatrix
        Mode = 'Track-Oriented'
        ConvergeThreshold = 10^(-6);
        MaxIterations = 1000;
    end
    
    methods
        function this = LoopyBeliefPropagationX(varargin)
        % LOOPYBELIEFPROPAGATIONX Constructor method
        % 
        % Parameters
        % ----------
        % ConvergeThreshold: scalar
        %   Convergence threshold. Lower value should yield more accurate
        %   approximations of the association prbabilities. 
        %   ( default: 10^-6, generally << 1, as per [1])
        % Mode: string
        %   Specifies the mode of operation of the EHM algorithm, which can
        %   be set to either:
        %       1) 'Track-Oriented'
        %       2) 'Measurement-Oriented'
        %   (default='Track-Oriented')
        % 
        % Usage
        % -----
        % * LoopyBeliefPropagationX('ConvergenceThreshold',c) creates a class instance of the
        %   algorithm, configured with ConvergenceThreshold c.
        % * LoopyBeliefPropagationX(__, 'Mode',mode) creates a class instance of the
        %   algorithm, configured to operate in Mode mode. Mode can be
        %   either 'Track-Oriented' (default) or 'Measurement-Oriented'.
        %
        % See also LoopyBeliefPropagationX/hypothesise
            % Return early if no arguments are supplied
            if(nargin==0)
                return;
            end
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'ConvergeThreshold'))
                        this.ConvergeThreshold  = config.ConvergeThreshold;
                    end
                    if (isfield(config,'Mode'))
                        this.Mode  = config.Mode;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            if (isfield(config,'ConvergeThreshold'))
                this.ConvergeThreshold  = config.ConvergeThreshold;
            end
            if (isfield(config,'Mode'))
                this.Mode  = config.Mode;
            end
        end
        
        function AssocWeightsMatrix = hypothesise2(this,LikelihoodMatrix)
        % HYPOTHESISE Generate hypotheses and calculate association weights
        % 
        % Parameters
        % ----------
        % LikelihoodMatrix: (NumMeasurements+1 x NumTracks) matrix
        %   A matrix whose (i,j)th element is the likelihood of the (i-1)th
        %   measurement, given the measurement prediction of the jth track.
        %   The row (1,:) represents the dummy (null) measurement likelihood
        %   for each track.
        %
        % Returns
        % -------
        % AssocWeightsMatrix: (NumMeasurements+1 x NumTracks) matrix
        %   A matrix whose (i,j)th element is the joint probability that the 
        %   (i-1)th measurement is associated to the jth track. The row (1,:) 
        %   represents the dummy (null) measurement association probability
        %   for each target.
        %
        % Usage
        % -----
        % * hypothesise(ehm) approximates all valid association hypotheses and  
        %   calculates the respective association weights, based on the
        %   likelihood matrix ehm.LikelihoodMatrix.
        % * hypothesise(ehm,lik) approximates all valid association hypotheses and  
        %   calculates the respective association weights, based on the
        %   likelihood matrix lik. ehm.LikelihoodMatrix is also updated
        %   with lik.
        %
        % See also LoopyBeliefPropagationX/hypothesise
            % First check to see if a structure was received
            if(nargin==2)
                this.LikelihoodMatrix = LikelihoodMatrix;
            end
            
            
            [numTargets,numMeasurements] = size(this.LikelihoodMatrix); 
            numMeasurements = numMeasurements-1;
            
            % Perform LBP as described in Appendix of [1]
            muba = ones(numTargets,numMeasurements);
            muba_tilde = zeros(numTargets,numMeasurements);
            om = ones(1,numMeasurements);
            on = ones(1,numTargets);
            iter = 0;
            while max(max(abs(muba-muba_tilde)))>this.ConvergeThreshold && iter<this.MaxIterations       
%                 disp(max(max(abs(muba-muba_tilde))));
                muba_tilde = muba;

                prodfact = muba .* this.LikelihoodMatrix(:,2:end);
                sumprod= this.LikelihoodMatrix(:,1)+ sum(prodfact,2);
                divider = sumprod(:,om) - prodfact;
                muab = this.LikelihoodMatrix(:,2:end) ./ divider;

                summuab= 0.01 + sum(muab,1);
                divider =  summuab(on,:) - muab;
                muba = 1 ./ divider;
                
                iter = iter + 1;
            end 
            AssocWeightsMatrix = zeros(numTargets,numMeasurements+1);
            if(iter>=this.MaxIterations)
                AssocWeightsMatrix(:,1) = 1;
            else
                for i = 1:numTargets
                    AssocWeightsMatrix(i,1) = this.LikelihoodMatrix(i,1);
                    AssocWeightsMatrix(i,2:end) = this.LikelihoodMatrix(i,2:end).*muba(i,:);
                end
            end
            AssocWeightsMatrix = AssocWeightsMatrix./sum(AssocWeightsMatrix,2);
        end
        
        function AssocWeightsMatrix = hypothesise(this,LikelihoodMatrix)
        % HYPOTHESISE Generate hypotheses and calculate association weights
        % 
        % Parameters
        % ----------
        % LikelihoodMatrix: (NumMeasurements+1 x NumTracks) matrix
        %   A matrix whose (i,j)th element is the likelihood of the (i-1)th
        %   measurement, given the measurement prediction of the jth track.
        %   The row (1,:) represents the dummy (null) measurement likelihood
        %   for each track.
        %
        % Returns
        % -------
        % AssocWeightsMatrix: (NumMeasurements+1 x NumTracks) matrix
        %   A matrix whose (i,j)th element is the joint probability that the 
        %   (i-1)th measurement is associated to the jth track. The row (1,:) 
        %   represents the dummy (null) measurement association probability
        %   for each target.
        %
        % Usage
        % -----
        % * hypothesise(ehm) approximates all valid association hypotheses and  
        %   calculates the respective association weights, based on the
        %   likelihood matrix ehm.LikelihoodMatrix.
        % * hypothesise(ehm,lik) approximates all valid association hypotheses and  
        %   calculates the respective association weights, based on the
        %   likelihood matrix lik. ehm.LikelihoodMatrix is also updated
        %   with lik.
        %
        % See also LoopyBeliefPropagationX/hypothesise
            % First check to see if a structure was received
            if(nargin==2)
                this.LikelihoodMatrix = LikelihoodMatrix;
            end
            this.LikelihoodMatrix = this.LikelihoodMatrix;
            [numTargetsp1,numMeasp1] = size(this.LikelihoodMatrix);
            numMeas = numMeasp1-1;
            numTargets = numTargetsp1-1;
            wnew = 0.01*ones(1,numMeas);
            % this.LikelihoodMatrix dimensions n x m+1
            % wnew dimensions m x 1
            % AssocWeightsMatrix, newTargets dimensions same

            eps_conv_threshold = 1e-4;

            mu = ones(numTargets,numMeas); % mu_ba
            mu_old = zeros(numTargets,numMeas);
            nu = zeros(numTargets,numMeas); % mu_ab

            AssocWeightsMatrix = zeros(numTargetsp1,numMeasp1);

            % Run LBP iteration
            numIter = 0;
            if(numTargets>0)
                while (max(abs(mu(:)-mu_old(:))) > eps_conv_threshold) ...
                        && numIter<this.MaxIterations
                    numIter = numIter+1;

                    mu_old = mu;

                    for i = 1:numTargets
                        prd = this.LikelihoodMatrix(i+1,2:end).*mu(i,:);
                        s = this.LikelihoodMatrix(i+1,1) + sum(prd);
                        nu(i,:) = this.LikelihoodMatrix(i+1,2:end) ./ (s - prd);
                    end

                    for j = 1:numMeas
                        s = this.LikelihoodMatrix(1,j+1) + sum(nu(:,j));
                        mu(:,j) = 1./(s - nu(:,j));
                    end
                end
            end

            % Calculate outputs--for existing tracks then for new tracks
            for i = 1:numTargets
                s = this.LikelihoodMatrix(i+1,1) + sum(this.LikelihoodMatrix(i+1,2:end).*mu(i,:));
                AssocWeightsMatrix(i+1,1) = this.LikelihoodMatrix(i+1,1)/s;
                AssocWeightsMatrix(i+1,2:end) = this.LikelihoodMatrix(i+1,2:end).*mu(i,:)/s;
            end

            for j = 1:numMeas
                s = this.LikelihoodMatrix(1,j+1) + sum(nu(:,j));
                AssocWeightsMatrix(1,j+1) = this.LikelihoodMatrix(1,j+1)/s;
            end
        end
    end        
 end

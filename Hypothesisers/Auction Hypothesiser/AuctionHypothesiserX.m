classdef AuctionHypothesiserX < HypothesiserX
% AuctionHypothesiserX Abstract class
%
% Summary of AuctionHypothesiserX:
% This is an implementation of the Efficient Hypothesis Management (EHM)
% algorithm described in [1]
%
% HypothesiserX Properties:
%   None
%
% HypothesiserX Methods:
%   + HypothesiserX - Constructor method
%
% (+) denotes puplic properties/methods
%
% [1] Maskell, Simon et al. “Fast Mutual Exclusion.” (2004).
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        NetObj
        LikelihoodMatrix
        AssocWeightsMatrix
        Mode = 'Track-Oriented'
        timeout = 5
    end
    
    methods
        function this = AuctionHypothesiserX(varargin)
        % AuctionHypothesiserX Constructor method
        % 
        % Parameters
        % ----------
        % Mode: string
        %   Specifies the mode of operation of the EHM algorithm, which can
        %   be set to either:
        %       1) 'Track-Oriented'
        %       2) 'Measurement-Oriented'
        %   (default='Track-Oriented')
        % 
        % Usage
        % -----
        % * AuctionHypothesiserX() creates a EHM object handle
        % * AuctionHypothesiserX('Mode',mode) creates an EHM
        %   object handle, configured to operate in Mode mode. Mode can be
        %   either 'Track-Oriented' (default) or 'Measurement-Oriented'.
        %
        % See also AuctionHypothesiserX/hypothesise
            % Return early if no arguments are supplied
            if(nargin==0)
                return;
            end
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
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
            config = parser.Results;
            if (isfield(config,'Mode'))
                this.Mode  = config.Mode;
            end
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
        % * hypothesise(obj) generates all valid association hypotheses and  
        %   calculates the respective association weights, based on the
        %   likelihood matrix ehm.LikelihoodMatrix.
        % * hypothesise(obj,lik) generates all valid association hypotheses and  
        %   calculates the respective association weights, based on the
        %   likelihood matrix lik. obj.LikelihoodMatrix is also updated
        %   with lik.
        %
        % See also AuctionHypothesiserX/hypothesise
            % First check to see if a structure was received
            if(nargin==2)
                this.LikelihoodMatrix = LikelihoodMatrix;
            end
            
            [NumTracks, NumMeasurements] = size(this.LikelihoodMatrix(:,2:end)); 
            
            assigns = auctionx(log(this.LikelihoodMatrix));
            AssocWeightsMatrix = zeros(NumTracks,NumMeasurements+1);
            a = (assigns-1)*(NumTracks)+[1:NumTracks]';
            
%             costMatrix = -log(this.LikelihoodMatrix);
%             costOfNoAssignment = costMatrix(1,1);
%             costMatrix = costMatrix(:,2:end);
%             [assignments, unassignedTracks, ~] = ...
%                 assignauction(costMatrix,costOfNoAssignment);
%             AssocWeightsMatrix = zeros(NumTracks,NumMeasurements+1);
%             AssocWeightsMatrix(unassignedTracks,1) = 1;
%             a = (assignments(:,2))*(NumTracks)+assignments(:,1);
            
            AssocWeightsMatrix(a) = 1;
            this.AssocWeightsMatrix = AssocWeightsMatrix;
        end
    end
 end

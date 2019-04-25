classdef (Abstract) DataAssociatorX < BaseX 
% DataAssociatorX Abstract class
%
% Summary of DataAssociatorX:
% This is the base class for all TrackingX data associators.
% Any custom defined data associator should be derived from this DataAssociatorX base class. 
%
% DataAssociatorX Properties:
%   None
%
% DataAssociatorX Methods:
%   + DataAssociatorX - Constructor method
%   + associate - Parforms data association
%   ~ performGating - Performs gating
%   ~ computeLikelihoods - Performs measurement likelihoood computations
%   ~ performClustering - Performs target/measurement clustering
%   ~ evaluateAssociations - Computes the target/measurement association weights
%
% (+) denotes puplic properties/methods
% (~) denotes abstract properties/methods
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        TrackList
        MeasurementList
        Gater
        Hypothesiser
        Clusterer
    end
    
    methods (Abstract, Access = protected)
        performGating(this);        % Performs gating
        computeLikelihoods(this);   % Performs measurement likelihoood computations
        performClustering(this);    % Performs target/measurement clustering
        evaluateAssociations(this); % Computes the target/measurement association weights
    end
    
    methods
        function this = DataAssociatorX(varargin)
        % DATAASSOCIATORX Constructor method
        %            
        end
        
        function associate(this)
        % ASSOCIATE Parforms data association
        %
        % A generic template for the operation of a data associator is
        % drawn. Subclasses can optionally be defined by simply
        % implementing the class abstract methods. If this is not desired,
        % the method can also be overriden.
        
            % Perform Gating
            this.performGating();
            
            % Perform Clustering
            this.performClustering();
            
            % Compute Likelihoods
            this.computeLikelihoods();

            % Evaluate the association probabilities
            this.evaluateAssociations();
        end
    end
end
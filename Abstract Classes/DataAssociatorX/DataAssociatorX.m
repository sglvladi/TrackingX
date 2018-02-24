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
%
% (+) denotes puplic properties/methods
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        TrackList
        MeasurementList
        Gater = []
        Hypothesiser = []
        Clusterer = [] 
    end
    
    methods (Abstract, Access = protected)
        performGating(this);
        computeLikelihoods(this);
        performClustering(this);
        evaluateAssociations(this);
    end
    
    methods
        function this = DataAssociatorX(varargin)
        % DATAASSOCIATORX Constructor method
        %   
        % DESCRIPTION: 
        % * DataAssociatorX() returns a "DataAssociatorX" object handle
            
        end
        
        function associate(this)
        % ASSOCIATE Parforms data association
        %
        % A generic template for the operation of a data associator is
        % drawn. By doing so, subclasses can be defined by simply
        % implementing the class abstract methods. If this is not desired,
        % the method can simply be overriden.
        
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
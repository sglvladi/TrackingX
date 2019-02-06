classdef BasicPlotterX < BaseX 
% BasicPlotterX class
%
% Summary of BasicPlotterX:
% This is a class implementation of a basic plotter object to provide basic
% configuration and plotting functionalities.
%
% BasicPlotterX Properties:
%   + Model       - A handle to the underlying state-space model
%   ~ gfxHandles_ - List of handles to all active axes/plots created 
%                   by and/or registered with the plotter.
%
% BasicPlotterX Methods:
%   + BasicPlotterX - Constructor method
%   + axis - Provides basic manipulation of available axes.
%   + figure - Provides basic manipulation of open figures.
%   + plot - Plot generic data to screen
%   + plotObject - Plot a given TrackingX object
%
% (+) denotes puplic properties/methods
% (~) denotes protected properties/methods
% 
% See also MeasurementSimulatorX, MetricGeneratorX
    
    properties %(Access = protected)
        gfxHandles_ = gobjects(3,1);
    end
    
    properties
        Model = StateSpaceModelX.empty();
    end
    
    methods
        
        function this = BasicPlotterX(varargin)
        % BasicPlotterX Constructor
        %
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   A handle to the underlying state-space model. This is required
        %   for the plotter to perform correct transformations between
        %   state and measurement space.
            
            this.Model = varargin{1};
        end
        
        function handles = plot(this, varargin)
        % PLOT Plot generic (2D) data using the standard PLOT command
        %
        % Ultimately, the function collects and forwards all inputs to an
        % internally invoked <a href="matlab:web('https://uk.mathworks.com/help/matlab/ref/plot.html')">PLOT</a> command. 
        %
        % Parameters
        % ----------
        % plotSettings: various
        %   A set of plot settings, as per PLOT's documentation
        %
        % Returns
        % -------
        % handles: matlab.graphics.Graphics 
        %   A list of MATLAB Graphics handles to the graphics components 
        %   used or generated during the call. They are given in the 
        %   following order:
        %       - handles(1): gcf
        %       - handles(2): gca
        %       - handles(3): plot
        %
        % See also PLOT, BasicPlotterX/figure 
            
            plotSettings = varargin{:};
            if isa(varargin{1}, 'Axes')
                ax = varargin{1};
                plotSettings = varargin{2:end};
            elseif isa(this.gfxHandles_(1),'Figure')
                fig = this.gfxHandles_(1);
                ax = fig.CurrentAxes;
            else
                ax = gcf.CurrentAxes;  
            end
            
            handles    = gobjects(3,1);
            handles(3) = plot(ax, plotSettings);
            handles(2) = ax;
            handles(1) = gcf;
            
            this.gfxHandles_ = handles;
        end
        
        function varargout = figure(this, varargin)
        % FIGURE Minor modification/extension to MATLAB's own <a href="matlab:web('https://uk.mathworks.com/help/matlab/ref/plot.html')">FIGURE</a> command.
        %
        % DESCRIPTION:
        %
        %   T = obj.FIGURE(cmd,[h,]....) will execute the user specified 
        %       figure command, where cmd is a character vector that can   
        %       take the following values:
        %       * {'n' or 'new'}: Generate a new window and return its
        %         handle T. Essentially equivalent to calling obj.figure(....)
        %       * {'c' or 'close'}: Attempt to close the figure with 
        %         handle h, if one is provided, or the (last) active 
        %         figure and return a locial T to indicate success/failure. 
        %       * 'clf': Clear the figure with handle h, if one
        %         is provided, or the (last) active figure and return the
        %         handle T of the cleared figure.
        %
        %  In all other aspects the commad should behave as per MATLAB's
        %  FIGURE command documentation.
        
            switch(nargin)
                case(1)
                    this.gfxHandles_(1) = figure;
                    varargout{1} = this.gfxHandles_(1);
                otherwise
                    command = varargin{1};
                    switch(command)
                        case {'n', 'new'}
                            this.gfxHandles_(1) = figure(varargin{2:end});
                            varargout{1} = this.gfxHandles_(1);
                        case {'c', 'close'}
                            varargout{1} = close(varargin{2:end});
                        case {'clf'}
                            this.gfxHandles_(1) = clf(varargin{2:end});
                            varargout{1} = this.gfxHandles_(1);
                        otherwise
                            this.gfxHandles_(1) = figure(varargin{:});  
                            varargout{1} = this.gfxHandles_(1);                          
                    end
            end
        end
        
        function handles = plotObject(this, jobject, varargin)
        % plotObject Plot a given TrackingX object.
        %
        % This function forwards all inputs to an internal call to a PLOT 
        % command, after which it additonally stores graphics handles.
        %
        % Parameters
        % ----------
        % object: any plottable TrackingX component
        %   The object to be visualised
        % plotSettings: various
        %   A set of plot settings, as per PLOT's documentation.
        %
        % Returns
        % -------
        % handles: matlab.graphics.Graphics 
        %   A list of MATLAB Graphics handles to the graphics components 
        %   used or generated during the call. They are given in the 
        %   following order:
        %       - handles(1): gcf
        %       - handles(2): gca
        %       - handles(3): plot
        %
        % See also BasicPlotterX/plot, BasicPlotterX/figure
            
            if nargin>2
                plotSettings = varargin{:};
                if isa(varargin{1}, 'Axes')
                    ax = varargin{1};
                    plotSettings = varargin{2:end};
                elseif isa(this.gfxHandles_(1),'Figure')
                    fig = this.gfxHandles_(1);
                    ax = fig.CurrentAxes;
                else
                    ax = gcf.CurrentAxes;  
                end
            else
                plotSettings = [];
                ax = gca;
            end
            
            handles    = gobjects(3,1);
            
            switch (class(jobject))
                case 'GroundTruthTrackX'
                    for track = jobject
                        path = [track.Trajectory.Vector];
                        handles(3) = plot(ax, path(1,:), path(3,:), '--');
                        hold on;
                    end
                case 'MeasurementListX'
                    for list=jobject
                        vectors = [list.Measurements.Vector];
                        s_vectors = this.Model.Measurement.finv(vectors);
                        handles(3) = plot(ax, s_vectors(1,:), s_vectors(3,:), 'r+');
                        hold on;
                    end
            end
            handles(2) = ax;
            handles(1) = gcf;
            
            this.gfxHandles_ = handles;
        end
                
    end
    
    methods (Access = protected)
        function success = addAxisHandle_(this, h)
            success = false;
            if(~sum(ismember(ax,this.gcaHandles_)))
                this.gcaHandles_(end+1) = gca;
            end
        end
        function success = addFigureHandle_(this, h)
            success = false;
            if(~sum(ismember(ax,this.gcaHandles_)))
                this.gcaHandles_(end+1) = gca;
            end
        end
    end
end
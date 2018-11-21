classdef OSPAX < MetricGeneratorX
% OSPAX Class 
%
% Summary of OSPAX:
% This is a class implementation of the Optimal Sub Pattern Assignment (OSPA)
% metric generator.
%
% OSPAX Properties:
%   + CutOffThreshold
%   + Order
%
% MultinomialResamplerX Methods:
%    OSPAX  - Constructor method
%    evaluate - Perform multinomial resampling 
% 
% See also SystematicResamplerX
    properties
        CutOffThreshold
        Order
    end
    
    methods
        function this = OSPAX(varargin)
        % OSPAX Constructor method
        %   
        % DESCRIPTION: 
        % * OSPAX(C,P) returns a "OSPAX" object handle configured with the
        %   cut-off threshold C and order paramater P.
        %
        % See also OSPAX/evaluate
            
        end
        
        function [metric,locDist,cardDist] = evaluate(this,TrackList,GroundTrouth)
        % EVALUATE Evaluate the OSPA metric of a set of tracks compared to
        % the groundtrouth.
        %   
        % DESCRIPTION: 
        % * [METRIC,LOCDIST,CARDDIST] = EVALUATE(THIS, TRACKLIST, GROUNDTRUTH) 
        %   evaluates the OSPA metric METRIC, including its localisation LOCDIST
        %   and cardinality CARDDIST components, of the TRACKLIST compared to the 
        %   GROUNDTRUTH.
        %
        % See also SystematicResamplerX/resample
            
            [metric,locDist,cardDist] = ospa_dist(GroundTrouth,TrackList,this.CutOffThreshold,this.Order);
        end
    end

    methods (Static)
    
        function [dist, varargout]= ospa_dist(X,Y,c,p)
            % This is the MATLAB code for OSPA distance proposed in
            % 
            % D. Schuhmacher, B.-T. Vo, and B.-N. Vo, "A consistent metric for performance evaluation in multi-object filtering," IEEE Trans. Signal Processing, Vol. 56, No. 8 Part 1, pp. 3447– 3457, 2008.
            % http://ba-ngu.vo-au.com/vo/SVV08_OSPA.pdf
            % ---BibTeX entry
            % @ARTICLE{OSPA,
            % author={D. Schuhmacher and B.-T. Vo and B.-N. Vo},
            % journal={IEEE Transactions on Signal Processing},
            % title={A Consistent Metric for Performance Evaluation of Multi-Object Filters},
            % year={2008},
            % month={Aug},
            % volume={56},
            % number={8},
            % pages={3447-3457}}  
            %---

            % Compute OSPA distance between two finite sets X and Y
            % Inputs: X,Y-   matrices of column vectors
            %        c  -   cut-off parameter
            %        p  -   p-parameter for the metric
            % Output: scalar distance between X and Y
            % Note: the Euclidean 2-norm is used as the "base" distance on the region

            if nargout ~=1 & nargout ~=3
               error('Incorrect number of outputs'); 
            end
            if isempty(X) & isempty(Y)
                dist = 0;
                if nargout == 3
                    varargout(1)= {0};
                    varargout(2)= {0};
                end
                return;
            end
            if isempty(X) | isempty(Y)
                dist = c;
                if nargout == 3
                    varargout(1)= {0};
                    varargout(2)= {c};
                end
                return;
            end


            %Calculate sizes of the input point patterns
            n = size(X,2);
            m = size(Y,2);

            %Calculate cost/weight matrix for pairings - fast method with vectorization
            XX= repmat(X,[1 m]);
            YY= reshape(repmat(Y,[n 1]),[size(Y,1) n*m]);
            D = reshape(sqrt(sum((XX-YY).^2)),[n m]);
            D = min(c,D).^p;

            %Compute optimal assignment and cost using the Hungarian algorithm
            [assignment,cost]= Hungarian(D);

            %Calculate final distance
            dist= ( 1/max(m,n)*( c^p*abs(m-n)+ cost ) ) ^(1/p);

            %Output components if called for in varargout
            if nargout == 3
                varargout(1)= {(1/max(m,n)*cost)^(1/p)};
                varargout(2)= {(1/max(m,n)*c^p*abs(m-n))^(1/p)};
            end

        end
    end
end


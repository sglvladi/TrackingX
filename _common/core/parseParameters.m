function varargout = parseParameters(Defaults,callerVarargin,varargin)
%   Parses parameter name-value pairs quickly. Does not handle options.
%
%   Within yourFunction(...), specify defaults in a struct: Defaults.x = 1;
%   Defaults.y = 1; Then call [x,y] = parseParameters(Defaults,varargin);
%   The user can then call yourFunction(...,'Y',2,'x',3); to change x and
%   y. See Examples section for more. This function combines brevity in
%   usage, high performance, and the convenience and clarity of parameter
%   name-value pairs. MATLAB's inputParser class is relatively slow.
%   Simpler methods are verbose and/or only handle positional arguments.
%   Many FEX parameter parsing solutions use assignin, which is convenient,
%   but also slow and unsafe.
%
%   Arguments
%   ---------
%   Defaults : struct
%       The default values for your parameters, specified as
%       Defaults.(paramName) = defaultValue. Order of declaration matters
%       when getting parsed output.
%   callerVarargin : cell array
%       Contains parameter name-value pairs. In normal usage, just pass
%       your function's varargin here.
%   isCaseSensitive : bool, optional
%       Default false. This optional positional argument allows you to 
%       specify case sensitivity. Useful if you have parameter names like
%       's' and 'S'. MATLAB parameters are usually not case sensitive.
%   expandStruct : bool, optional
%       Default true. This optional positional argument by default expands
%       the output struct, enabling syntax such as [x,y,z] =
%       parseParameters(...). If you would rather use the output struct
%       itself (i.e. Results = parseParameters(...) ), set expandStruct to
%       true.
%
%   Returns
%   -------
%   varargout : cell
%       The values of the parameters, returned in the same order as the
%       Defaults struct fields were initialized. If expandStruct is false,
%       returns a struct formatted like the Defaults struct.
% 
%   See Also
%   --------
%   inputParser
%
%   Notes
%   -----
%   In my tests this function is 9-20x faster than inputParser. To keep
%   your programs responsive, use inputParser for functions called less
%   than 1000 times and parseParameters for functions called less than
%   10,000 times. I recommend positional argument parsing with nargin
%   instead of name-value pairs for functions called less than 100,000
%   times. Any more and consider not using optional arguments.
%   
%   Positional argument parsing methods include using exist, cell arrays,
%   and nargin. parseParameters is nearly on par with the exist method,
%   2-3x slower than the cell method, and 7-10x slower than the nargin
%   method. 
%
%   I sacrificed features for performance. Users may be interested in other
%   solutions with more features. I tested the performance of several FEX
%   parameter name-value pair parsing functions. The fastest was
%   loadOptions by N. Brahms. The syntax is verbose, but it is only a tad
%   (1.2-1.5x) slower than parseParameters and allows for type checking.
%   Other fairly fast FEX solutions are getargs, parseargs, and
%   parse_pv_pairs.
%
%   Here's an extra tip. To make MATLAB's inputParser more concise, use::
%
%       >> results = struct2cell(p.Results);
%       >> [param1, param2, param3] = results{:};
%
%   Note: This documentation adheres to NumPy-style (see link below). It is
%   compatible with Sphinx using the matlabdomain and napoleon extensions.
%   https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
%
%   Examples
%   --------
%   .. highlight:: matlab
%
%   Basic usage within a function::
%
%       function foo(varargin)
%           Defaults.a = 'a';
%           Defaults.b = 'b';
%           [a,b] = parseParameters(Defaults,varargin);
%           disp(['a: ' a ', b: ' b '.']);
%       end
%       >> foo('a','A')
%       a: A, b: b.
%
%   Order of Defaults initialization matters. Case insensitive by default::
%
%       >> Defaults.x = 3.14159;
%       >> Defaults.msg = 'hello';
%       >> [x,msg] = parseParameters(Defaults,{'MsG','bye','x',1.61803})
%       x = 1.61803
%       msg = 'bye'
%
%   Using case sensitivity and returning a structure::
%
%       function foo(varargin)
%           Defaults.FOOL = 1;
%           Defaults.food = 'apple';
%           Fooey = parseParameters(Defaults,varargin,true,false);
%           disp(Fooey);
%       end
%       >> foo('food','orange','fool',0)
%       FOOL: 1
%       food: 'orange'
%
%   Notice that Fooey.FOOL was not changed. The options support using empty
%   arrays for defaults.
%
%   Copyright 2015 Jeffrey Chiou and everyone else. Feel free to copy,
%   distribute, and modify.

defaultParamNames = fieldnames(Defaults);
numInputs = length(callerVarargin);

%% Check for optional inputs and assign defaults.
% Using the nargin method for speed, sacrificing some flexibility
isCaseSensitive = false;
expandStruct = true;
if nargin > 2
    if ~isempty(varargin{1}), isCaseSensitive = varargin{1}; end;
end
if nargin > 3
    if ~isempty(varargin{2}), expandStruct = varargin{2}; end;
end

%% Check for an even number of arguments
if mod(numInputs,2)==1
    throw( MException('PARSEPARAMETERS:missingValue',...
        'A parameter seems to be missing a value.') );
end

%% Assign values. Repeated code for speed.
% Repeated code is bad for maintanability and readibility, but this
% function is optimized for speed. It runs faster if the conditionals are
% outside the loop.

if isCaseSensitive
    %% Case sensitive block.  
    for i = 1:2:numInputs
        % Make sure parameter names are strings
        try
            % Counterintuitively, using find() here makes overall
            % performance faster. ~isempty() is faster than any(), and
            % access into defaultParamNames with a number is faster than
            % with a logical array.
            match = find( strcmp(callerVarargin{i},defaultParamNames), 1 );
        catch ME
            causeException = MException('PARSEPARAMETERS:nameNotString',...
                'Parameter names must be strings.');
            ME = addCause(ME,causeException);
            rethrow(ME);
        end
        
        % Overwrite defaults
        if ~isempty(match)
            Defaults.(defaultParamNames{match}) = callerVarargin{i+1};
        end        
    end
else
    %% Case insensitive block. Nearly identical repeat of above.
    % Differences are find(...) instead of find(...,1), strcmpi() instead
    % of strcmp(), length()==1 instead of ~isempty(), and an error check.
    for i = 1:2:numInputs
        try
            match = find( strcmpi(callerVarargin{i},defaultParamNames) );
        catch ME
            causeException = MException('PARSEPARAMETERS:nameNotString',...
                'Parameter names must be strings.');
            ME = addCause(ME,causeException);
            rethrow(ME);
        end
        
        % Check for double default parameter names (ex. Defaults.x and
        % Defaults.X). This is not a problem with the case sensitive
        % version. In the case sensitive version, ~isempty() is just barely
        % faster than length()==1.
        if length(match)==1
            Defaults.(defaultParamNames{match}) = callerVarargin{i+1};
        else
            throw( MException('PARSEPARAMETERS:badDefaultNames',...
                ['Two default parameter names are identical, besides '...
                'case. Set case sensitivity to true in parseParameters']));
        end        
    end
end % of isCaseSensitive if-else

%% Choose output format
% Speed of both options are about the same.
if expandStruct
    varargout = struct2cell(Defaults);
else
    varargout = {Defaults};
end

end % of function
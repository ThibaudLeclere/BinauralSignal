classdef BinauralSignal
    % Defining class properties
    properties (SetAccess = private)
        leftChannel double
        rightChannel double% Binaural signal values (left/right as columns)
        label  char ; % Signal label
        fs  double; % Sampling frequency
    end
    properties (Dependent)
       isDiotic logical 
       signalLength double 
    end
    
    % -- Dynamic Methods ------------------------------
    methods
        % Defining Constructor -----------------
        function obj = BinauralSignal(audio, fs, label) % Signal, Sampling Frequency, Label
            
            if nargin < 1
                audio = [];
            end
            
            if nargin < 2
                fs = [];
            end
            if nargin < 3
                label = '';
            end
            
            
            % Construciton from signal vectors
            if isnumeric(audio)
                if isempty(audio)
                    left = [];
                    right = [];
                    
                elseif size(audio,2) == 2
                    left = audio(:,1);
                    right = audio(:,2);
                elseif size(audio,2) == 1
                    left = audio(:,1);
                    right = left;
                else
                    % error
                end
            end
            
            % Construction from audio filename
            if ischar(audio)
                [~, label,~] = fileparts(audio);
                [signals, fs] = audioread(audio);
                
                if size(signals,2) == 2
                    left = signals(:,1);
                    right = signals(:,2);
                elseif size(signals,2) == 1
                    left = signals(:,1);
                    right = left;
                else
                    % error
                end
                
            end
            
            obj.fs = fs;
            obj.label = label;
            obj.leftChannel = left;
            obj.rightChannel = right;
            
            
        end
    
        % Overload operators ------------------
        function obj = mtimes(x,binauralObj) % Multiplication using '*' operator
            N = numel(binauralObj);
            obj(N) = BinauralSignal;
            
            
            if isscalar(x) && isnumeric(x)
                for n = 1:N
                    obj(n).leftChannel = x*binauralObj(n).leftChannel;
                    obj(n).rightChannel = x*binauralObj(n).rightChannel;
                    obj(n).fs = binauralObj(n).fs;
                    obj(n).label = binauralObj(n).label;
                end
            end
        end
        
        function obj = times(obj1,obj2) % multiplication using '.*' operator
            obj = BinauralSignal;
            if isnumeric(obj2) && isscalar(obj2)
                obj.leftChannel = obj2*obj1.leftChannel;
                obj.rightChannel = obj2*obj1.rightChannel;
                obj.fs = obj1.fs;
            elseif isa(obj2,'BinauralSignal') && isa(obj1,'BinauralSignal')
                if obj1.fs == obj2.fs
                    obj.fs = obj1.fs;
                    if length(obj1.leftChannel) == length(obj2.leftChannel) &&...
                       length(obj1.rightChannel) == length(obj2.rightChannel)
                        
                        
                        obj.leftChannel = obj2.leftChannel.*obj1.leftChannel;
                        obj.rightChannel = obj2.rightChannel.*obj1.rightChannel;
                    end
                end
            end
        end
        
        function obj = plus(obj1,obj2)
            N1 = numel(obj1);
            N2 = numel(obj2);
            
            if N1 == N2
                obj(N1) = BinauralSignal;
                for n = 1:N1
                    if isa(obj1(n),'BinauralSignal') && isa(obj2(n),'BinauralSignal')
                        lengthDifference = obj1(n).signalLength - obj2(n).signalLength;
                        if obj1(n).fs == obj2(n).fs
                            if lengthDifference > 0 % If obj1 is longer than obj2
                                % Add zeros at the end of obj2
                                obj2(n).leftChannel = vertcat(obj2(n).leftChannel, zeros(lengthDifference));
                                obj2(n).rightChannel = vertcat(obj2(n).rightChannel, zeros(lengthDifference));
                                
                            elseif lengthDifference < 0 % If obj2 is longer than obj1
                                % Add zeros at the end of obj1
                                obj1.leftChannel = vertcat(obj1(n).leftChannel, zeros(-lengthDifference,1));
                                obj1.rightChannel = vertcat(obj1(n).rightChannel, zeros(-lengthDifference,1));
                                
                            end
                            
                            
                            obj(n).leftChannel = obj1(n).leftChannel + obj2(n).leftChannel;
                            obj(n).rightChannel = obj1(n).rightChannel + obj2(n).rightChannel;
                            
                            
                            obj(n).fs = obj1(n).fs;
                        else
                            error('Binaural signals have different sample rates')
                        end
                    end
                end
                obj = reshape(obj, size(obj1));
            else
                error('Binaural signals have different sizes')
            end
        end
        
       
        
        % Defining convolution
        function obj = conv(obj,x)
            if isa(x,'BinauralSignal') % Convolve two BinauralSignal objects between them
                if obj.fs ~= x.fs
                    error('Sampling rates do not match')
                else
                    obj.leftChannel = conv(obj.leftChannel,x.leftChannel);
                    obj.rightChannel = conv(obj.rightChannel,x.rightChannel);
                end
            elseif isnumeric(x) && isvector(x) % Convolve one BinauralSignal with a single signal (diotic)
                obj.leftChannel = conv(obj.leftChannel,x);
                obj.rightChannel = conv(obj.rightChannel,x);
            else
                error('Wrong type of data for x')
            end
        end
        
        %-------------------------------------
        %--------------------------------------
        function diotic = get.isDiotic(obj)
            diotic = isequal(obj.leftChannel, obj.rightChannel);
        end
        
        function L = get.signalLength(obj)
            Nobj = numel(obj);
            L = zeros(Nobj,1);
            for n = 1:Nobj
                L(n) = length(obj.leftChannel);
            end
            L = reshape(L,size(obj));
        end
        % Other methods ------------------------        
        % Filtering Binaural Channels with the same filter
        function obj = filter(obj,filterObject)
            if isa(filterObject,'digitalFilter')
                obj.leftChannel = filter(filterObject,obj.leftChannel);
                obj.rightChannel = filter(filterObject,obj.rightChannel);
            end             
        end
        
        % Compute Interaural Correlation Function
        function [ICF, lags] = computeICF(obj,varargin)
            % computeICF Interaural Coherence Function between left and
                % right channels of a BinauralSignal object by computing
                % the crosscorrelation operator.
                
                % [ICF, lags] = computeICF(BinauralSignal,varargin) returns
                % intercorrelation function and the associated lag (in
                % samples).
                
            if isa(obj,'BinauralSignal')                
                if ~isempty(obj.leftChannel) && ~isempty(obj.rightChannel)
                    if ~isempty(varargin)
                        maxlags = round(varargin{1}*obj.fs);
                        [ICF, lags] = xcorr(obj.leftChannel,obj.rightChannel,maxlags,'coeff');
                    else
                        [ICF, lags] = xcorr(obj.leftChannel,obj.rightChannel,'coeff');
                    end
                else
                    ICF = NaN;
                    lags = NaN;
                end
            else
                error('Wrong entry type. Obj must be a BinauralSignal object')
            end
        end
        
        % Compute Mean ILD
        function ILD = computeMeanILD(obj)
            ILD = zeros(numel(obj),1);
            for n = 1:numel(obj)                
                ILD(n) = 10*log10(sum(obj(n).leftChannel.^2) ./ sum(obj(n).rightChannel.^2));
            end
            ILD = reshape(ILD,size(obj));
        end
        
        % Compute Interaural Coherence
        function ICFanalysis  = computeIC(obj,varargin)
            % computeIC performs an intercorrelation operation on the left
                % and right channels of a BinauralSignal object.
                
                % The function returns a structure of the same size as the
                % input object with the following fields: 'ICF' (Interaural
                % Coherence Function), 'IC' (the maximum of ICF), 'Delay'
                % (the location of IC in seconds), 'Label' (unchanged from
                % the original object).
                
                
            if isempty(varargin)
                maxlag = [];
            else
                maxlag = varargin{1}; 
            end
            
            Nelements = numel(obj);
            ICFanalysis(Nelements) = struct('InterCorrelationFunction',[],'InterauralCoherence',[],'Delay',[],'Label','');
            h = waitbar(0,'','Name','Computing Interaural Coherence');
            for n = 1:numel(obj)
                waitbar(n/Nelements,h, sprintf('%s (%d/%d)',strrep(obj(n).label,'_','\_'),n,Nelements))
                [ICF, lag] = computeICF(obj(n),maxlag);
                [IC, idx] = max(abs(ICF));
                delay = lag(idx)/obj(n).fs;
                
                ICFanalysis(n).InterCorrelationFunction = ICF;
                ICFanalysis(n).InterauralCoherence = IC;
                ICFanalysis(n).Delay = delay;
                ICFanalysis(n).Label = obj(n).label;
            end
            delete(h)
            ICFanalysis = reshape(ICFanalysis,size(obj));
        end
        
        % Gammatone filtering
        function filteredSignal = gammatoneFiltering(obj,analyzer)
            Nbands = length(analyzer.center_frequencies_hz);
            filteredSignal(numel(obj),Nbands) = BinauralSignal();
            
            for n = 1:numel(obj)
                filteredSignal_left = gfb_analyzer_process(analyzer,obj(n).leftChannel);
                filteredSignal_right =  gfb_analyzer_process(analyzer,obj(n).rightChannel);
                for k = 1:Nbands
                    filteredSignal(n,k).leftChannel = real(filteredSignal_left(k,:))';
                    filteredSignal(n,k).rightChannel = real(filteredSignal_right(k,:))';
                    filteredSignal(n,k).fs = analyzer.fs;
                    filteredSignal(n,k).label = [obj(n).label ' (' num2str(analyzer.center_frequencies_hz(k)) ' Hz)'];
                end
            end
            filteredSignal = squeeze(reshape(filteredSignal,[size(obj) Nbands]));
            if iscolumn(filteredSignal)
                filteredSignal = filteredSignal';
            end
        end
        
        % Envelope extraction
        function envelopes = envelopeExtraction(obj)
            envelopes(numel(obj),1) = BinauralSignal();
            for n = 1:numel(envelopes)
               
                leftEnv = abs(hilbert(real(obj(n).leftChannel)));
                rightEnv = abs(hilbert(real(obj(n).rightChannel)));
                
                envelopes(n).leftChannel = leftEnv;
                envelopes(n).rightChannel = rightEnv;
                envelopes(n).fs = obj(n).fs;
                envelopes(n).label = [obj(n).label '_env'];
            end
            envelopes = reshape(envelopes,size(obj));
        end
        % Compute Interaural Level Difference (ILD)
        function ILD_distribution = computeShortTermILD(obj,analyzer)            
            
            obj_filtered_Env = obj.gammatoneFiltering(analyzer).envelopeExtraction;  
            
            S = size(obj_filtered_Env);
            Nelements = numel(obj_filtered_Env);
            ILD_distribution(numel(obj)) = struct('Distribution',[],'Label','','CenterFrequencies',analyzer.center_frequencies_hz,'SampleRate',48000);
            for n = 1:Nelements
                [quotient, remainder] = euclideanDivision(n-1,prod(S(1:end-1)));
                frequencyBand = quotient +1;
                position = remainder + 1;
                
                ILDvector = 10*log10(obj_filtered_Env(n).leftChannel.^2 ./ obj_filtered_Env(n).rightChannel.^2);
                
                ILD_distribution(position).Distribution(:,frequencyBand) = ILDvector;
                ILD_distribution(position).Label = obj(position).label;
                ILD_distribution(position).CenterFrequencies = analyzer.center_frequencies_hz;
                ILD_distribution(position).SampleRate = obj(position).fs;
            end
           ILD_distribution = reshape(ILD_distribution,size(obj));
        end
        
        % Short term signals
        function shortTermSignals = shortTerm(obj,timeFrame,overlap,win) % enter the times in milliseconds
                
            % shortTerm(obj, timeFrame, overlap, win) returns a
                % BinauralSignal object array where each element correspond
                % to the segmented signal at a specific time frame.
                %
                % Inputs: 
                %   - obj : BinauralSignal object
                %   - timeFrame : integer to set the length of each time frame (in ms)
                %   - overlap : integer to set the overlap between two
                %   contiguous time frames (in ms)
                %   - win : window type applied to each frame.
            if isscalar(obj)
                if nargin < 4
                    win = 'rectwin';
                end
                if nargin < 3
                    win = 'rectwin';
                    overlap = 0;
                end
              
                   % Converting durations from ms into samples
                   timeFrame_samples = round(obj.fs*timeFrame*1e-3);
                   overlap_samples = round(obj.fs*overlap*1e-3);
                   
                   left = buffer(obj.leftChannel,timeFrame_samples,overlap_samples);
                   right = buffer(obj.rightChannel,timeFrame_samples,overlap_samples);
                   
                   Nframes = size(left,2);
                   shortTermSignals(Nframes,1) = BinauralSignal();
                   
                   % Windowing
                   w = feval(win,timeFrame_samples);
                   
                   w = repmat(w,1,Nframes);
                   leftWindowed = left.*w;
                   rightWindowed = right.*w;
                   
                   
                   for frame = 1:Nframes
                       shortTermSignals(frame) = BinauralSignal([leftWindowed(:,frame) rightWindowed(:,frame)],obj.fs,[obj.label '/Frame ' num2str(frame)]);
                   end
                                      
            else
               error('The shortTerm method only apply on scalar BinauralSignal. Segmenting multidimensional array has not been developed yet.')
            end
        end
         
        function normObj = normalize(obj)
            N = numel(obj);
            normObj(N) = BinauralSignal;
            
            for n = 1:N
                
                normFactor = max(abs([obj(n).leftChannel; obj(n).rightChannel]));
                normObj(n).leftChannel = obj(n).leftChannel / normFactor;
                normObj(n).rightChannel = obj(n).rightChannel / normFactor;
                normObj(n).label = obj(n).label;
                normObj(n).fs = obj(n).fs;                
            end
        end
        
        % Plotting signals
        function plot(obj)
            figure, hold on
            plot(obj.leftChannel)
            plot(obj.rightChannel)
             
            hold off
        end
        
        % Playing binaural signal
        function play(obj)
           disp(['Playing binaural signal: ' obj.label]) 
           soundsc([obj.leftChannel obj.rightChannel],obj.fs)
        end
        
        % Writing binaural signal as audio file
        function write(obj, filename)
            
            N = numel(obj);
            for n = 1:N
                if nargin < 1
                    fname = [obj(n).label '.wav'];
                else                    
                    if strcmp(filename{n,1}(end),'\') % If only the destination folder has been passed
                        fname = [filename{n} obj(n).label '.wav'];
                    else                        
                        fname = [filename{n} '.wav'];
                    end
                    
                end
                
                y = horzcat(obj(n).leftChannel, obj(n).rightChannel);
                audiowrite(fname, y, obj(n).fs)
            end
        end
    end
    % Static methods -------------------------------
    methods (Static)
        % Get the spectrum of a signal
        function spectrum = computeSpectrum(temporalSignal)
            
                spectrum = fft(temporalSignal);

        end
        
        
        
            
       
    end
end

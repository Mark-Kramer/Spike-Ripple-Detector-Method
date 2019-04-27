% Spike-ripple detector.
% Developed by Catherine Chu, Arthur Chan, and Mark Kramer.
%
% INPUTS:
% data = the time series data, in this case EEG from one electrode.
% time = the time axis for the data, in units of seconds.
% varargin = set this paramater to a value between 0 and 1 to choose the
%   envelope threshold. If this parameter is not specified, the default
%   envelope threshold of 0.85 is used.
%
% OUTPUT:
% res = structure that holds each candidate spike-ripple event.
%   res.INPOS                %Start time [s]
%   res.FIPOS                %End time [s]
%   res.LEN                  %Duration [s]
%   res.freq                 %Frequency of ripple [Hz].
%   res.zc                   %Zero-crossing of ripple [Hz]
%   res.fano                 %Fano factor of ripple.
%   res.Lhite                %Difference between spike peak and start of interval.
%   res.Rhite                %Difference between spike peak and end of interval.
%   res.Ctime                %Time difference between start of ripple and start of spike.
%   res.Vpeak                %Peak voltage of spike.
% diagnostics = structure that holds method diagnostics (see code).
%
% DEPENDENCIES:
% findseq.m developed by Oleg Komarov.

function [res,diagnostics] = spike_ripple_detector(data,time, varargin)

      Fs = 1/(time(2)-time(1));                             %Sampling frequency.
      res = {};                                             %Structure to store results.

      fprintf('Design the filter ... \n')
      fNQ = Fs/2;                                        	%Define Nyquist frequeuncy.
      d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
        (60)/fNQ,...                                        %Freq @ edge first stop band.
        (100)/fNQ,...                                       %Freq @ edge of start of passband.
        (300)/fNQ,...                                       %Freq @ edge of end of passband
        (350)/fNQ,...                                       %Freq @ edge second stop band
        80,...                                              %Attenuation in the first stop band in decibels
        0.1,...                                         	%Amount of ripple allowed in the pass band.
        40);                                                %Attenuation in the second stop band in decibels
      Hd = design(d,'equiripple');                         	%Design the filter
      [num, den] = tf(Hd);                               	%Convert filter to numerator, denominator expression.
      order = length(num);                                	%Get filter order.
      diagnostics.num = num;
      diagnostics.den = den;
      diagnostics.order = order;
  
      dfilt = filter(num, den, data);                       %Filter the data.
      dfilt = [dfilt(floor(order/2+1):end); zeros(floor(order/2),1)];     %Shift after filtering.
      amp   = abs(hilbert(dfilt));                          %Compute amplitude envelope.
      
      if ~isempty(varargin)                                 %Choose envelope threshold,
          percentile_envelope = varargin{1};                %... as input,
      else                                                  %... or,
          percentile_envelope = 0.85;                       %... as default envelope threshold.
      end
      threshold = quantile(amp, percentile_envelope);       %Set amplitude threshold.
      fprintf(['Percentile envelope = ' num2str(percentile_envelope) ' ... \n'])
      diagnostics.threshold = threshold;                    %... save as diagnoistic to return.
      
                                                            %Get a sampling of max-start values.
      n_max_min = 10000;                                    %For 10,000 resamples,
      win_max_min = round(0.050*Fs);                        %... and window interval of 50 ms,
      max_min_distribution = zeros(n_max_min,1);            %... create a surrogate distribution,
      for n=1:n_max_min                                     %... for each surrogate,
          istart=floor(rand*(length(data)-win_max_min));    %... choose a random time index.
          if istart+win_max_min > length(data)              %... make sure the time interval is not too big,
              istart = length(data)-win_max_min-1;
          end
          if istart==0                                      %... make sure the time interval does not start at 0,
              istart=1;
          end                                               %... compute max value - value @ start of interval.
          max_min_distribution(n) = max(data(istart:istart+win_max_min-1))-data(istart);
      end
      

      percentile_max_and_peaks = 0.95;                      %Set max & peak threshold,
                                                            %Get threshold for max-start values.
      max_min_threshold = quantile(max_min_distribution, percentile_max_and_peaks);
                                                            %Get threshold for max voltage values.
      peak_threshold = quantile(data, percentile_max_and_peaks);

      above = find(amp > threshold);                        %Find amp's above threshold.
      t_separation = 0.005;                                 %Set small time seperation to 5 ms,
                                                            %... and merge small separations.
      small_separation = find(diff(above) > 1 & diff(above) < round(t_separation*Fs));
      for js=0:length(small_separation)-1
          ileft  = small_separation(end-js);
          iright = ileft+1;
          nleft  = above(ileft)+1;
          nright = above(iright)-1;
          above = [above(1:ileft); (nleft:nright)'; above(iright:end)];
      end
      
      [VALUES, INPOS, FIPOS, LEN] = findseq(diff(above));
      i_values=find(VALUES==1);                             %Locate sequences of value=1.
      diagnostics.sequences_above = length(i_values);
      
      if ~isempty(i_values)                                 %If we find sequences of 1's,
          INPOS = INPOS(i_values);                          %...save start index of each sequence,
          FIPOS = FIPOS(i_values);                          %...save end index of each sequence,
          LEN   = LEN(i_values);                            %...and save length.
          
          INPOS = time(above(INPOS));                       %Convert detections to TIME.
          FIPOS = time(above(FIPOS));
          LEN   = LEN/Fs;
          
          long_enough = find(LEN > 0.02);                   %Find intervals > 20 ms.
          diagnostics.number_long_enough = length(long_enough);
          
          if ~isempty(long_enough)                          %If we find intervals > 20 ms,
              fprintf(['Found sequences > 20 ms ... \n' ])
              INPOS=INPOS(long_enough);                     %... get those intervals.
              FIPOS=FIPOS(long_enough);
              LEN  =LEN(long_enough);
                                                            %Find intervals away from first/last time index.
              away_from_edges = find(INPOS > 1 & FIPOS < time(end)-1);
              INPOS=INPOS(away_from_edges);
              FIPOS=FIPOS(away_from_edges);
              LEN  =LEN(away_from_edges);
              
              % For each interval, compute zero-crossings, frequency, fano, left-right-max value, and time of max.
              zc   = zeros(length(INPOS),1);
              freq = zeros(length(INPOS),1);
              fano = zeros(length(INPOS),1);
              Lhite = zeros(length(INPOS),1);
              Rhite = zeros(length(INPOS),1);
              Ctime = zeros(length(INPOS),1);
              Vpeak = zeros(length(INPOS),1);
              for k=1:length(INPOS)                         % Find candidate HFO interval.
                  good = find(time >= INPOS(k) & time < FIPOS(k));
                  d0 = dfilt(good);                         % Get filtered data.
                  d0 = d0 - mean(d0);                       % Subtract mean.
                  d0((d0>0))=1;                             % Set values > 0 equal to 1.
                  d0((d0<0))=0;                             % Set values < 0 equal to 0.
                  zc0 = find(diff(d0)==1);                  % ZC when transition from 0 to 1.
                  zc(k)=length(zc0);                        % Count ZC.
                  ISI0 = diff(zc0);                         % Distance between ZC.
                  freq(k)=mean(1/(mean(ISI0)/Fs));          % Approx freq.
                  fano(k)=var(ISI0)/mean(ISI0);             % Fano factor of of zero crossing times.
                  center_time = mean(time(good));           % Center time of window.
                  good = find(time >= center_time-0.05 & time < center_time+0.05);  %+/- 50 ms around center.
                  dorig = data(good);                       % Get unfiltered data.
                  dorig = smooth(dorig,11);                 % Smooth it.
                  [mx, imx] = max(dorig);                   % Find max.
                  Lhite(k) = mx-dorig(1);                   % Difference between max & left (or start) of interval.
                  Rhite(k) = mx-dorig(end);                 % Difference betweem max & right (or end) of interval.
                  Ctime(k) = INPOS(k) - time(good(imx));    % Time from ripple start to max
                  Vpeak(k) = mx;                            % Max values.
              end
              
              % Find intervals that pass tests.
              threshold_fano = 1;                           %Fix Fano threshold.

              %To classify as a spike-ripple detection, must have:
              good = find(zc >= 3 ...                       % At least 3 ZC.
                  & fano < threshold_fano ...               % Fano < 1.
                  & Lhite > max_min_threshold ...           % Max - start value > threshold.
                  & Rhite > max_min_threshold ...           % Max - end value > threshold.
                  & Ctime < 0 ...                           % Ripple begin before peak.
                  & Vpeak > peak_threshold);                % Max > threshold.
              end
              % Save candidate spike-ripple events.
              INPOS=INPOS(good);
              FIPOS=FIPOS(good);
              LEN  =LEN(good);
              zc   = zc(good);
              freq = freq(good);
              fano = fano(good);
              Lhite = Lhite(good);
              Rhite = Rhite(good);
              Ctime = Ctime(good);
              Vpeak = Vpeak(good);
              diagnostics.number_detections = length(good);
              
              %Sort detections by starting time.
              [~, isort] = sort(INPOS, 'ascend');
              INPOS = INPOS(isort);
              FIPOS = FIPOS(isort);
              LEN   = LEN(isort);
              freq  = freq(isort);
              zc    = zc(isort);
              fano  = fano(isort);
              Lhite  = Lhite(isort);
              Rhite  = Rhite(isort);
              Ctime  = Ctime(isort);
              Vpeak  = Vpeak(isort);
              fprintf(['Candidate spike-ripple events = ' num2str(length(good)) ' ... \n' ])
              
              %Save the results for each candidate spike-ripple event.
              res.INPOS = INPOS;                %Start time [s]
              res.FIPOS = FIPOS;                %End time [s]
              res.LEN   = LEN;                  %Duration [s]
              res.freq  = freq;                 %Frequency of ripple [Hz].
              res.zc    = zc;                   %Zero-crossing of ripple [Hz]
              res.fano  = fano;                 %Fano factor of ripple.
              res.Lhite  = Lhite;               %Difference between spike peak and start of interval.
              res.Rhite  = Rhite;               %Difference between spike peak and end of interval.
              res.Ctime  = Ctime;               %Time difference between start of ripple and start of spike.
              res.Vpeak = Vpeak;                %Peak voltage of spike.
          end
end
% Ripple detector.
%
% INPUTS:
% data = the time series data, in this case EEG from one electrode.
% time = the time axis for the data, in units of seconds.
% ADVANCED INPUTS:
% varargin = 'PercentileEnvelope', value
%   set 'value'  between 0 and 1 to choose the
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
% diagnostics = structure that holds method diagnostics (see code).
%
% DEPENDENCIES:
% findseq.m developed by Oleg Komarov.

function [res,diagnostics] = ripple_detector(data,time, varargin)

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
      if any(isnan(dfilt))                                  % If there's a nan in the data,
          nan_dfilt = find(isnan(dfilt));                   % ... find nans,
          temp_for_hilbert = dfilt;                         % ... make dfilt for hilbert,
          temp_for_hilbert(nan_dfilt)=0;                    % ... and put 0s at nans
          amp   = abs(hilbert(temp_for_hilbert));                                                                                       %Compute amplitude envelope.
      else
          amp   = abs(hilbert(dfilt));                      %Compute amplitude envelope.
      end
      
      if any(strcmp(varargin, 'PercentileEnvelope'))        %Choose envelope threshold (ADVANCED)
          i0 = 1+find(strcmp(varargin, 'PercentileEnvelope'));
          percentile_envelope = varargin{i0};               %... as input,
      else                                                  %... or,
          percentile_envelope = 0.85;                       %... as default envelope threshold.
      end
      threshold = quantile(amp, percentile_envelope);       %Set amplitude threshold.
      fprintf(['Percentile envelope = ' num2str(percentile_envelope) ' ... \n'])
      diagnostics.threshold = threshold;                    %... save as diagnoistic to return.

      binary_above = amp > threshold;
      above = find(amp > threshold);                        %Find amp's above threshold.
      t_separation = 0.005;                                 %Set small time seperation to 5 ms,
                                                            %... and merge small separations.
      small_separation = find(diff(above) > 1 & diff(above) < round(t_separation*Fs));
      for js=0:length(small_separation)-1
          ileft  = small_separation(end-js);
          iright = ileft+1;
          binary_above(above(ileft):above(iright))=1;
      end
      above=find(binary_above);
      
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
              parfor k=1:length(INPOS)                         % Find candidate HFO interval.
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
              end
              
              % Find intervals that pass tests.
              threshold_fano = 1;                           %Fix Fano threshold.

              %To classify as a spike-ripple detection, must have:
              good = find(zc >= 3 ...                       % At least 3 ZC.
                  & fano < threshold_fano);                 % Fano < 1.
  
              end
              % Save candidate spike-ripple events.
              INPOS=INPOS(good);
              FIPOS=FIPOS(good);
              LEN  =LEN(good);
              zc   = zc(good);
              freq = freq(good);
              fano = fano(good);
              diagnostics.number_detections = length(good);
              
              %Sort detections by starting time.
              [~, isort] = sort(INPOS, 'ascend');
              INPOS = INPOS(isort);
              FIPOS = FIPOS(isort);
              LEN   = LEN(isort);
              freq  = freq(isort);
              zc    = zc(isort);
              fano  = fano(isort);
              fprintf(['Candidate spike-ripple events = ' num2str(length(good)) ' ... \n' ])
              
              %Save the results for each candidate spike-ripple event.
              res.INPOS = INPOS;                %Start time [s]
              res.FIPOS = FIPOS;                %End time [s]
              res.LEN   = LEN;                  %Duration [s]
              res.freq  = freq;                 %Frequency of ripple [Hz].
              res.zc    = zc;                   %Zero-crossing of ripple [Hz]
              res.fano  = fano;                 %Fano factor of ripple.
          end
end
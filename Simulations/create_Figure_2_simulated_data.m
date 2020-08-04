% Create the simulated data in Figure 2 of [Chu et al, J Neurosci Methods, 2017]

function [d0,t0] = create_Figure_2_simulated_data(simulation)

  Fs = 2035;              %Sampling rate.
  T  = 10*60;             %Total s.
    
  %%%%%%%%%% PINK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if strcmp(simulation, 'Pink')
      fprintf(['Running Pink simulation ... \n'])
      d0 = make_pink_noise(0.5,T*Fs,1/Fs);
  end
  
  %%%%%%%%%% PINK+PULSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if strcmp(simulation, 'Pink+Pulse')
      fprintf(['Running Pink+Pulse simulation ... \n'])
      d0 = make_pink_noise(0.5,T*Fs,1/Fs);
      t = 1/Fs : 1/Fs : 600;
      d = 1/Fs : 1    : 600;
      y = pulstran(t,d,'tripuls',0.05,0);
      y = y / max(y);
      [~, i0] = findpeaks(y);
      for i=1:length(i0)
          y(i0(i)-100:i0(i)+100) = 20*std(d0)*y(i0(i)-100:i0(i)+100);
      end
      d0 = d0 + y';
  end

  %%%%%%%%%% PINK+SPIKE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if strcmp(simulation, 'Pink+Spike')
      fprintf(['Running Pink+Spike simulation ... \n'])
      load('average_spike.mat')
      Fs = round(Fs);
      spike = spike - mean(spike);
      spike = spike / max(spike);
      d0 = make_pink_noise(0.5,T*Fs,1/Fs);
      std0 = std(d0);
      for i=1:600
          i0 = round(Fs/2) + round((i-1)*Fs);
          d0(i0:i0+length(spike)-1) = d0(i0:i0+length(spike)-1) + 20*std0*spike;
      end
  end
  
  %%%%%%%%%% PINK+SPIKE+HFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if strcmp(simulation, 'Pink+Spike+HFO')
      fprintf(['Running Pink+Spike+HFO simulation ... \n'])
      load(['average_spike.mat']);
      Fs = round(Fs);
      spike = spike - mean(spike);
      spike = spike / max(spike);
      spike = spike .*hann(length(spike));
      spike0 = spike;
      d0 = make_pink_noise(0.5,T*Fs,1/Fs);
      std0 = std(d0);
      
      for i=1:600
          spike = spike0;
          HFO = randn(Fs,1);
          Wn = [110, 120]/(Fs/2);	%...set the passband,
          n  = 100;					%...and filter order,
          b  = fir1(n,Wn);			%...build bandpass filter.
          HFO = filtfilt(b,1,HFO);	%...and apply filter.
          HFO = HFO(1000:1000+round(0.05*Fs)-1);
          HFO = hann(length(HFO)).*HFO;
          HFO = HFO/std(HFO);
          HFO = 0.05*std(spike)*HFO;
          istart = 190;
          spike(istart:istart+length(HFO)-1) = spike(istart:istart+length(HFO)-1) + 1*HFO;
          i0 = round(Fs/2) + round((i-1)*Fs);
          d0(i0:i0+length(spike)-1) = d0(i0:i0+length(spike)-1) + 20*std0*spike;
      end
  end
  
  %%%%%%%%%% PINK+SPIKE+HFO+30% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if strcmp(simulation, 'Pink+Spike+HFO+30%')
      fprintf(['Running Pink+Spike+HFO+30% simulation ... \n'])
      
      load(['average_spike.mat']);
      Fs = round(Fs);
      spike = spike - mean(spike);
      spike = spike / max(spike);
      spike = spike .*hann(length(spike));
      spike0 = spike;
      d0 = make_pink_noise(0.5,T*Fs,1/Fs);
      std0 = std(d0);
      
      for i=1:600
          spike = spike0;
          i0 = round(Fs/2) + round((i-1)*Fs);
          
          if mod(i,3)==0
              HFO = randn(Fs,1);
              Wn = [110, 120]/(Fs/2);			%...set the passband,
              n  = 100;                         %...and filter order,
              b  = fir1(n,Wn);                  %...build bandpass filter.
              HFO = filtfilt(b,1,HFO);          %...and apply filter.
              HFO = HFO(1000:1000+round(0.05*Fs)-1);
              HFO = hann(length(HFO)).*HFO;
              HFO = HFO/std(HFO);
              HFO = 0.05*std(spike)*HFO;
              istart = 190;
              spike(istart:istart+length(HFO)-1) = spike(istart:istart+length(HFO)-1) + 1*HFO;
          end
          
          d0(i0:i0+length(spike)-1) = d0(i0:i0+length(spike)-1) + 20*std0*spike;
      end
  end
  
  t0 = (1:length(d0))/Fs;
  
end

function [x1new] = make_pink_noise(alpha,L,dt)

  x1 = randn(L,1);
  xf1 = fft(x1);
  A = abs(xf1);
  phase = angle(xf1);

  df = 1.0 / (dt*length(x1));
  faxis = (0:length(x1)/2)*df;
  faxis = [faxis, faxis(end-1:-1:2)];  %(end-1:-1:2)
  oneOverf = 1.0 ./ faxis.^alpha;
  oneOverf(1)=0.0;

  Anew = A.*oneOverf';
  xf1new = Anew .* exp(i*phase);
  x1new = real(ifft(xf1new));
  
end

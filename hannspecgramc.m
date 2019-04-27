function [S,times,freq]=hannspecgramc(d1,movingwin,params)

  Fs = params.Fs;
  tEDF = params.tEDF;
  flo  = params.fpass(1);
  fhi  = params.fpass(2);

  iStart = 1;
  iWin   = round(movingwin(1)*Fs/2)*2;
  iStep  = round(movingwin(2)*Fs/2)*2;
  T      = iWin/Fs;
  df     = 1/T;
  fNQ    = Fs/2;
  freq   = (0:iWin/2)*df;
  counter=1;
  S      = zeros(ceil(length(tEDF)/iStep), iWin/2+1);
  times  = zeros(ceil(length(tEDF)/iStep), 1);
  while iStart+iWin < length(tEDF)

      dnow    = d1(iStart:iStart+iWin-1);
      dnow    = dnow - mean(dnow);
      dnow    = hann(length(dnow)).*dnow;
      spectrum= fft(dnow).*conj(fft(dnow));
      S(counter,:)   = spectrum(1:iWin/2+1);
      times(counter) = tEDF(iStart+iWin/2);
      
      counter=counter+1;
      iStart = iStart + iStep;
  end
  S = S(1:counter-1,:);
  times = times(1:counter-1);
  fgood = find(freq >= flo & freq <= fhi);
  S = S(:,fgood);
  freq = freq(fgood);

end
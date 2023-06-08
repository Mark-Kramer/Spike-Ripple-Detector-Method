% Function to visualize and classify the candidate spike gamma events
% detected in spike_gamma_detector.m
%
% This function produces a figure showing the (1) the original data, (2)
% the filtered data, and (3) the spectrogram surrounding each candidate spike gamma
% event.
%
% For each candidate spike gamma event, the user enters "y" for YES or "n"
% for NO at the command line. The results of these classifications are then
% returned in the output.

% INPUTS.
% data = the LFP data (same as for spike_gamma_detector.m)
% time = the time axis (same as for spike_gamma_detector.m)
% res0 = the 1st output of spike_gamma_detector.m
% diagnostics = the 2nd output of spike_gamma_detector.m

% OUTPUTS.
% expert_classify = the classification of each candidate spike gamma
%   event as 'y' or 'n'. 

function [expert_classify] = spike_gamma_visualizer(data, time, res0, diagnostics)

    Fs = 1/(time(10)-time(9));
    num = diagnostics.num;
    den = diagnostics.den;
    order = diagnostics.order;
    
    dfilt = filter(num, den, data);                                    %Filter it.
    dfilt = [dfilt(floor(order/2)+1:end); zeros(floor(order/2),1)];   %Shift after filtering.
    
    n_detections = length(res0.INPOS);
    counter = 1;
    for k=1:n_detections
        INPOS = res0.INPOS(k);
        FIPOS = res0.FIPOS(k);
        LEN   = res0.LEN(k);
        freq  = res0.freq(k);
        zc    = res0.zc(k);
        fano  = res0.fano(k);
        Lhite = res0.Lhite(k);
        Rhite = res0.Rhite(k);
        Ctime = res0.Ctime(k);
        Vpeak = res0.Vpeak(k);
        
        igood = find(time > INPOS-0.5 & time < FIPOS+0.5);
        t     = time(igood);
        dat   = data(igood);
        datf  = dfilt(igood);
        
        dspec = data(igood)-mean(data(igood));
        params.Fs    = Fs;             % Sampling frequency [Hz]
        params.fpass = [30 100];       % Frequencies to visualize in spectra [Hz]
        movingwin    = [0.200,0.005];  % Window size, Step size [s]
        params.tEDF  = t;
        [S,S_times,S_freq] = hannspecgramc(dspec,movingwin,params);
        %Smooth the spectra.
        t_smooth = 11;
        dt_S     = S_times(2)-S_times(1);
        myfilter = fspecial('gaussian',[1 t_smooth], 1);
        if k==1
            fprintf(['Smooth spectra over +/- ' num2str((t_smooth-1)/2*(dt_S)*1000,3) ' ms \n'])
        end
        S_smooth = imfilter(S, myfilter, 'replicate');   % Smooth the spectrum.
        
        res{counter}.INPOS = INPOS;
        res{counter}.FIPOS = FIPOS;
        res{counter}.LEN   = LEN;
        res{counter}.freq  = freq;
        res{counter}.zc    = zc;
        res{counter}.fano  = fano;
        res{counter}.Lhite = Lhite;
        res{counter}.Rhite = Rhite;
        res{counter}.Ctime = Ctime;
        res{counter}.Vpeak = Vpeak;
        res{counter}.threshold = diagnostics.threshold;
        res{counter}.t     = t;
        res{counter}.data  = dat;
        res{counter}.dfilt = datf;
        res{counter}.S_smooth = S_smooth;
        res{counter}.S_times  = S_times;
        res{counter}.S_freq   = S_freq;
        res{counter}.k     = k;
        
        counter = counter + 1;
    end

    %% Visualize

    expert_classify = cell(length(res),1);
    
    ek = 1;
    while ek <= length(res)
        
        INPOS = res{ek}.INPOS;
        FIPOS = res{ek}.FIPOS;
        LEN   = res{ek}.LEN;
        freq  = res{ek}.freq;
        zc    = res{ek}.zc;
        fano  = res{ek}.fano;
        Lhite = res{ek}.Lhite;
        Rhite = res{ek}.Rhite;
        Ctime = res{ek}.Ctime;
        Vpeak = res{ek}.Vpeak;
        
        data  = res{ek}.data;
        dfilt = res{ek}.dfilt;
        amp   = abs(hilbert(dfilt));
        t     = res{ek}.t;
        threshold = res{ek}.threshold;
        
        S_times   = res{ek}.S_times;
        S_freq    = res{ek}.S_freq;
        S_smooth  = res{ek}.S_smooth;
                
        t0 = 0;
        t = t-t0;
        INPOS = INPOS - t0;
        FIPOS = FIPOS - t0;
        S_times = S_times - t0;
        
        subplot(3,1,1)
        plot(t, data)
        ax = axis;
        axis tight
        hold on
        plot([INPOS INPOS], [ax(3) ax(4)], 'k')
        plot([FIPOS FIPOS], [ax(3) ax(4)], 'k')
        hold off
        title(['# ' num2str(length(INPOS)) ', @ ' num2str(k) ...
            ' DUR ' num2str(LEN*1000,3)  ', FQ ' num2str(freq,3) ...
            ', LHT ' num2str(Lhite,3)  ', RHT ' num2str(Rhite,3) ', CT ' num2str(Ctime,3)...
            ', VPK ' num2str(Vpeak,3) ...
            ', ZC ' num2str(zc)     ', FF ' num2str(fano,2)])
        
        subplot(3,1,2)
        plot(t, dfilt)
        hold on
        plot(t, threshold*ones(size(t)), 'k')
        axis tight
        ax = axis;
        plot([INPOS INPOS], [ax(3) ax(4)], 'k')
        plot([FIPOS FIPOS], [ax(3) ax(4)], 'k')
        hold off
        
        %Compute the spectrogram.
        subplot(3,1,3)
        colormap(jet)
        imagesc(S_times,S_freq,log10(S_smooth)')  %, [-3 -0.5])
        axis xy
        ax = axis;
        hold on
        plot([INPOS INPOS], [ax(3) ax(4)], 'k')
        plot([FIPOS FIPOS], [ax(3) ax(4)], 'k')
        hold off
        xlim([INPOS-0.5, FIPOS+0.5])
        
        input0 = input(['Event ' num2str(ek) ' of ' num2str(length(res)) ', Is there an spike-gamma? y/[n]/b: '], 's');
        switch input0
            case 'b'
               ek = ek-1; fprintf(['going back one to ' num2str(ek) '\n'])
            case ''                                       %If you enter nothing,
               expert_classify{ek,1}         = 'n';       %... it's "no".
               ek = ek+1;
            otherwise
               expert_classify{ek,1}         = input0;
               ek = ek+1;
        end
    end
    
end
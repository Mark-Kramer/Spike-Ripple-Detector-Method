# Spike-Gamma-Detector-Method

Analysis code for detection of candidate spike-gamma events.


---
---
## Beta-version only. Not fully debugged.
---
---

**Extend** the spike-ripple detector to detect spike-**gamma** events.

The gamma events (see [Ren et. al., Neurology 2015](https://pubmed.ncbi.nlm.nih.gov/25589669/) and [Thomas et. al., Ann Neurology, 2023](https://pubmed.ncbi.nlm.nih.gov/36373178/)):
- lower frequency (30-100 Hz), 
- precede spikes.

Modifications from the original spike-ripple code include:
- Filtering approximately 30-100 Hz.
- Larger interval (400 ms) around events.
- Most (75%) of the gamma event must precede the spike maximum.
- Fano threshold = 2

The analysis code is located in the file:  <code>spike_gamma_detector.m</code>.  This function accepts two inputs which specify the data (i.e., a single channel of EEG) and a time axis in units of seconds.  A third optional input may also be included to specify the envelope threshold (default value is 0.85).  This function returns two outputs. The first output is a structure, in which each entry contains the features of a candidate spike-gamma event.  The second output returns additional diagnostics. This output may then be used as input to visualize and classify the candidate spike gamma events in the file: <code>spike_gamma_visualizer.m</code>

For example:

<pre><code> data = EEG;  %where "EEG" corresponds to voltage data recorded from one scalp electrode.
 Fs = 1000;  %the sampling frequency in Hz.
 t  = (1:length(EEG))/Fs; %corresponding time axis.
 [res,diagnostics] = spike_gamma_detector(EEG,t);  %Call the function, and return the candidate spike-gamma events.
 [expert_classify] = spike_gamma_visualizer(EEG,t,res,diagnostics)	% Visualize and classify the candidate spike-gamma events.
 </code></pre>

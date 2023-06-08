# Spike-Ripple-Detector-Method

Analysis code for detection of candidate spike-ripple events.

This MATLAB repository contains the data analysis method developed in:

<i>C. J. Chu et al, "A semi-automated method for rapid detection of ripple events on interictal voltage discharges in the scalp electroencephalogram", Journal of Neuroscience Methods, Volume 277, 1 February 2017, Pages 46–55.</i>

The analysis code is located in the file:  <code>spike_ripple_detector.m</code>.  This function accepts two inputs which specify the data (i.e., a single channel of EEG) and a time axis in units of seconds.  A third optional input may also be included to specify the envelope threshold (default value is 0.85).  This function returns two outputs. The first output is a structure, in which each entry contains the features of a candidate spike-ripple event.  The second output returns additional diagnostics. This output may then be used as input to visualize and classify the candidate spike ripple events in the file: <code>spike_ripple_visualizer.m</code>

For example:

<pre><code> data = EEG;  %where "EEG" corresponds to voltage data recorded from one scalp electrode.
 Fs = 1000;  %the sampling frequency in Hz.
 t  = (1:length(EEG))/Fs; %corresponding time axis.
 [res,diagnostics] = spike_ripple_detector(EEG,t);  %Call the function, and return the candidate spike-ripple events.
 [expert_classify] = spike_ripple_visualizer(EEG,t,res,diagnostics)	% Visualize and classify the candidate spike-ripple events.
 </code></pre>

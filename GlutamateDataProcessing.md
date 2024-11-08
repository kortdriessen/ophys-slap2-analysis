## Glutamate data preprocessing

The top-level functions that process glutamate imaging data are `processSLAP2.m` and `processBergamo.m` depending on which microscope the data were acquired on. These functions generate output in a shared format. They can be called with no arguments, in which case a GUI will gather user input. The pipeline consists of three steps, which the top-level function calls in sequence:
1. Defining the trial structure
    * `buildTrialTableSLAP2` or `builtTrialTableBergamo.m`
2. Registering the images within each recording:  
    * ``multiROIRegSLAP2`` or ``stripRegBergamo.m``
3. Extracting signals and generating a "Summary" file for the experiment:
    * `summarize_NoLoCo.m`
  
If processing aborts, you can call the top-level function again and it will skip over any completed steps (unless you set an option flag to repeat them).

When you call these functions without parameters, they pop up a GUI with default parameters that you can change. The default parameters for each function are set in the function `setParams.m`. There are tooltips that explain each parameter when you hover over it.

### 1. Defining the trial structure
This step asks the user to identify the files that are part of the experiment, and builds a **trial table** that keeps track of files generated during processing. The experiment is assumed to consist of one or more **epochs**, each consisting of one or more **trials**. All recordings in an experiment are assumed to image the same field of view, and must be in the same folder. The recordings will be aligned to each other and share a common set of extracted synapses; if you image different fields of view in a single session, these should be split into separate folders and processed separately. 

For SLAP2 recordings, the corresponding SLAP2 reference stack should exist in the same folder as the activity recording, or in a subfolder. Do not include multiple reference stacks in the folder; this will cause an error.

When you run this function, a list of files will be shown, and you will be asked to select files belonging to each epoch, in sequence. Files not assigned to an epoch will not be analyzed.

**Bergamo recordings:** If all of your files belong to a single epoch you can automatically generate the trial table by directly calling `stripRegBergamo.m`. This simplifies batch analysis without GUI input.

### 2. Registering Images
Image registration aligns the images within each trial to each other and generates a downsampled, aligned copy of each recording as a `.tif` file. It also generates an associated `..._ALIGNMENTDATA.mat` file for each `.tif`, which contains data related to the alignment process and is used by later steps.  

### 3. Extracting signals
We represent signals with a matrix factorization model (_i.e._ as a set of spatial patterns that each vary in brightness over time). This approach is shared with common approaches to analyze somatic calcium imaging data (_e.g._ using CaImAn or Suite2P). Our approach has been tailored specifically for synaptic imaging with the following assumptions:
* Sources are small and round (with characteristic size `sigma_px`)
* Sources flash with positive-going transients that rise quickly and have exponential decay time constant `tau_s`.
* Sources frequently overlap with neighboring sources, requiring a superresolution approach to separate them (_e.g._ SOFI or localization microscopy).

Usin the downsampled recordings from step 2, for each trial, we filter the movie in space and time to emphasize sources that match the priors above, subtract out any slow fluctuations, and  compute a high-order statistic (similar to SOFI) to produce an image that sharpens and enhances fluctuating point sources. This image `actIM`, is the activity image for each trial. This computation is done by the function `localizeSources.m` 

We then align all trials in the experiment, averaging together the activity images to create a single activity image for the experiment. Putative synapses are identified as local maxima in this image.

Still using the downsampled recordings, we then form a [pixels x time] matrix from the pixels surrounding the putative synapses across all recordings, and factorize this using constrained nonnegative matrix factorization while enforcing our spatial and temporal priors. During this process, we merge sources that overlap and do not show sufficiently distinct activity (using the Akaike Information Criterion) and discard sources that do not explain sufficient variance in the recordings. 

Finally, we use the resulting footprints of these sources to extract activity from the original, high-framerate trial recordings. We extract different versions of the signal with different temporal structure:  
* dFraw: Fluctuations in fluorescence with the temporal prior enforced
* matchFiltered: Fluctuations in fluorescence, with an optimal linear filter for detection applied
* spikes: Fluctuations in fluorescence, deconvolved, so release events are brief spikes
* dFls: Fluctuations in Fluorescence with no temporal prior enforced

We also compute F0 (the time-varying fluorescence baseline for each source) and dF/F0 (_i.e._ fractional fluorescence change, for dFraw and dFls.

## Tips
**Parallel Processing.** The pipeline uses matlab's parallel computing toolbox to greatly speed up processing. The optimal number of parallel workers to use depends on the computer's CPU and memory capacity and the length of your recordings, is different for the two major steps (Registration and Signal Extraction), and is also different for SLAP2 vs Bergamo analysis because of the amount of data in memory at once. Errors and significant slowdowns can occur if you run out of memory due to too many parallel workers. When you start processing, you may want to monitor memory and CPU usage, adjust the number of workers accordingly, and then restart processing.

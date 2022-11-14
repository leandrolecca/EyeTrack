# EyeTrack
EyeTrack is an R-based toolbox intended to facilitate research analysis of eyetracking data (see [How it works](#how-it-works)). It has been programmed in a modular and scalable fashion to open collaborative work between researchers interested in implementing different algorithms for saccade, blink and fixation detection while respecting a main architecture (see [Pipeline architecture](#pipeline-architecture)).

## How it works
From the x-y gaze position of a trial, EyeTrack provides event information about saccades, blinks and fixations that are taking place. In the main.md notebook there is an example of how to run the toolbox after giving a minimal structure format to the eyetracking data, aka *etr.db* (see [DataFrame Structure](#dataframe-structure)). Before, the user must provide a set of parameters to make work the detection, which are:

* MINDUR: minimum duration of a saccade, in millisecons (ms)
* VFAC: velocity threshold function
* SAMPLING: sampling rate at which signal was collected in Hertz
* blink.threshold: minimum duration of a blink in ms
* blink.window: time window to remove around the blinks in ms

``` r
MINDUR <- 10
VFAC <- seq(6,10)
SAMPLING <- 1000
blink.threshold <- 20
blink.window <- 100
```

The for loop passes through each participant along the *etr.db*. This loop can be more complex, such as including different VFAC levels, runs or conditions, but ensure that a subset of trials is sent to the **run.participant()** function and the iterative information is stored in the *result.db* variable. Upon completion, the results should be visually examined to roughly assess the quality of the analysis (see [EyeTrack Visualizer](#eyetrack-visualizer)). Here is the dataset to run the code sample below.

``` r
SUBJ <- c("subS01","subS02")

# define an empty list to save the summary, fixations, saccades and blinks
#-------------------------------------------------------------------------------
result.db <- list(summary=c(),fixations=c(),blinks=c(),saccades=c())

# Now I use the run.participant function to get the events (blinks, fixations,
# saccades) on every trial for each eye
#-------------------------------------------------------------------------------
for (sbj in 1:length(SUBJ)){
  for (vfc in 1:length(VFAC)){
    participant <- etr.db[etr.db$subID == SUBJ[sbj],]
    right.eye <- run.participant(participant,eye='right',SAMPLING,
                                 blink.threshold,blink.window,MINDUR,VFAC[vfc])
    left.eye <- run.participant(participant,eye='left',SAMPLING, 
                                blink.threshold,blink.window,MINDUR,VFAC[vfc])
    
    result.eyes <- mapply(rbind, left.eye, right.eye,SIMPLIFY=FALSE)
    result.eyes <- mapply(cbind, result.eyes, "VFAC"=VFAC[vfc], SIMPLIFY=FALSE)
    result.eyes <- mapply(cbind, result.eyes, "MINDUR"=MINDUR, SIMPLIFY=FALSE)
    result.eyes <- mapply(cbind, result.eyes, "subID"=SUBJ[sbj], SIMPLIFY=FALSE)
    
    result.db <- mapply(rbind,result.db,result.eyes,SIMPLIFY=FALSE)
  } 
}

save(list=c("result.db"), file=file.path(basedir, 'eyetrack-report.RData'))
```

## EyeTrack Visualizer
The **plt.sbj.analysis()** function plots the time series of the x-y position of a trial and marks the events detected along the signal in black, blue, yellow and red for fixations, saccades, blinks and artifacts, respectively. 

``` r
sbj <- "sub009"
plt.trial <- plt.sbj.analysis(sbj,vfac=6,SAMPLING=SAMPLING,etr.db=etr.db,
                              trl.idx=390,EYE='right',result.db=result.db)
# show raw signal
plt.trial$plt.raw
# show events
plt.trial$plt.analysis
```

## Pipeline architecture
The pipeline begins with a call to **run.participant()**, which takes a subset of trials for a single subject and calls **saccfixblink()** therein for each trial (if the trial does not account for at least 25 % of its total signal, this is interpreted as a high eyetracking loss, and the analysis for that trial is not performed). The saccfixblink() function first calls the **microsacc()** function for saccade detection, taken and adapted from [Engbert’s microsaccade toolbox](http://read.psych.uni-potsdam.de/index.php?option=com_content&view=article&id=140:engbert-et-al-2015-microsaccade-toolbox-for-r&catid=26:publications&Itemid=34) [1], and whose algorithm is based on a velocity threshold (see [here](https://reader.elsevier.com/reader/sd/pii/S0042698903000841?token=D920381623BEBD3293EFA0C66393604FA29032371144D8C9E4AEBA121ED09967D2BE5A4A9209C85430377A11CE466C18&originRegion=eu-west-1&originCreation=20221114081753) for more details). Besides, blinks are detected with the original **blink.detect()** function, which looks for a series of NA values (e.g., loss of eye tracking due to eyelid closure) that is a minimum time in milliseconds specified in the blink.threshold parameter, and then removes a time window specified in the blink.window parameter around the ends of the signal loss, since eyelid closure and reopening typically cause spurious movements that lead to artifacts in the signal. After both saccades and blinks are detected, saccades inside the blink window range, e.g., at the beginning and end, are removed. If the onset of the saccade is outside but the offset is inside the blink range, the start of the blink becomes the onset of the saccade and the saccade is removed. If the onset of the saccade is inside but the offset is outside the blink range, the start of the blink becomes the offset of the saccade and the saccade is removed. Finally, fixations and artifacts are detected and labeled with the **aggregate.fixations()** and **label.fixations()** functions taken from [Tmalsburg's toolboox](https://github.com/tmalsburg/saccades), respectively, from the remaining signal without event detection. Therefore, each time step of the signal should be labeled with a unique event. The information about the beginning and end of each event is stored in results.db$summary, while more detailed information about saccades, blinks and fixations can be found in result.db$saccades, result.db$blinks and result.db$fixations.

The workflow between [core functions](#core-functions-information) is depicted below.

## Core functions information
The core functions can be edited or even new functions can be added with different detection heuristics as long as they adhere to the output structure (see full code in eyetrack_toolbox.R).

## DataFrame Structure
Running the toolbox requires to give specific column names to the “etr.db” dataframe, which are:

* lxdeg: left eye x position in degs
* lydeg: left eye y position in degs
* rxdeg: right eye x position in degs
* rydeg: right eye y position in degs
* TRIAL_INDEX: int number for the trial
* RECORDING_SESSION_LABEL: label for the trial
* subID: label of the subject

This data can be taken from the .edf files that EyeLink saves after running your task script. To do this, open the .edf file with [DataViewer](https://www.sr-research.com/data-viewer/) and save an Excel file in which you select the variables to be exported under Analysis -> Report -> Sample report. The main variables to export are RECORDING_SESSION_LABEL, TRIAL_INDEX, SAMPLE_MESSAGE, LEFT_GAZE_X, LEFT_GAZE_Y, RIGHT_GAZE_X and RIGHT_GAZE_Y. Since the left and right gaze are in pixels, they must be converted to visual angles and then renamed to lxdeg and lydeg for the left eye and rxdeg and rydeg for the right eye, corresponding to the X-Y components. Finally, the column RECORDING_SESSION_LABEL must contain a unique label for each trial, which must then be extracted from the SAMPLE_MESSAGE, as it must contain specific information about the trial at the correct time (this is programmed by the user when building the stimuli presentation script, e.g., in PsychToolbox with *Eyelink('Message', 'Your message')* command). See [here](#sample) a 10 ms sample.

### Sample
```
        subID	RECORDING_SESSION_LABEL	                                TRIAL_INDEX	    lxdeg       lydeg           rxdeg           rydeg
6000000	subS02	58_seq2_HML_L_run6_Q3_right_x=-3.2361º_y=-2.3511º.png	    1018	    -0.1586261	0.4833037	-0.0029375	0.2354599
6000001	subS02	58_seq2_HML_L_run6_Q3_right_x=-3.2361º_y=-2.3511º.png	    1018	    -0.1630324	0.4874343	-0.0029375	0.2327060
6000002	subS02	58_seq2_HML_L_run6_Q3_right_x=-3.2361º_y=-2.3511º.png	    1018	    -0.1645011	0.4956954	-0.0029375	0.2299521
6000003	subS02	58_seq2_HML_L_run6_Q3_right_x=-3.2361º_y=-2.3511º.png	    1018	    -0.1615636	0.5163483	 0.0220315	0.2258213
6000004	subS02	58_seq2_HML_L_run6_Q3_right_x=-3.2361º_y=-2.3511º.png	    1018	    -0.1512823	0.5328705	 0.0484692	0.2230674
6000005	subS02	58_seq2_HML_L_run6_Q3_right_x=-3.2361º_y=-2.3511º.png	    1018	    -0.1410010	0.5480158	 0.0514067	0.2216904
6000006	subS02	58_seq2_HML_L_run6_Q3_right_x=-3.2361º_y=-2.3511º.png	    1018	    -0.1615636	0.5604073	 0.0264377	0.2203135
6000007	subS02	58_seq2_HML_L_run6_Q3_right_x=-3.2361º_y=-2.3511º.png	    1018	    -0.1880011	0.5700452	 0.0088126	0.2161827
6000008	subS02	58_seq2_HML_L_run6_Q3_right_x=-3.2361º_y=-2.3511º.png	    1018	    -0.2217822	0.5590305	 0.0132189	0.2024132
6000009	subS02	58_seq2_HML_L_run6_Q3_right_x=-3.2361º_y=-2.3511º.png	    1018	    -0.2276572	0.5383779	 0.0146876	0.1776281
6000010	subS02	58_seq2_HML_L_run6_Q3_right_x=-3.2361º_y=-2.3511º.png	    1018	    -0.2203135	0.5177252	 0.0117501	0.1555968
...
```

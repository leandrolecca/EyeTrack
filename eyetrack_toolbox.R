#===============================================================================
# Blink detection function
# OUTPUT:
# This function returns a dataframe indicating:
# 1. start: index of the onset
# 2. end : index of the offset
# 3. mindur : duration of the blink (in samples) without the window
# 4. dur : duration of the blink (in samples) with the window
# 5. t_mindur : minimum duration of the saccade (in ms) without the window
# 6. t_dur : duration of the saccade(in ms) with the window
# INPUTS:
# For that, user must provide the following inputs:
# 1. sample : dataframe with x-y etr series
# 2. SAMPLING : sampling frequency at which data were collected
# 3. blink.threshold : duration of a blink in ms
# 4. blink.window : window to remove around the blink in ms
#===============================================================================

blink.detect <- function(sample,SAMPLING=1000,blink.threshold=70,blink.window=70) {
  secsamp <- 1000/SAMPLING
  sample$idx <- 1:nrow(sample)
  
  trl <- sample[,c("x","y")] # get x-y position
  
  blinkidx <- sample[which(rowSums(is.na(trl))>0),]
  idx_blink_sep <- which(diff(blinkidx$idx) != 1)
  
  df_blinks <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, 
                          c("start", "end", "mindur","dur","t_mindur","t_dur"))))
  if (empty(blinkidx)){
    nblinks <- 0
    df_blinks <- NULL
  } else{
    nblinks <- length(idx_blink_sep) + 1
    
    # the first onset
    onset <- blinkidx[1,'idx']
    if (nblinks==1) {
      offset <- tail(blinkidx$idx,1)
      minduration <- (offset - onset)*secsamp
      duration <- minduration + blink.window*2
      df_blinks[1,] <- c(onset - blink.window/secsamp,
                         offset + blink.window/secsamp,
                         minduration/secsamp,
                         duration/secsamp, 
                         minduration,
                         duration)
    } else{
      offset <- blinkidx[idx_blink_sep[1],'idx']
      minduration <- (offset - onset)*secsamp
      duration <- minduration + blink.window*2
      df_blinks[1,] <- c(onset - blink.window/secsamp,
                         offset + blink.window/secsamp,
                         minduration/secsamp,
                         duration/secsamp,
                         minduration,
                         duration)
      for(nb in 2:nblinks) {
        onset <- blinkidx[idx_blink_sep[nb-1]+1,'idx']
        if (nb == nblinks){
          offset <- tail(blinkidx$idx,1)
        } else {
          offset <- blinkidx[idx_blink_sep[nb],'idx']
        }
        minduration <- (offset - onset)*secsamp
        duration <- minduration + blink.window*2
        df_blinks[nb,] <- c(onset - blink.window/secsamp,
                            offset + blink.window/secsamp,
                            minduration/secsamp,
                            duration/secsamp,
                            minduration,
                            duration)
      }
    }
    
    df_blinks <- df_blinks[which(df_blinks$t_mindur > blink.threshold),]
    if (dim(df_blinks)[1] == 0) {
      df_blinks <- NULL
    } else{
      df_blinks <- cbind(label=unique(sample$label),df_blinks)
      df_blinks <- cbind(trial=unique(sample$trial),df_blinks)
      # if the start index is negative, means that the blink window covers
      # beyond of the signal beginning, therefore I have to reset the onset to 1
      # and then correct the durations  
      idxcorrect <- df_blinks$start < 1
      df_blinks[idxcorrect,'start'] <- 1
      df_blinks[idxcorrect,'dur'] <- df_blinks[idxcorrect,'end'] - df_blinks[idxcorrect,'start']
      df_blinks[idxcorrect,'t_dur'] <- df_blinks[idxcorrect,'dur']*secsamp
      
      # If the end index goes beyond the length of the sample means that we have to short
      # the blink duration
      idxcorrect <- df_blinks$end > nrow(sample)
      df_blinks[idxcorrect,'end'] <- nrow(sample)
      df_blinks[idxcorrect,'dur'] <- df_blinks[idxcorrect,'end'] - df_blinks[idxcorrect,'start']
      df_blinks[idxcorrect,'t_dur'] <- df_blinks[idxcorrect,'dur']*secsamp
    } 
  }
  return(df_blinks)
}

#===============================================================================
# aggregate.fixations
# Function adapted from tmalsburg/saccades repository
# OUTPUT:
# 1. A dataframe containing what is not a saccade neither a fixation
# INPUTS:
# 1. samples: a dataframe containing these columns:
#   - "x": x position
#   - "y": y position 
#   - "trial": number of the trial
#   - "label": a label for the trial, usually "RECORDING_SESSION_LABEL" 
#   - "saccblink": a column indicating with 1s saccades and blinks events
# 2. SAMPLING: sampling frequency at which data were collected
#===============================================================================
aggregate.fixations <- function(samples,SAMPLING) {
  
  # In saccade.events a 1 marks the start of a saccade and a -1 the
  # start of a fixation.
  
  saccblink.events <- sign(c(0, diff(samples$saccblink)))
  
  trial.numeric  <- as.integer(factor(samples$trial))
  trial.events   <- sign(c(0, diff(trial.numeric)))
  
  # New fixations start either when a saccade / blink ends or when a trial
  # ends:
  samples$fixation.id <- cumsum(saccblink.events==-1|trial.events==1)
  samples$idx <- 1:nrow(samples)
  samples$idx <- ifelse(trial.events==1, NA, samples$idx)
  
  # Discard samples that occurred during saccades and blinks:
  samples <- samples[!samples$saccblink,,drop=FALSE]
  
  fixations <- with(samples, data.frame(
    trial   = tapply(trial, fixation.id, function(x) x[1]),
    label   = tapply(label, fixation.id, function(x) x[1]),
    start   = tapply(idx,  fixation.id, min),
    end     = tapply(idx,    fixation.id, function(x) max(x, na.rm=TRUE)),
    x       = tapply(x,     fixation.id, stats::median),
    y       = tapply(y,     fixation.id, stats::median),
    mad.x   = tapply(x,     fixation.id, stats::mad),
    mad.y   = tapply(y,     fixation.id, stats::mad),
    # peak.vx = tapply(vx,    fixation.id, function(x) x[which.max(abs(x))]),
    # peak.vy = tapply(vy,    fixation.id, function(x) x[which.max(abs(x))]),
    stringsAsFactors=FALSE))
  
  fixations$dur <- fixations$end - fixations$start
  fixations$t_dur <- fixations$dur*(1000/SAMPLING)
  
  fixations
}

#===============================================================================
# label.blink.artifacts
# Function adapted from tmalsburg/saccades repository
# OUTPUT:
# 1. It appends a column named "event" indicating to each event whether it is a
#     a fixation or a blink.
# INPUT:
# 2. A dataframe containing what is not a saccade neither a fixation (e.g.: the
#     dataframe returned by the aggregte.fixation function)
#===============================================================================
label.blinks.artifacts <- function(fixations) {
  
  # Blink and artifact detection based on dispersion:
  lsdx <- log10(fixations$mad.x)
  lsdy <- log10(fixations$mad.y)
  median.lsdx <- stats::median(lsdx, na.rm=TRUE)
  median.lsdy <- stats::median(lsdy, na.rm=TRUE)
  mad.lsdx <- stats::mad(lsdx, na.rm=TRUE)
  mad.lsdy <- stats::mad(lsdy, na.rm=TRUE)
  
  # Dispersion too low -> blink:
  threshold.lsdx <- median.lsdx - 4 * mad.lsdx
  threshold.lsdy <- median.lsdy - 4 * mad.lsdy
  event <- ifelse((!is.na(lsdx) & lsdx < threshold.lsdx) &
                    (!is.na(lsdy) & lsdy < threshold.lsdy),
                  "blink", "fixation")
  
  # Dispersion too high -> artifact:
  threshold.lsdx <- median.lsdx + 4 * mad.lsdx
  threshold.lsdy <- median.lsdy + 4 * mad.lsdy
  event <- ifelse((!is.na(lsdx) & lsdx > threshold.lsdx) &
                    (!is.na(lsdy) & lsdy > threshold.lsdy),
                  "too dispersed", event)
  
  # Artifact detection based on duration:
  dur <- 1/fixations$dur
  median.dur <- stats::median(dur, na.rm=TRUE)
  mad.dur <- stats::mad(dur, na.rm=TRUE)
  
  # Duration too short -> artifact:
  threshold.dur <- median.dur + mad.dur * 5
  event <- ifelse(event!="blink" & dur > threshold.dur, "too short", event)
  
  factor(event, levels=c("fixation", "blink", "too dispersed", "too short"))
}

#===============================================================================
# Function vecvel() -- Microsaccade Toolbox 0.9
# (R-language Version)
# Authors: Ralf Engbert, Petra Sinn, Konstantin Mergenthaler, 
# and Hans Trukenbrod
# Date: February 20th, 2014
#===============================================================================
#-------------------------------------------------------------------------------
# Compute velocity times series from position data
#-------------------------------------------------------------------------------
vecvel <- function(x,SAMPLING=500,TYPE=2) {
  d <- dim(x)
  N <- d[1]
  v <- matrix(rep(0,2*N),ncol=2)
  
  if ( TYPE==2 ) {
    v[3:(N-2),] <- SAMPLING/6*(x[5:N,] + x[4:(N-1),] - x[2:(N-3),] - x[1:(N-4),])
    v[2,] = SAMPLING/2*(x[3,] - x[1,])
    v[(N-1),] = SAMPLING/2*(x[N,] - x[(N-2),])   
  }  else  {
    v[2:(N-1),] <- SAMPLING/2*(x[3:N,] - x[1:(N-2),])
  }
  return(v)
}

#===============================================================================
# microsacc
# Function adapted from Engbert 2015
# OUTPUT:
#   1. Dataframe containing in each row a saccade detected
# INPUTS:
#   1. sample: dataframe containing x and y columns
#   2. VFAC : velocity threshold factor for saccade detection
#   3. MINDUR : minimum duration of saccade in ms
#   4. SAMPLING: sampling frequency at which data were collected
#===============================================================================
microsacc <- function(sample,VFAC=5,MINDUR=3,SAMPLING=500) {
  
  # get NA indx
  # indxna <- rowSums(is.na(sample[,c('x','y')])) > 0
  # omit the NA values  
  # sample <- na.omit(sample)
  x <- cbind(sample$x,sample$y)
  # Compute velocity (this include NA)
  vna <- vecvel(x,SAMPLING=SAMPLING, TYPE=2)
  
  # this is the hampel filter to delete peak-outliers  
  # v[,1]  <- hampel(v[,1],10,t0=7)$y #check this line
  # v[,2]  <- hampel(v[,2],10,t0=7)$y #check this line
  
  # remove NA from velocity vector  
  v <- vna[!rowSums(is.na(vna)) > 0,] #check this line
  
  # Compute threshold
  medx <- median(v[,1])
  msdx <- sqrt( median((v[,1]-medx)^2) )
  medy <- median(v[,2]) 
  msdy <- sqrt( median((v[,2]-medy)^2) )
  
  if (msdx<1e-10 ) {
    msdx <- sqrt( mean(v[,1]^2) - (mean(v[,1]))^2 )
    if ( msdx<1e-10 ) {
      stop("msdx<realmin in microsacc.R")
    }
  }  
  if ( msdy<1e-10 ) {
    msdy <- sqrt( mean(v[,2]^2) - (mean(v[,2]))^2 )
    if ( msdy<1e-10 ) {
      stop("msdy<realmin in microsacc.R")
    }
  }  
  radiusx <- VFAC*msdx
  radiusy <- VFAC*msdy
  radius <- c(radiusx,radiusy)
  
  # Apply test criterion: elliptic treshold
  test <- (vna[,1]/radiusx)^2 + (vna[,2]/radiusy)^2
  
  indx <- which(test>1)
  # indx <- sample[which(test>1),'time'] 
  
  # Determine saccades
  N <- length(indx) 
  nsac <- 0
  sac <- NULL
  dur <- 1
  a <- 1
  k <- 1
  
  secsamp <- 1000/SAMPLING
  
  # Loop over saccade candidates
  while ( k<N ) {
    if ( indx[k+1]-indx[k]==1 ) {
      dur <- dur + 1
    } else {
      # Minimum duration criterion (exception: last saccade)
      if ( dur*secsamp>=MINDUR ) {
        nsac <- nsac + 1
        b <- k
        sac <- rbind(sac,c(indx[a],indx[b],dur,dur*secsamp,rep(0,6)))
      }
      a <- k+1
      dur <- 1
    }
    k <- k + 1
  }
  
  # Check minimum duration for last microsaccade
  if  ( dur*secsamp>=MINDUR ) {
    nsac <- nsac + 1
    b <- k
    sac <- rbind(sac,c(indx[a],indx[b],dur,dur*secsamp,rep(0,6)))
  }
  
  if ( nsac>0 ) {
    # Compute peak velocity, horizontal and vertical components
    for ( s in 1:nsac ) {
      # Onset and offset for saccades
      a <- sac[s,1] 
      b <- sac[s,2]
      idx <- a:b
      # Saccade peak velocity (vpeak) / 
      # magnitude of the vlocity vecotr |v| = v = sqrt(Vx^2 + Vy^2)
      vpeak <- max( sqrt( vna[idx,1]^2 + vna[idx,2]^2 ) )
      sac[s,5] <- vpeak
      # Saccade vector (dx,dy)
      dx <- x[b,1]-x[a,1] 
      dy <- x[b,2]-x[a,2] 
      sac[s,6:7] <- c(dx,dy)
      # Saccade amplitude (dX,dY)
      minx <- min(x[idx,1])
      maxx <- max(x[idx,1])
      miny <- min(x[idx,2])
      maxy <- max(x[idx,2])
      ix1 <- which.min(x[idx,1])
      ix2 <- which.max(x[idx,1])
      iy1 <- which.min(x[idx,2])
      iy2 <- which.max(x[idx,2])
      dX <- sign(ix2-ix1)*(maxx-minx)
      dY <- sign(iy2-iy1)*(maxy-miny)
      sac[s,8:9] = c(dX,dY)
      sac[s,10]  <- sqrt(dX^2 + dY^2)
    }
    colnames(sac) <- c('start','end','dur','t_dur','peakv',
                       'horizontalcomp', 'verticalcomp',
                       'horizontalamp','verticalamp','amp')
    # sac <- list(table=as.data.frame(sac),radius=radius)
    sac <- as.data.frame(sac)
    sac$radiusx <- radiusx
    sac$radiusy <- radiusy
    sac <- cbind(label=unique(sample$label),sac)
    sac <- cbind(trial=unique(sample$trial),sac)
    # sac$trial <- unique(sample$trial)
  } else  {
    sac <- NULL
  }
  
  return(sac)
}

#===============================================================================
# saccfixblink
# OUTPUT:
# 1. A list containing:
#   - summary: all events that ocurred during the trial
#   - fixations: information about fixations (start, end, dispersion ...)
#   - saccades: information about saccades (start, end, duration, amplitude...)
#   - blink: information about blinks (start, end, duration, minimumduration...)
# INPUTS:
# 1. sample: a dataframe containing these columns:
#   - "x": x position
#   - "y": y position 
#   - "trial": number of the trial
#   - "label": a label for the trial, usually "RECORDING_SESSION_LABEL"
# 2. SAMPLING: sampling frequency at which data were collected
# 3. blink.threshold : duration of a blink in ms
# 4. blink.window : window to remove around the blink in ms
# 5. MINDUR: minimum duration threshold for saccade detection (10 ms)
# 6. VFAC: velocity threshold factor for saccade detection
#===============================================================================
saccfixblink <- function(sample,SAMPLING,blink.threshold,blink.window,MINDUR,VFAC){
  
  # sample$time <- seq(1,nrow(sample),1000/SAMPLING)
  sample$saccblink <- 0 # column to fill with 1s when a saccade or blink takes place 
  
  # detect saccades
  # ------------------------------------------------------------------
  saccades <- microsacc(sample,VFAC,MINDUR,SAMPLING)
  
  # detect blinks
  # ------------------------------------------------------------------
  blinks <- blink.detect(sample,SAMPLING,blink.threshold,blink.window)
  
  # ... if blinks detected ...
  # ------------------------------------------------------------------
  if(!is.null(dim(blinks))) {
    # 1. get blink series 
    del.blinks <- c()
    for (nb in 1:nrow(blinks)) {
      del.blinks <- c(del.blinks,blinks[nb,'start']:blinks[nb,'end'])
    }
    blinks$event <- "blink"
    
    # 2. mark with 1 blinks
    sample[del.blinks,'saccblink']  <- 1
    
    # 3. while saccades detected --> get saccades series present inside the blink range
    if(!is.null(saccades)){
      del.sac <- c()
      for (nb in 1:nrow(blinks)){
        blink <- blinks[nb,c('start','end')]
        blink.series <- blink$start:blink$end
        for (sac in 1:nrow(saccades)){
          if (saccades[sac,'start'] %in% blink.series && saccades[sac,'end'] %in% blink.series){
            del.sac <- c(del.sac,sac)
          # if the start is inside but the end is outside, remove the saccade  
          } else if(saccades[sac,'start'] %in% blink.series && !(saccades[sac,'end'] %in% blink.series)){
            del.sac <- c(del.sac,sac)
            # fill with 1 until the end of the saccade
            sample[blinks[nb,'end']:saccades[sac,'end'] ,'saccblink'] <- 1 
            # ... and take the end of the saccade as the end of the blink
            blinks[nb,'end'] <- saccades[sac,'end']  
            blinks[nb,'dur'] <- blinks[nb,'end'] - blinks[nb,'start']
            blinks[nb,'t_dur'] <- blinks[nb,'dur']*(1000/SAMPLING)
          # if the start is outside but the end is inside, remove the saccade   
          } else if(!(saccades[sac,'start'] %in% blink.series) && saccades[sac,'end'] %in% blink.series){
            del.sac <- c(del.sac,sac)
            # fill with 1 until the start of the blink
            sample[saccades[sac,'start']:blinks[nb,'start'],'saccblink'] <- 1 
            # ... and take the start of the saccade as the start of the blink
            blinks[nb,'start'] <- saccades[sac,'start']  
            blinks[nb,'dur'] <- blinks[nb,'end'] - blinks[nb,'start']
            blinks[nb,'t_dur'] <- blinks[nb,'dur']*(1000/SAMPLING)
          } else{
            # fill with 1s saccades that survive (i.e., are outside the blink)
            sample[saccades[sac,'start']:saccades[sac,'end'],'saccblink'] <- 1 
          } 
        }
      }
      # 4. remove saccades inside the blink range
      if(!is.null(del.sac)){
        saccades <- saccades[-del.sac,] 
      } 
      # saccades <- saccades[-del.sac,]
      # ... if all saccades removed, the dataframe will be NULL ...  
      if (dim(saccades)[1] == 0){
        saccades <- NULL
      } else{
        saccades$event <- "saccade"
      } 
    }
    # ... if blinks not detected but saccades ...
    # ------------------------------------------------------------------
  } else if(!is.null(dim(saccades))){
    for (sac in 1:nrow(saccades)){
      # fill with 1s saccades 
      sample[saccades[sac,'start']:saccades[sac,'end'],'saccblink'] <- 1 
    }
    saccades$event <- "saccade"
  } 
  
  # get fixations after blinks and saccades were detected
  # ------------------------------------------------------------------
  fixations <- aggregate.fixations(sample, SAMPLING=SAMPLING)
  if(dim(fixations)[1] != 0){
    fixations$event <- label.blinks.artifacts(fixations)
    rownames(fixations) <- NULL
  } else{
    fixations <- NULL
  } 
  
  labels <- c('trial','label','start','end','dur','event')
  summary <- rbind(fixations[,labels], saccades[,labels], blinks[,labels])
  
  # order by start-event  
  summary <- summary[order(summary$start),]
  rownames(summary) <- NULL
  # summary <- cbind(trial=unique(sample$trial),summary)
  # summary$trial <- unique(sample$trial)
  
  # return(result)
  return(list(summary=summary,fixations=fixations,
              saccades=saccades,blinks=blinks))
} 

#===============================================================================
# plt.etr.series
# OUTPUT:
# 1. A plot showing x-y time series and events detected for a trial in shaded
#    colors:
#   - gray: fixations
#   - blue: saccades
#   - yellow: blinks
#   - red: artifacts
# INPUTS:
# 1. sample: a dataframe containing these columns:
#   - "x": x position
#   - "y": y position
# 2. SAMPLING: sampling frequency at which data were collected
# 3. summary: all events that ocurred during the trial
# 4. fixations: information about fixations (start, end, dispersion ...)
# 5. saccades: information about saccades (start, end, duration, amplitude...)
# 6. blink: information about blinks (start, end, duration, minimumduration...)
#===============================================================================
plt.etr.series <- function(sample,SAMPLING,summary,fixations,blinks,saccades){
  secsamp <- (1000/SAMPLING)
  sample$time <- c(0,seq(1,nrow(sample)-1)*secsamp)
  # sample$idx <- seq(1,nrow(sample))
  sample.melted <- reshape2::melt(sample, id = colnames(sample)[c(-1,-2)])
  plt <- ggplot(sample.melted, aes(x = time, y = value, color = variable)) +
    geom_point(size=0.1)
  
  # plot start-end of each event (they can be saccades, blinks, fixations, 
  # artifacts...) with vertical lines:  
  if(!missing(summary) && !is.null(summary)){
    plt <- plt +
      geom_vline(data=summary,aes(xintercept=start*secsamp),
                 size=0.1,color="black",aes(linetype=2))+
      geom_vline(data=summary,aes(xintercept=end*secsamp),
                 size=0.1,color="black",aes(linetype=2))
  }
  
  # mark x-y average position for each fixation and shade them
  # also shade artifacts in red:
  if(!missing(fixations) && !is.null(fixations)){
    artifacts <- subset(fixations, event == "too short" | event == "too dispersed")
    fixations <- subset(fixations, event == "fixation")
    
    if(dim(fixations)[1] !=0 ){
      plt <- plt +
        geom_segment(data=fixations,inherit.aes = FALSE, 
                     aes(x=start*secsamp,xend=end*secsamp,y=y,yend=y)) +
        geom_segment(data=fixations,inherit.aes = FALSE, 
                     aes(x=start*secsamp,xend=end*secsamp,y=x,yend=x)) +
        geom_rect(data=fixations,inherit.aes = FALSE, 
                  aes(xmin=start*secsamp,xmax=end*secsamp,ymin=-Inf,ymax=Inf),
                  alpha=0.1,fill="black")
      
    } 
    
    if(dim(artifacts)[1] !=0 ){
      plt <- plt +
        geom_rect(data=artifacts,inherit.aes = FALSE,aes(xmin=start*secsamp, 
                                                         xmax=end*secsamp, 
                                                         ymin=-Inf, ymax=Inf),
                  alpha=0.5,fill="red")
    } 
  } 
  
  # shade in yellow blinks  
  if (!missing(blinks) && !is.null(blinks)){
    if(dim(blinks)[1] !=0 ){
      plt <- plt +
        geom_rect(data=blinks,inherit.aes = FALSE,
                  aes(xmin=start*secsamp, xmax=end*secsamp, ymin=-Inf, ymax=Inf), 
                  alpha=0.1,fill="yellow")
      
    }
  } 
  
  # shade in blue saccades    
  if (!missing(saccades) && !is.null(saccades)){
    if(dim(saccades)[1] !=0 ){
      plt <- plt +
        geom_rect(data=saccades,inherit.aes = FALSE, 
                  aes(xmin=start*secsamp, xmax=end*secsamp, ymin=-Inf, ymax=Inf),
                  alpha=0.5,fill="blue")
      
    } 
  }
  
  plt <- plt + theme(panel.background = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     legend.key = element_rect(fill = "white"),
                     legend.title = element_text(size=20),
                     legend.text = element_text(size=18),
                     plot.title = element_text(hjust = 0.5),
                     text = element_text(size=20),
                     axis.text=element_text(size=18),
                     axis.title=element_text(size=18),
                     strip.text = element_text(size=20)) +
    scale_color_discrete(name='Position',labels= c('x','y')) +
    scale_y_continuous(limits=c(-6,6)) +
    scale_x_continuous(breaks = round(seq(min(sample$time), max(sample$time), 
                                          by = 400*secsamp),1)) +
    guides(color=guide_legend(override.aes = list(size=2))) +
    xlab("TIME [MS]") + ylab("AMPLITUDE [DEG]") 
  return(plt)
}  

#===============================================================================
# run.participant
# OUTPUT:
# 1. A list containing saccades, blinks and fixations for all the participant's
#   dataset
# INPUTS:
# 1. stimulidata: whole ETR dataset for a single participant containing these
#    columns typically returned by the edf converter:
#       - 'rxdeg'
#       - 'rydeg'
#       - 'lxdeg'
#       - 'lydeg'
#       - 'TRIAL_INDEX'
#       - 'RECORDING_SESSION_LABEL'
# 2. eye: aye to analyze, 'right' or 'left'
# 3. SAMPLING: sampling frequency at which data were collected
# 4. blink.threshold : duration of a blink in ms
# 5. blink.window : window to remove around the blink in ms
# 6. MINDUR: minimum duration threshold for saccade detection (10 ms)
# 7. VFAC: velocity threshold factor for saccade detection
#===============================================================================
run.participant <- function(stimulidata,eye,SAMPLING,blink.threshold,
                            blink.window,MINDUR,VFAC){
  if(missing(eye)){
    stop('Please, type right or left in the eye variable to start the analysis')
  } else if(eye=='right'){
    xy <- c('rxdeg','rydeg')
  } else if(eye=='left'){
    xy <- c('lxdeg','lydeg')
  } else{
    stop('Eye must be right or left')
  } 
  
  sample <- stimulidata[,c(xy,'TRIAL_INDEX', 'RECORDING_SESSION_LABEL')]
  names(sample) <- c("x", "y", "trial", "label") 
  
  summary <- c()
  fixations <- c()
  blinks <- c()
  saccades <- c()
  
  trials <- unique(sample$trial)
  for (trl in 1:length(trials)){
    trial <- subset(sample, trial == trials[trl])
    
    # check if non-null values exceed more than the 75 % of its length, if so,
    # it passes to the next trial  
    if(sum(rowSums(is.na(trial[,c('x','y')])) == 2) > nrow(trial)*0.75){
      result <- NULL
    } else {
      result <- saccfixblink(trial,SAMPLING,blink.threshold,blink.window,MINDUR,VFAC)
    } 
    
    summary <- rbind(summary,result$summary)
    fixations <- rbind(fixations,result$fixations)
    blinks <- rbind(blinks,result$blinks)
    saccades <- rbind(saccades,result$saccades)
  }
  
  if (is.null(summary)){
    names <- c("trial","label","start","end","dur","event")
    summary <- data.frame(matrix(ncol=length(names),nrow=1))
    colnames(summary) <- names
  } 
  
  if (is.null(fixations)){
    names <- c("trial","label","start","end","x","y","mad.x","mad.y","dur",
               "t_dur","event")
    fixations <- data.frame(matrix(ncol=length(names),nrow=1))
    colnames(fixations) <- names
  } 
  
  if (is.null(blinks)){
    names <- c("trial","label","start", "end", "mindur","dur","t_mindur","t_dur",
               "event")
    blinks <- data.frame(matrix(ncol=length(names),nrow=1))
    # ,dimnames=list(NULL, 
    #                                                             c("start", "end", "mindur","dur","t_mindur","t_dur"))))
    colnames(blinks) <- names 
  }
  
  if (is.null(saccades)){
    names <- c("trial","label",'start','end','dur','t_dur','peakv',
               'horizontalcomp', 'verticalcomp','horizontalamp','verticalamp',
               'amp',"radiusx","radiusy","event")
    saccades <- data.frame(matrix(ncol=length(names),nrow=1))
    colnames(saccades) <- names
  }
  
  result <- list(summary=summary, fixations=fixations,
                 blinks=blinks,saccades=saccades)
  
  return(mapply(cbind, result, "eye"=eye, SIMPLIFY=FALSE))
}

#===============================================================================
# plt.sbj.analysis
# This function plots a random trial for a subject
# OUTPUT: for a random trial given the subID
# 1. A plot showing x-y time series and detected events in shaded colors
# 2. A plot showing x-y time series
# INPUTS:
# 1. sbj: subject ID to plot a random trial
# 2. vfac: velocity threshold factor to plot
# 3. SAMPLING: sampling frequency at which data were collected
# 4. etr.db: whole EyeTracker database containing these columns:
#     - "RECORDING_SESSION_LABEL"
#     - "TRIAL_INDEX"
#     - "lxdeg"
#     - "lydeg"
#     - "rxdeg"
#     - "rydeg"
#     - "subID" <- subject ID
# 5. EYE: aye to analyze, 'right' or 'left'
# 6. result.db: a list containing saccades, blinks and fixations for all the 
#    participant's dataset
#===============================================================================
plt.sbj.analysis <- function(sbj,vfac,SAMPLING,etr.db,trl.idx,EYE,result.db,Run=NULL){
  # select a random subject
  if(missing(sbj) || is.null(sbj) || is.na(sbj)){
    stop('Please, type the subject ID to plot a trial')
  } 
  
  sbj.etr.db <- subset(etr.db, subID==sbj)
  trials <- unique(sbj.etr.db$TRIAL_INDEX)
  
  if(missing(trl.idx) || is.null(trl.idx) || is.na(trl.idx)){
    trl.idx <- sample(1:length(trials),1)
  } 
  
  sbj.trial <- subset(sbj.etr.db,TRIAL_INDEX==trl.idx)
  
  if(missing(EYE) || is.null(EYE) || is.na(EYE)){
    stop('Please, type right or left in the EYE variable to plot the analysis')
  } else if(EYE=='right'){
    xy <- c('rxdeg','rydeg')
  } else if(EYE=='left'){
    xy <- c('lxdeg','lydeg')
  } else{
    stop('Eye must be right or left')
  } 
  
  sbj.blinks <- subset(result.db$blinks,subID==sbj & trial==trl.idx & eye==EYE)
  sbj.summary <- subset(result.db$summary,subID==sbj & trial==trl.idx & eye==EYE)
  sbj.saccades <- subset(result.db$saccades,subID==sbj & trial==trl.idx & eye==EYE)
  sbj.fixations <- subset(result.db$fixations,subID==sbj & trial==trl.idx & eye==EYE)
  
  if(!missing(vfac) || !is.null(vfac)){
    sbj.blinks <- subset(sbj.blinks,VFAC==vfac)
    sbj.summary <- subset(sbj.summary,VFAC==vfac)
    sbj.saccades <- subset(sbj.saccades,VFAC==vfac)
    sbj.fixations <- subset(sbj.fixations, VFAC==vfac)
  } 
  
  if(!missing(Run) || !is.null(Run)){
    sbj.blinks <- subset(sbj.blinks,run == Run)
    sbj.summary <- subset(sbj.summary,run == Run)
    sbj.saccades <- subset(sbj.saccades,run == Run)
    sbj.fixations <- subset(sbj.fixations, run == Run)
    
    sbj.trial <- subset(sbj.trial, run == Run)
  }
  
  sbj.trial <- sbj.trial[,c(xy,'RECORDING_SESSION_LABEL')]
  names(sbj.trial) <- c("x", "y", "label") 
  
  plt.trial <- plt.etr.series(sbj.trial, SAMPLING = SAMPLING,
                              fixations=sbj.fixations ,
                              blinks=sbj.blinks,
                              saccades=sbj.saccades,
                              summary=sbj.summary)
  
  plt.raw.trial <- plt.etr.series(sbj.trial,SAMPLING=SAMPLING)
  
  return(list(plt.analysis=plt.trial,
              plt.raw=plt.raw.trial,
              trial=trl.idx,
              label=unique(sbj.trial$label)))
}


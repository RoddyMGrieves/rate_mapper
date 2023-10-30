#rate_mapper map spike and position data in 2D                                

  rmap = rate_mapper(pos,spk) maps the positions in pos to a 2D dwellmap
  and the spikes in spk to a 2D spikemap, computes a 2D firing rate map

  rmap = rate_mapper(pos,spk,rmset) uses additional settings specified by 
  an rmset structure (see below)

  rmap = rate_mapper(pos,spk,[],speedlift) uses a speedlift input to
  decrease computation time (see note 1 below)

  [rmap,dmap,smap,rmset,speedlift] = rate_mapper(pos,spk,[],speedlift)
  also returns the dwell time map in 'dmap', the spike map in 'smap',
  the settings used to generate the map (and some additional info, see
  below) and a speedlift output which can be passed to a later iteration
  of rate_mapper to decrease computation time (see note 1 below)

  main input options include:

  'pos'           -   [Nx2] Numeric matrix, the position data x,y coordinates
                      Units are in mm

  'spk'           -   [Nx2] Numeric matrix, the spike data x,y coordinates
                      Units are in mm

  rmset optional fields include:

  'method'        -   String or character vector or numeric scalar that
                      specifies the mapping method to be used on position
                      and spike data. 

                      'histogram'
                      Spike and position data are binned seperately using
                      the bivariate histogram method (histcounts2). 
                      Smoothing is performed using imgaussfilt or nanconv
                      depending on the value of rmset.smethod

                      Default value is 'histogram'.

  'binsize'       -   Scalar, positive integer that specifies the pixel 
                      or bin size to use for mapping, units are in mm.

                      Default value is 2mm.

  'ssigma'        -   Scalar, positive integer that specifies the sigma 
                      or standard deviation of the smoothing kernel, units 
                      are in mm except if smethod is set to 4 (boxcar smoothing)
                      in which case the units are in bins.

                      Default value is 4mm.

  'ssize'         -   Scalar, positive integer that specifies the size of 
                      the smoothing kernel, units are in mm.

                      Default is 2*ceil(2*(rmset.ssigma./rmset.binsize))+1

  'maplims'       -   1x4 vector [xmin ymin xmax ymax] specifies the desired 
                      outer boundaries of the rate map. Coordinates should
                      be in the pos and spk reference frame and in mm.

                      Default value is: [min(pos) max(pos)]

  'padding'       -   Scalar, positive integer that specifies the amount 
                      of space or padding to add around the position data 
                      when mapping. Units are in mm.

                      Default value is 0mm.

  'mindwell'      -   Positive scalar that specifies the duration of time 
                      a rat must spend in a bin for it to be considered 
                      visited. Unvisited bins are set to NaN. Units are in 
                      seconds.

                      Default value is 0s.

  'mindist'       -   Positive scalar that specifies the minimum distance a 
                      bin must be from some position data for it to be  
                      considered valid. Invalid bins are set to NaN. Units are 
                      in mm.

                      Default value is 40mm.

  'maxdist'       -   Positive scalar that specifies the maximum distance a 
                      bin can be from some position data before it is considered 
                      invalid. Invalid bins are set to NaN. Units are in mm.

                      Default value is 640mm.

  'srate'         -   Positive scalar that specifies the sampling rate of 
                      the position data. Units are in Hz.

                      Default value is 50Hz.

  'smethod'       -   Scalar, positive integer that specifies the smoothing
                      method to be used. 1 = smooth spikemap & dwellmap before
                      division using imgaussfilt, 2 = smooth after division using 
                      imfilter and nanconv, 3 = no smoothing, 4 = smooth with an
                      average boxcar of size 2*floor(rmset.ssigma/2)+1

                      Default value is 1 (smooth before)

  'ash'           -   Scalar, positive integer that specifies the number of
                      bin subdivisions that should be used for the average shifted
                      histogram or 'ash' method. Higher = greater smoothing but
                      increased computation time.

                      Default value is 6

  'kern'          -   String that specifies the kernel shape to use for the average
                      shifted histogram or 'ash' method. Options include 'biweight',
                      'epanechnikov', 'triangular', 'uniform' or 'box' and 'gaussian'
                      or 'normal'.

                      Default value is 'biweight'

  'steps'         -   Scalar, positive integer that specifies the number of
                      convolution filter sizes that should be used for the kernel
                      accelerated methods such as 'kadaptive'. Higher = greater
                      accuracy but increased computation time.

                      Default value is 32

  additional outputs included in rmset output:

  'map_pos'       -   [Nx2] Position data converted to match the firing rate map.
                      Each row specifies an x,y coordinate and matches the order
                      given by the pos input.

  'map_spk'       -   [Nx2] Spike data converted to match the firing rate map.
                      Each row specifies an x,y coordinate and matches the order
                      given by the spk input.

  'xgrid'         -   [1xN] Location of the firing rate map bin centers along the
                      x-axis (i.e. the columns of the rate map)

  'ygrid'         -   [1xN] Location of the firing rate map bin centers along the
                      y-axis (i.e. the columns of the rate map)

  'points_to_map' -   Anonymous function that can be used to convert x,y coordinates
                      to match the maps. Example:
                      mapXY = rmset.points_to_map(XY);
                      where XY is an Nx2 set of input coordinates in the same reference
                      frame as pos and spk, mapXY are the coordinates converted to
                      the firing rate map reference frame. See example below.

  Class Support
  -------------
  The input matrices pos and spk must be a real, non-sparse matrix of
  the following classes: uint8, int8, uint16, int16, uint32, int32,
  single or double.


  Notes
  -----
  1. If cells are recorded in the same session their spikemaps are 
     expected to differ while the position data and dwellmap will
     remain the same. To save time this function can accept a precomputed 
     output, skipping that computation (which is usually more time 
     consuming than the spikemap). Simply provide the 'speedlift' output 
     from rate_mapper to the next iteration of rate_mapper. For some 
     methods a dwellmap is not computed or would not speed up 
     computation if provided, so instead a matrix or similar form of 
     data output is used instead. Thus in some cases (i.e. 'histogram' 
     method) the speedlift matrix will be identical to the dwellmap, 
     but in other cases (i.e. 'kadaptive') it will take a different 
     form and contain different data.

  2. For the 'histogram' method and smoothing method 1, smoothing is 
     achieved using imgaussfilt. The FilterDomain is set to 'spatial'
     to ensure convolution in the spatial domain. 'Padding' is set to 
     a scalar value of 0 as there should be no position or spike data
     in bins outside the map limits.

  2. For the 'histogram' method and smoothing method 2, smoothing is 
     achieved using a Gaussian filter via nanconv with the 'nanout'
     option. This will smooth the data while ignoring nans.


  Example
  ---------
 % create dummy data and then generate a firing rate map
    np = 1000; % number of position points to simulate
    pos = rand(np,2)*1000; % simulate position data
    spk = normrnd(500,100,ceil(np/10),2); % simulate spike data

    [ratemap,dwellmap,spikemap,rmset] = rate_mapper(pos,spk); % generate map

    figure
    subplot(2,2,1)
    imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on; % plot ratemap
    daspect([1 1 1]); axis xy; title('Ratemap');

    subplot(2,2,2)
    imagesc(dwellmap,'alphadata',~isnan(dwellmap)); hold on; % plot dwellmap
    daspect([1 1 1]); axis xy; title('Dwellmap');

    subplot(2,2,3)
    imagesc(spikemap,'alphadata',~isnan(spikemap)); hold on; % plot spikemap
    daspect([1 1 1]); axis xy; title('Spikemap');

    subplot(2,2,4)
    imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on; % plot ratemap
    mapXY = rmset.points_to_map(pos); % convert position data to map coordinates
    plot(mapXY(:,1),mapXY(:,2),'w'); % plot map-scaled position data
    daspect([1 1 1]); axis xy; title('Ratemap + positions');

  See also imgaussfilt, nanconv

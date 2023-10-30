function [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset,speedlift)
%rate_mapper map spike and position data in 2D
%    __     ___ __           __  __  __ __ 
%   |__) /\  | |_   |\/| /\ |__)|__)|_ |__)
%   | \ /--\ | |__  |  |/--\|   |   |__| \                                      
%
%   rmap = rate_mapper(pos,spk) maps the positions in pos to a 2D dwellmap
%   and the spikes in spk to a 2D spikemap, computes a 2D firing rate map
%
%   rmap = rate_mapper(pos,spk,rmset) uses additional settings specified by 
%   an rmset structure (see below)
%
%   rmap = rate_mapper(pos,spk,[],speedlift) uses a speedlift input to
%   decrease computation time (see note 1 below)
%
%   [rmap,dmap,smap,rmset,speedlift] = rate_mapper(pos,spk,[],speedlift)
%   also returns the dwell time map in 'dmap', the spike map in 'smap',
%   the settings used to generate the map (and some additional info, see
%   below) and a speedlift output which can be passed to a later iteration
%   of rate_mapper to decrease computation time (see note 1 below)
%
%   main input options include:
%
%   'pos'           -   [Nx2] Numeric matrix, the position data x,y coordinates
%                       Units are in mm
% 
%   'spk'           -   [Nx2] Numeric matrix, the spike data x,y coordinates
%                       Units are in mm
%
%   rmset optional fields include:
%
%   'method'        -   String or character vector or numeric scalar that
%                       specifies the mapping method to be used on position
%                       and spike data. 
%
%                       'histogram'
%                       Spike and position data are binned seperately using
%                       the bivariate histogram method (histcounts2). 
%                       Smoothing is performed using imgaussfilt or nanconv
%                       depending on the value of rmset.smethod
%
%                       Default value is 'histogram'.
%
%   'binsize'       -   Scalar, positive integer that specifies the pixel 
%                       or bin size to use for mapping, units are in mm.
%
%                       Default value is 2mm.
%
%   'ssigma'        -   Scalar, positive integer that specifies the sigma 
%                       or standard deviation of the smoothing kernel, units 
%                       are in mm except if smethod is set to 4 (boxcar smoothing)
%                       in which case the units are in bins.
%
%                       Default value is 4mm.
%
%   'ssize'         -   Scalar, positive integer that specifies the size of 
%                       the smoothing kernel, units are in mm.
%
%                       Default is 2*ceil(2*(rmset.ssigma./rmset.binsize))+1
%
%   'maplims'       -   1x4 vector [xmin ymin xmax ymax] specifies the desired 
%                       outer boundaries of the rate map. Coordinates should
%                       be in the pos and spk reference frame and in mm.
%
%                       Default value is: [min(pos) max(pos)]
% 
%   'padding'       -   Scalar, positive integer that specifies the amount 
%                       of space or padding to add around the position data 
%                       when mapping. Units are in mm.
%
%                       Default value is 0mm.
% 
%   'mindwell'      -   Positive scalar that specifies the duration of time 
%                       a rat must spend in a bin for it to be considered 
%                       visited. Unvisited bins are set to NaN. Units are in 
%                       seconds.
%
%                       Default value is 0s.
% 
%   'mindist'       -   Positive scalar that specifies the minimum distance a 
%                       bin must be from some position data for it to be  
%                       considered valid. Invalid bins are set to NaN. Units are 
%                       in mm.
%
%                       Default value is 40mm.
% 
%   'maxdist'       -   Positive scalar that specifies the maximum distance a 
%                       bin can be from some position data before it is considered 
%                       invalid. Invalid bins are set to NaN. Units are in mm.
%
%                       Default value is 640mm.
% 
%   'srate'         -   Positive scalar that specifies the sampling rate of 
%                       the position data. Units are in Hz.
%
%                       Default value is 50Hz.
% 
%   'smethod'       -   Scalar, positive integer that specifies the smoothing
%                       method to be used. 1 = smooth spikemap & dwellmap before
%                       division using imgaussfilt, 2 = smooth after division using 
%                       imfilter and nanconv, 3 = no smoothing, 4 = smooth with an
%                       average boxcar of size 2*floor(rmset.ssigma/2)+1
%
%                       Default value is 1 (smooth before)
% 
%   'ash'           -   Scalar, positive integer that specifies the number of
%                       bin subdivisions that should be used for the average shifted
%                       histogram or 'ash' method. Higher = greater smoothing but
%                       increased computation time.
%
%                       Default value is 6
% 
%   'kern'          -   String that specifies the kernel shape to use for the average
%                       shifted histogram or 'ash' method. Options include 'biweight',
%                       'epanechnikov', 'triangular', 'uniform' or 'box' and 'gaussian'
%                       or 'normal'.
%
%                       Default value is 'biweight'
% 
%   'steps'         -   Scalar, positive integer that specifies the number of
%                       convolution filter sizes that should be used for the kernel
%                       accelerated methods such as 'kadaptive'. Higher = greater
%                       accuracy but increased computation time.
%
%                       Default value is 32
%
%   additional outputs included in rmset output:
%
%   'map_pos'       -   [Nx2] Position data converted to match the firing rate map.
%                       Each row specifies an x,y coordinate and matches the order
%                       given by the pos input.
%
%   'map_spk'       -   [Nx2] Spike data converted to match the firing rate map.
%                       Each row specifies an x,y coordinate and matches the order
%                       given by the spk input.
%
%   'xgrid'         -   [1xN] Location of the firing rate map bin centers along the
%                       x-axis (i.e. the columns of the rate map)
%
%   'ygrid'         -   [1xN] Location of the firing rate map bin centers along the
%                       y-axis (i.e. the columns of the rate map)
%
%   'points_to_map' -   Anonymous function that can be used to convert x,y coordinates
%                       to match the maps. Example:
%                       mapXY = rmset.points_to_map(XY);
%                       where XY is an Nx2 set of input coordinates in the same reference
%                       frame as pos and spk, mapXY are the coordinates converted to
%                       the firing rate map reference frame. See example below.
%
%   Class Support
%   -------------
%   The input matrices pos and spk must be a real, non-sparse matrix of
%   the following classes: uint8, int8, uint16, int16, uint32, int32,
%   single or double.
%
%
%   Notes
%   -----
%   1. If cells are recorded in the same session their spikemaps are 
%      expected to differ while the position data and dwellmap will
%      remain the same. To save time this function can accept a precomputed 
%      output, skipping that computation (which is usually more time 
%      consuming than the spikemap). Simply provide the 'speedlift' output 
%      from rate_mapper to the next iteration of rate_mapper. For some 
%      methods a dwellmap is not computed or would not speed up 
%      computation if provided, so instead a matrix or similar form of 
%      data output is used instead. Thus in some cases (i.e. 'histogram' 
%      method) the speedlift matrix will be identical to the dwellmap, 
%      but in other cases (i.e. 'kadaptive') it will take a different 
%      form and contain different data.
%
%   2. For the 'histogram' method and smoothing method 1, smoothing is 
%      achieved using imgaussfilt. The FilterDomain is set to 'spatial'
%      to ensure convolution in the spatial domain. 'Padding' is set to 
%      a scalar value of 0 as there should be no position or spike data
%      in bins outside the map limits.
% 
%   2. For the 'histogram' method and smoothing method 2, smoothing is 
%      achieved using a Gaussian filter via nanconv with the 'nanout'
%      option. This will smooth the data while ignoring nans.
%
%
%   Example
%   ---------
%  % create dummy data and then generate a firing rate map
%     np = 1000; % number of position points to simulate
%     pos = rand(np,2)*1000; % simulate position data
%     spk = normrnd(500,100,ceil(np/10),2); % simulate spike data
% 
%     [ratemap,dwellmap,spikemap,rmset] = rate_mapper(pos,spk); % generate map
% 
%     figure
%     subplot(2,2,1)
%     imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on; % plot ratemap
%     daspect([1 1 1]); axis xy; title('Ratemap');
% 
%     subplot(2,2,2)
%     imagesc(dwellmap,'alphadata',~isnan(dwellmap)); hold on; % plot dwellmap
%     daspect([1 1 1]); axis xy; title('Dwellmap');
% 
%     subplot(2,2,3)
%     imagesc(spikemap,'alphadata',~isnan(spikemap)); hold on; % plot spikemap
%     daspect([1 1 1]); axis xy; title('Spikemap');
% 
%     subplot(2,2,4)
%     imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on; % plot ratemap
%     mapXY = rmset.points_to_map(pos); % convert position data to map coordinates
%     plot(mapXY(:,1),mapXY(:,2),'w'); % plot map-scaled position data
%     daspect([1 1 1]); axis xy; title('Ratemap + positions');
% 
%   See also imgaussfilt, nanconv

% HISTORY:
% version 1.0.0, Release 01/06/19 Initial release
% version 1.1.0, Release 01/06/19 Histogram method borrowed from previous code
% version 1.1.1, Release 01/06/19 Common code moved to start
% version 1.2.0, Release 20/06/19 Leutgeb method added from previous code
% version 1.2.1, Release 22/06/19 Improvements to Leutgeb method
% version 1.3.0, Release 15/07/19 KSDE method working
% version 1.3.1, Release 17/07/19 Bug fix in KSDE
% version 1.4.0, Release 07/08/19 Yartsev adaptive method added from 3D code
% version 1.5.0, Release 08/09/19 Adaptive method added
% version 1.5.1, Release 12/09/19 Bsplines method added
% version 2.0.0, Release 13/05/21 Kernel accelerated adaptive method working
% version 2.1.0, Release 21/07/21 Kernel accelerated Yartsev adaptive method working
% version 3.0.0, Release 21/10/22 Bootstrap method added
% version 4.0.0, Release 14/11/22 Renamed from graphDATA to rate_mapper
% version 4.0.1, Release 16/11/22 Bug fix in inputs
% version 4.0.2, Release 17/11/22 Bug fix in grid vectors
% version 4.1.0, Release 17/11/22 Added considerable comments
% version 4.1.1, Release 17/11/22 Added fixed bin grid option to histogram
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Defaults
    defset              = struct;
    defset.method       = 'histogram';
    defset.binsize      = 20;
    defset.ssigma       = 80;
    defset.ash          = 16;
    defset.maplims      = [min(pos(:,1)) min(pos(:,2)) max(pos(:,1)) max(pos(:,2))];     
    defset.padding      = 0; % in mm
    defset.mindwell     = 0; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
    defset.mindist      = 40; % (mm, default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    defset.maxdist      = 640; % (mm, default 50) used by kadaptive, adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    defset.srate        = 50; % (default 50) sampling rate of data in Hz, used to calculate time
    defset.steps        = 32; % the number of convolution size steps to use for kadaptive    
    defset.kern         = 'biweight'; % kernel
    defset.smethod      = 1; % smoothing method, 1 = before division, 2 = after, 3 = no smoothing
    defset.bmethod      = 0; % binning method, 0 = use rmset.binsize and create a grid with binsize x binsize pixels, N = create a grid that is NxN and ends at the data limits
    defset.twindow      = 0.25; % fyhn method, time window over which to estimate instantaneous firing rate

    %% check that all parameters are included in rmset 
    % Fill in missing inputs in rmset using defset
    f1 = fieldnames(defset);
    if ~exist('rmset','var')
        rmset = struct;
    end
    f2 = fieldnames(rmset);
    for i = 1:size(f1,1)
        if ~ismember(f1{i},f2) % if this defset field does not exist in rmset
            rmset.(f1{i}) = defset.(f1{i}); % add it to rmset
        end
    end
    
    % Check to see if a speedlift input was provided, if not set to empty
    if ~exist('speedlift','var')
        speedlift = [];
    end    
    
    % check if kernel size was set or not
    % if not, scale according to smoothing sigma
    if ~isfield(rmset,'ssize') || isempty(rmset.ssize)
        rmset.ssize = 2*ceil(2*(rmset.ssigma./rmset.binsize))+1;
    end
        
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PREPARE INFO COMMON TO ALL METHODS
    % get spike and position data
    pox = double(pos(:,1));
    poy = double(pos(:,2)); 
    if ~isempty(spk)
        spx = double(spk(:,1));
        spy = double(spk(:,2)); 
    else
        spx = NaN;
        spy = NaN;
    end

    % centre data on the origin (i.e. centre on the mid point of the data convex hull)
    mid_x = mean(rmset.maplims([1 3])); % centre in X
    mid_y = mean(rmset.maplims([2 4])); % centre in Y
    pox = pox - mid_x; 
    poy = poy - mid_y;   
    spx = spx - mid_x; 
    spy = spy - mid_y;  
    mid_x = rmset.maplims([1 3]) - mean(rmset.maplims([1 3]));
    mid_y = rmset.maplims([2 4]) - mean(rmset.maplims([2 4]));
  
    % generate vectors for bin edges
    if rmset.bmethod > 0  
        % generate an NxN grid of bins where N = rmset.bmethod
        % there is no padding, maps will always be NxN
        xvec = linspace(-max(abs(mid_x)),max(abs(mid_x)),rmset.bmethod+1); 
        yvec = linspace(-max(abs(mid_y)),max(abs(mid_y)),rmset.bmethod+1); 
        xvec = unique(sort([-xvec xvec],'ascend')); % mirror these vectors and then sort them in ascending order
        yvec = unique(sort([-yvec yvec],'ascend'));     
        rmset.xgrid = xvec;
        rmset.ygrid = yvec;

        % add function for converting data to match map  
        fx = @(x) ( (x-[mean(rmset.maplims([1 3])) mean(rmset.maplims([2 4]))]) ./ [max(abs(mid_x)) max(abs(mid_y))] .* (rmset.bmethod./2) + ([length(xvec) length(yvec)]./2) );
        rmset.points_to_map = fx;     
    else
        % generate vectors for bin edges, these span from the middle to the extreme X and Y, then we mirror them across each axis
        % the +binsize is to ensure the vectors always encapsulate all the data
        % i.e. if the data span from 0:10 but the desired binsize is 6, the vector 0:6:10 is only [0 6] so any data between 6 and 10 would be lost.
        xvec = 0 : rmset.binsize : ( max(abs(mid_x)) + rmset.binsize + rmset.padding ); 
        yvec = 0 : rmset.binsize : ( max(abs(mid_y)) + rmset.binsize + rmset.padding );    
        xvec = unique(sort([-xvec xvec],'ascend')); % mirror these vectors and then sort them in ascending order
        yvec = unique(sort([-yvec yvec],'ascend'));
        rmset.xgrid = xvec;
        rmset.ygrid = yvec;

        % add function for converting data to match map  
        fx = @(x) ( ( ( x-[mean(rmset.maplims([1 3])) mean(rmset.maplims([2 4]))] ) ./ rmset.binsize ) + ([length(xvec) length(yvec)]./2) );
        rmset.points_to_map = fx;        
    end
    % project the position and spike data on to the map coordinates
    rmset.map_pos = fx(pos(:,[1 2]));
    if ~isempty(spk)
        rmset.map_spk = fx(spk(:,[1 2]));  
    else
        rmset.map_spk = [];
    end
        
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> IMPLEMENTED MAPPING METHODS
    switch rmset.method
        
%% >>>>>>>>>>>>>>>>>>>> Bivariate histogram                               
        case {'histogram'}
% REF
% O'Keefe (1983) Spatial memory within and without the hippocampal system
% Muller, Kubie, and Ranck (1987) Spatial firing patterns of hippocampal complex-spike cells in a fixed environment
% https://doi.org/10.1523/jneurosci.07-07-01935.1987

% DESCRIPTION IN REF
% The X and Y positions are used as indices into a 64 x 64 time-in-location array and a 64 x 64 spikes-in-location array. The indexed
% element in the time array is incremented by 1; the same element in the spike array i incremented by the number of spikes fired during the 1/60th
% sec sample interval. This sequence is repeated for all the samples in the time series. A 64 x 64 rate array is then filed by dividing the time
% array into the spike array on an element-by-element basis. Unvisited pixels are marked by setting the appropriate element of the rate array to -1.

% SIMPLE DESCRIPTION            
% A bivariate histogram is generated for the position data and the spike data seperately. Firing rate in each pixel is then calculated by dividing
% the spike histogram by the position histogram multiplied by the sampling interval of the position data so that the final product is in units of Hz.
% The original papers did not include any smoothing, but these histograms are generally smoothed nowadays before the division occurs. This smoothing 
% is usually Gaussian but may also be an average/uniform boxcar kernel. Histograms are extremely fast to compute and combined with smoothing provide
% suprisingly accurate maps.
            % generate spikemap and dwellmap if necessary
            % coordinates are backwards so the resulting map has the correct orientation in IJ
            spikemap = histcounts2(spy,spx,yvec,xvec);

            % the dwellmap only needs to be computed if one wasn't provided
            if ~isempty(speedlift) && ~all(isnan(speedlift(:)))
                dwellmap = speedlift; 
                expected_size = [length(yvec)-1,length(xvec)-1];                
                if ~all(size(dwellmap) == expected_size)
                    error(sprintf('Size of provided dwellmap (%s) does not match the expected size {%s)... exiting',mat2str(size(dwellmap)),mat2str(expected_size)));        
                end                
            else
                dwellmap = histcounts2(poy,pox,yvec,xvec) .* (1/rmset.srate);
                if rmset.smethod==1 && rmset.ssigma>0
                    dwellmap = imgaussfilt(dwellmap,rmset.ssigma./rmset.binsize,'FilterSize',rmset.ssize,'FilterDomain','spatial','Padding',0);
                elseif rmset.smethod==4 && rmset.ssigma>0 % boxcar smoothing
                    box_size = 2*floor(rmset.ssigma/2)+1; % make sure boxcar is an odd size
                    H = ones(box_size,box_size) ./ (box_size^2); % an average boxcar
                    dwellmap = imfilter(dwellmap,H,0,'same','conv');
                end
            end                
            dwellmap(dwellmap < rmset.mindwell | isinf(abs(dwellmap))) = NaN;            

            % generate and smooth ratemap
            if rmset.smethod==1 && rmset.ssigma>0 % if we want to smooth spikes and time then calculate ratio
                spikemap = imgaussfilt(spikemap,rmset.ssigma./rmset.binsize,'FilterSize',rmset.ssize,'FilterDomain','spatial','Padding',0);
                ratemap = spikemap ./ dwellmap;    
            elseif rmset.smethod==2 && rmset.ssigma>0 % if we want to calculate ratio then smooth result
                H = fspecial('gaussian',rmset.ssize,rmset.ssigma./rmset.binsize); % kernel for smoothing                
                ratemap = nanconv(spikemap./dwellmap,H,'nanout');                       
            elseif rmset.smethod==3 || rmset.ssigma==0 % if we want no smoothing at all 
                ratemap = spikemap ./ dwellmap;  
            elseif rmset.smethod==4 && rmset.ssigma>0 % boxcar smoothing
                box_size = 2*floor(rmset.ssigma/2)+1; % make sure boxcar is an odd size
                H = ones(box_size,box_size) ./ (box_size^2); % an average boxcar
                spikemap = imfilter(spikemap,H,0,'same','conv');  
                ratemap = spikemap ./ dwellmap;      
            end

            % Make sure low dwell times and spurious results are removed for consistency
            spikemap(isnan(dwellmap) | isinf(abs(spikemap))) = NaN;
            ratemap(isnan(dwellmap) | isinf(abs(ratemap))) = NaN;
            speedlift = dwellmap;
            
%% >>>>>>>>>>>>>>>>>>>>  Averaged Shifted Histogram (ASH)                 
        case {'ash'}        
% REF
% David Scott (1992) Multivariate Density Estimation, John Wiley, (chapter 5)
% https://onlinelibrary.wiley.com/doi/book/10.1002/9781118575574

% DESCRIPTION IN REF
% A simple device has been proposed for eliminating the bin edge problem of the frequency polygon while retaining many of the computational advantages of a density
% estimate based on bin counts. Scott (1983, 1985b) considered the problem of choosing among the collection of multivariate frequency polygons, each with the same
% smoothing parameter but differing bin origins. Rather than choosing the “smoothest” such curve or surface, he proposed averaging several of the shifted frequency polygons. 
% As the average of piecewise linear curves is also piecewise linear, the resulting curve appears to be a frequency polygon as well. If the weights are nonnegative and
% sum to 1, the resulting “averaged shifted frequency polygon” (ASFP) is nonnegative and integrates to 1. A nearly equivalent device is to average several shifted histograms, 
% which is just as general but simpler to describe and analyze. The result is the “averaged shifted histogram” (ASH). Since the average of piecewise constant functions such as the
% histogram is also piecewise constant, the ASH appears to be a histogram as well. In practice, the ASH is made continuous using either of the linear interpolation schemes
% described for the frequency polygon in Chapter 4, and will be referred to as the frequency polygon (FP) of the (ASH). The ASH is the practical choice for computationally 
% and statistically efficient density estimation. Algorithms for its evaluation are described in detail.

% SIMPLE DESCRIPTION            
% Histograms are very fast to compute, but kernel smoothed density estimates are more accurate, one tradeoff is to instead calculate a lot of different histograms with slightly 
% offset bins (i.e. the bin edges are moved with respect to the data by a fraction of 1 bin length). The result is a histogram with a much lower MISE than a single histogram.
% In practice the result should not need to be smoothed and an added benefit is that the resulting histograms still sum to the number of data points.
% In reality we don't make multiple histograms but instead make one very fine scale one and then use a special weighting kernel to average across them.
            rmset.ash = max([ceil(rmset.ash) 1]);
            h = rmset.binsize;
            m = rmset.ash;
            delta = h/m; % fine grid increment
            
            % prepare fine bin centres
            xf = min(xvec):delta:max(xvec);
            yf = min(yvec):delta:max(yvec);
            rmset.xgrid_ash = xf;
            rmset.ygrid_ash = yf;

            % create biweight (quartic) filter and apply it to both maps
            s = rmset.ssigma;
            switch rmset.kern
                case {'biweight'}
                    kern = @(x) (15/16)*(1-x.^2).^2;
                case {'epanechnikov'}
                    kern = @(x) (3/4)*(1-x.^2);
                case {'triangular'}
                    kern = @(x) (1-abs(x));                    
                case {'uniform','box'}
                    kern = @(x) ones(size(x)).*.5;
                case {'gaussian','normal'}
                    kern = @(x) 1/(s*sqrt(2*pi)) .* exp(-(x.^2./(2*s.^2)));     
            end
            w = zeros(2*m-1,2*m-1);
            w(ceil(size(w,1)/2),ceil(size(w,2)/2)) = 1; % set center pixel to 1
            d = bwdist(w,'euclidean'); % distance within kernel from center
            d = d ./ max([d(:); eps]); % normalized
            wm = kern(double(d)); % create kernel
            wm = wm ./ sum(wm(:),'omitnan');

            % generate spikemap and dwellmap if necessary using fine grid
            % coordinates are backwards so the resulting map has the correct orientation in IJ
            s1 = histcounts2(spy,spx,yf,xf);
            spikemap = imfilter(s1,wm,0,'same');
            
            % the dwellmap only needs to be computed if one wasn't provided
            if ~isempty(speedlift) && ~all(isnan(speedlift(:)))
                dwellmap = speedlift; 
                expected_size = [length(yf)-1,length(xf)-1];
                if ~all(size(dwellmap) == expected_size)
                    error(sprintf('Size of provided dwellmap (%s) does not match the expected size {%s)... exiting',mat2str(size(dwellmap)),mat2str(expected_size)));        
                end                
            else
                d1 = histcounts2(poy,pox,yf,xf) .* (1/rmset.srate);
                dwellmap = imfilter(d1,wm,0,'same');              
            end  
            dwellmap(dwellmap < rmset.mindwell | isinf(abs(dwellmap))) = NaN;                   
            speedlift = dwellmap;     

            % generate ratemap
            ratemap = spikemap ./ dwellmap;   

            % make sure there are no weird values in the ratemap
            spikemap(isnan(dwellmap) | isinf(abs(spikemap))) = NaN;
            ratemap(isnan(dwellmap) | isinf(abs(ratemap))) = NaN;
                 
%% >>>>>>>>>>>>>>>>>>>> Adaptive                             
        case {'adaptive','yadaptive'}
% REF
% Skaggs and McNaughton (1998) Spatial Firing Properties of Hippocampal CA1 Populations in an Environment Containing Two Visually Identical Regions
% https://doi.org/10.1523/JNEUROSCI.18-20-08455.1998
% Skaggs, McNaughton, Wilson and Barnes (1996) Theta phase precession in hippocampal neuronal populations and the compression of temporal sequences
% https://doi.org/10.1002/(SICI)1098-1063(1996)6:2%3C149::AID-HIPO6%3E3.0.CO;2-K

% DESCRIPTION IN REF
% The 64 x 64 pixel firing rate maps were constructed using an “adaptive smoothing” method that has been described in previous publications
% (Skaggs et al., 1996). Briefly, the method is designed to optimize the tradeoff between blurring error (attributable to averaging together 
% data from locations with different true firing rates) and sampling error (the statistical error attributable to the limited number of samples available).
% To calculate the firing rate at a given point, a circle centered on the point is gradually expanded until the following criterion is met:
% r >= a / (n * sqrt(s))
% where a is a constant, r is the radius of the circle in pixels, n is the number of 50-msec-long occupancy samples lying within the circle, and s is the
% total number of spikes contained in those occupancy samples. Once this criterion was met, the firing rate assigned to the point was equal to s/n.
% For this experiment, a was set to the value 1000.

            % prepare bin centres
            [X1,Y1] = meshgrid(movmean(xvec,2,'EndPoints','discard'),movmean(yvec,2,'EndPoints','discard')); 
            xcen = X1(:);
            ycen = Y1(:);           
            ratemap = NaN(size(X1));
            spikemap = NaN(size(X1));
            dwellmap = NaN(size(X1));
            radmap = NaN(size(X1));
            
            k = sort(unique([linspace(rmset.binsize,rmset.maxdist,rmset.steps) rmset.mindist]));
            mindx = find(k==rmset.mindist,1);
            
            % Run through every bin
            for bb = 1:numel(xcen) % for every bin  
                % logical index which we can use to cut position data to within a square with side length config.maxdist of the bin
                rindx = pox>xcen(bb)-rmset.maxdist & pox<xcen(bb)+rmset.maxdist & poy>ycen(bb)-rmset.maxdist & poy<ycen(bb)+rmset.maxdist;                  
                if ~any(rindx) % if there are no position data within this distance of the bin, the ratemap should stay NaN
                    continue
                else
                    dp = ( sum(([pox(rindx) poy(rindx)]-[xcen(bb) ycen(bb)]).^2,2) ); % calculate the squared distance to every position data point                    
                end                         

                % logical index which we can use to cut spike data to within a square with side length config.maxdist of the bin                
                rindx = spx>xcen(bb)-rmset.maxdist & spx<xcen(bb)+rmset.maxdist & spy>ycen(bb)-rmset.maxdist & spy<ycen(bb)+rmset.maxdist;  
                if ~any(rindx) % if there are no spike data within this distance of the bin, the ratemap should equal 0
                    ratemap(bb) = 0;
                    continue
                else        
                    ds = ( sum(([spx(rindx) spy(rindx)]-[xcen(bb) ycen(bb)]).^2,2) ); % calculate the squared distance to every spike data point 
                end

                % count the spikes and position data points falling within each unique distance
                n = sum(dp(:) <= k(:)'.^2,1); % count the number of position data falling within each distance                   
                s = sum(ds(:) <= k(:)'.^2,1); % count the number of spikes falling within each distance   
   
                if n(mindx)*(1/rmset.srate) < rmset.mindwell
                    ratemap(bb) = NaN;
                    continue
                end                 
                
                % >>>>>>>>>>>>>>>>>>>> Adaptive (Skaggs method)                                             
                if strcmp(rmset.method,'adaptive')
                    a = ones(size(k)) .* rmset.ssigma; % this is the "a" in the adaptive equation
                    r = k ./ rmset.binsize; % this is the "r" in the adaptive equation
                    rindx = find( r(:) >= a(:) ./ (n(:) .* sqrt(s(:))) , 1, 'first' ); % this is the full "r >= a / (n * sqrt(s))"

                    if isempty(rindx)
                        spikemap(bb) = s(end);
                        dwellmap(bb) = n(end)*(1/rmset.srate);                    
                        ratemap(bb) = spikemap(bb) ./ dwellmap(bb);
                        radmap(bb) = k(end);
                    else
                        spikemap(bb) = s(rindx);
                        dwellmap(bb) = n(rindx)*(1/rmset.srate);                      
                        ratemap(bb) = spikemap(bb) ./ dwellmap(bb);
                        radmap(bb) = k(rindx);
                    end
                    rmset.radmap = single(radmap); % map of bin radii (mm)
                    
                % >>>>>>>>>>>>>>>>>>>> Adaptive (Yartsev 1s method)                                              
                elseif strcmp(rmset.method,'yadaptive')
                    n = n .* (1/rmset.srate); % convert to time

                    % calculate the adaptive equation and find the minimum radius
                    rindx = find( n >= 1 , 1, 'first' ); % find the first radius containing more than 1 second of data
                    if isempty(rindx)
                        spikemap(bb) = s(end);
                        dwellmap(bb) = n(end);
                        ratemap(bb) = spikemap(bb) ./ dwellmap(bb);
                        radmap(bb) = k(end);
                    else
                        spikemap(bb) = s(rindx);
                        dwellmap(bb) = n(rindx);                    
                        ratemap(bb) = spikemap(bb) ./ dwellmap(bb);
                        radmap(bb) = k(rindx);
                    end  
                    rmset.radmap = single(radmap); % map of bin radii (mm)                    
                   
                end                
            end
            
%% >>>>>>>>>>>>>>>>>>>> Adaptive (kernel accelerated method)                                   
        case {'kadaptive','kyadaptive'}  
% SIMPLE DESCRIPTION            
% There is no Ref for this method as I don't think it has been proposed before. This mixes the benefits and methods underlying both the average shifted histogram
% and adaptive binning. We create a fine grid (suggested is 1cm or less binsize) and bin the spike and position data into this, like we do with ASH
% Then we apply a kernel to this like the ASH but instead of weighting this to average over histograms we sum all the spikes and position data within a series of distances
% of each bin, like with the adaptive method. Then we apply the same formula as in the adaptive method to every bin, finding the radius for each that meets the criteria
% Because we bin the data into a fine grid before this process we don't have to compute distances which saves a lot of time and the kernel we use is very simple
% and can be applied using convolution methods for an overall very fast process which should approximate a KDE

% REF for original adaptive method
% Skaggs and McNaughton (1998) Spatial Firing Properties of Hippocampal CA1 Populations in an Environment Containing Two Visually Identical Regions
% https://doi.org/10.1523/JNEUROSCI.18-20-08455.1998
% Skaggs, McNaughton, Wilson and Barnes (1996) Theta phase precession in hippocampal neuronal populations and the compression of temporal sequences
% https://doi.org/10.1002/(SICI)1098-1063(1996)6:2%3C149::AID-HIPO6%3E3.0.CO;2-K

% DESCRIPTION IN REF for original adaptive method
% The 64 x 64 pixel firing rate maps were constructed using an “adaptive smoothing” method that has been described in previous publications
% (Skaggs et al., 1996). Briefly, the method is designed to optimize the tradeoff between blurring error (attributable to averaging together 
% data from locations with different true firing rates) and sampling error (the statistical error attributable to the limited number of samples available).
% To calculate the firing rate at a given point, a circle centered on the point is gradually expanded until the following criterion is met:
% r >= a / (n * sqrt(s))
% where a is a constant, r is the radius of the circle in pixels, n is the number of 50-msec-long occupancy samples lying within the circle, and s is the
% total number of spikes contained in those occupancy samples. Once this criterion was met, the firing rate assigned to the point was equal to s/n.
% For this experiment, a was set to the value 1000.
            
            % generate spikemap and dwellmap if necessary
            max_bin_distance = rmset.maxdist ./ rmset.binsize;
            kvect = sort(unique([ceil(linspace(1,max_bin_distance,rmset.steps)) ceil(rmset.mindist/rmset.binsize)]));
            mindx = find(kvect==ceil(rmset.mindist/rmset.binsize),1);
            if ~isempty(speedlift) && ~all(isnan(speedlift(:)))
                dmaps = speedlift(:,:,2:end);
                dwellmap = speedlift(:,:,1);                
            else
                % coordinates are backwards so the resulting map has the correct orientation in IJ
                dwellmap = histcounts2(poy,pox,yvec,xvec);            
                dmaps = NaN([size(dwellmap) length(kvect)]);
                for kk = 1:length(kvect)
                    SE = strel('disk',kvect(kk));
                    kern = double(SE.Neighborhood);
                    % kern = double(fspecial('disk',kvect(kk))>0); % faster and possibly better                       
                    dmaps(:,:,kk) = imfilter(dwellmap,kern,0,'same');
                end  
                speedlift = cat(3,dwellmap,dmaps);
            end
            
            % coordinates are backwards so the resulting map has the correct orientation in IJ
            spikemap = histcounts2(spy,spx,yvec,xvec);
            smaps = NaN([size(spikemap) length(kvect)]);
            karea = NaN(size(kvect));
            for kk = 1:length(kvect)
                SE = strel('disk',kvect(kk));
                kern = double(SE.Neighborhood);
                % kern = double(fspecial('disk',kvect(kk))>0); % faster and possibly better                                                    
                karea(kk) = sum(kern(:)).^2;
                smaps(:,:,kk) = imfilter(spikemap,kern,0,'same');  
            end

            % >>>>>>>>>>>>>>>>>>>> Adaptive (Skaggs method)  
            if strcmp(rmset.method,'kadaptive')            
                % adaptive equation
                amaps = ones(size(dmaps)) .* rmset.ssigma; % this is the "a" in the adaptive equation
                rmaps = amaps ./ (dmaps .* sqrt(smaps)); % this is the "a / (n * sqrt(s))" in the adaptive equation
                distmap = ones(size(dmaps)) .* permute(kvect,[1 3 2]) + 0.5; % this is the "r" in the adaptive equation
                emaps = distmap>=rmaps; % this is the full "r >= a / (n * sqrt(s))"

                % fill in the ratemap values
                ratemap = NaN(size(dwellmap),'double');
                radmap = NaN(size(dwellmap),'double'); % map of bin radii (mm)  
                for pp = 1:numel(dwellmap) % for every bin in the dwellmap
                    [i,j] = ind2sub(size(dwellmap),pp); % get its row column index
                    if ~any(emaps(i,j,:)) % if none of the bins satisfied the adaptive equation, use the maximum radius
                        indx = size(emaps,3);                                        
                    else % otherwise, use the first radius that satisfied the equation
                        indx = find(emaps(i,j,:),1,'first');                    
                    end
                    if (dmaps(i,j,indx)*(1/rmset.srate)) < rmset.mindwell || isnan(dmaps(i,j,mindx)) % if the bin was unvisited according to min dwell duration
                        ratemap(i,j) = NaN; % set it to NaN regardless of its value                 
                    else
                        ratemap(i,j) = squeeze(smaps(i,j,indx) ./ (dmaps(i,j,indx).*(1/rmset.srate))); % otherwise, set the firing rate to be "s/n"
                        radmap(i,j) = (kvect(indx)+0.5)*rmset.binsize;
                    end
                end
                rmset.radmap = single(radmap); % map of bin radii (mm)  
                
            % >>>>>>>>>>>>>>>>>>>> Adaptive (Yartsev method)                                              
            elseif strcmp(rmset.method,'kyadaptive')
                tmaps = (dmaps .* (1/rmset.srate)) > rmset.ssigma;                                
                
                % fill in the ratemap values
                ratemap = NaN(size(dwellmap),'double');
                radmap = NaN(size(dwellmap),'double'); % map of bin radii (mm)  
                for pp = 1:numel(dwellmap) % for every bin in the dwellmap
                    [i,j] = ind2sub(size(dwellmap),pp); % get its row column index            
                    if ~any(tmaps(i,j,:))
                        indx = size(tmaps,3);                                        
                    else
                        indx = find(tmaps(i,j,:),1,'first');                    
                    end
                    if (dmaps(i,j,indx)*(1/rmset.srate)) < rmset.mindwell || isnan(dmaps(i,j,mindx)) % if the bin was unvisited according to min dwell duration
                        ratemap(i,j) = NaN; % set it to NaN regardless of its value                       
                    else
                        ratemap(i,j) = squeeze(smaps(i,j,indx) ./ (dmaps(i,j,indx))); % otherwise, set the firing rate to be spikes / time
                        radmap(i,j) = (kvect(indx)+0.5)*rmset.binsize;
                    end
                end  
                rmset.radmap = single(radmap); % map of bin radii (mm)  
              
%             % >>>>>>>>>>>>>>>>>>>> K Nearest Neighbors                                              
%             elseif strcmp(rmset.method,'knn')   
%                                 
%                 % fill in the ratemap values
%                 ratemap = NaN(size(dwellmap),'double');
%                 radmap = NaN(size(dwellmap),'double'); % map of bin radii (mm)  
%                 d = NaN(size(dwellmap),'double');
%                 s = NaN(size(dwellmap),'double');
%                 
%                 
%                 tmap1 = smaps > rmset.ssigma;   
%                 tmap2 = dmaps > rmset.ssigma;                   
%                 for pp = 1:numel(dwellmap) % for every bin in the dwellmap
%                     [i,j] = ind2sub(size(dwellmap),pp); % get its row column index 
% 
%                     % position data
%                     if ~any(tmap2(i,j,:))
%                         indx1 = size(tmap2,3);                                        
%                     else
%                         indx1 = find(tmap2(i,j,:),1,'first');                    
%                     end
%                     if (dmaps(i,j,indx1)*(1/rmset.srate)) < rmset.mindwell || isnan(dmaps(i,j,mindx)) % if the bin was unvisited according to min dwell duration
%                         v2 = NaN; % set it to NaN regardless of its value                       
%                     else
% %                         v2 = rmset.ssigma ./ ( (pi*(kvect(indx1)+0.5).^2).*numel(pox) );
%                         v2 = rmset.ssigma ./ ( karea(indx1).*numel(pox) );
%                     end
%                         
%                     % spikes
%                     if ~any(tmap1(i,j,:))
%                         indx2 = size(tmap1,3);                                        
%                     else
%                         indx2 = find(tmap1(i,j,:),1,'first');                    
%                     end
%                     if (dmaps(i,j,indx1)*(1/rmset.srate)) < rmset.mindwell || isnan(dmaps(i,j,mindx)) % if the bin was unvisited according to min dwell duration
%                         v1 = NaN; % set it to NaN regardless of its value                       
%                     else
%                         v1 = rmset.ssigma ./ ( karea(indx2).*numel(spx) );
%                     end                    
% 
%                     
%                     % rate
%                     ratemap(i,j) = v1 ./ v2;
%                     d(i,j) = v2;
%                     s(i,j) = v1;
%                     
% %                     keyboard
%                 end
% figure
% subplot(2,2,1)
% plot(pox,poy,'k'); hold on;
% plot(spx,spy,'r.','MarkerSize',20)
% daspect([1 1 1])
% axis xy
% 
% subplot(2,2,2)
% imagesc(s)
% daspect([1 1 1])
% colorbar
% axis xy
%                 
% subplot(2,2,3)
% imagesc(d)
% daspect([1 1 1])
% colorbar
% axis xy
% 
% subplot(2,2,4)
% imagesc(ratemap)
% daspect([1 1 1])
% colorbar
% axis xy
% 
% keyboard
            end
              
%% >>>>>>>>>>>>>>>>>>>> Kernel smoothed density estimate (KSDE)                       
        case {'ksde'}
% This method uses a multivariate kernel density estimate approach to estimate probability density functions for the position data 
% and spike data separately. The firing rate map is then calculated as the ratio of these two maps. Here smoothing is achieved using 
% the bandwidth of the KDE 
            
            % prepare bin centres
            [X1,Y1] = meshgrid(movmean(xvec,2,'EndPoints','discard'),movmean(yvec,2,'EndPoints','discard')); 
            xcen = X1(:);
            ycen = Y1(:);
                        
            % the dwellmap only needs to be computed if one wasn't provided
            if ~isempty(speedlift) && ~all(isnan(speedlift(:)))
                dwellmap = speedlift; 
                expected_size = [length(yvec)-1,length(xvec)-1];
                if ~all(size(dwellmap) == expected_size)
                    error(sprintf('Size of provided dwellmap (%s) does not match the expected size {%s)... exiting',mat2str(size(dwellmap)),mat2str(expected_size)));        
                end                
            else
                pos_ksde = mvksdensity([pox,poy],[xcen,ycen],'bandwidth',rmset.ssigma,'Function','pdf','Kernel','normal','BoundaryCorrection','reflection');
                pos_ksde = pos_ksde .* (rmset.binsize.^2) .* numel(pox) .* (1/rmset.srate); % convert to probability and then time      
                dwellmap = reshape(pos_ksde,size(X1));
            end  
            dwellmap(dwellmap < rmset.mindwell | isinf(abs(dwellmap))) = NaN;            
            speedlift = dwellmap;
            
            % Generate spikemap and ratemap      
            spk_ksde = mvksdensity([spx,spy],[xcen,ycen],'bandwidth',rmset.ssigma,'Function','pdf','Kernel','normal','BoundaryCorrection','reflection');        
            spk_ksde = spk_ksde .* (rmset.binsize.^2) .* numel(spx); % convert to probability and then spikes           
            spikemap = reshape(spk_ksde,size(X1));
            ratemap = spikemap ./ dwellmap;

            % make sure there are no weird values in the ratemap
            spikemap(isnan(dwellmap) | isinf(abs(spikemap))) = NaN;
            ratemap(isnan(dwellmap) | isinf(abs(ratemap))) = NaN;
            
%% >>>>>>>>>>>>>>>>>>>> Leutgeb kernel smoothed density estimate                            
        case {'leutgeb'}
% REF
% Leutgeb, Leutgeb, Barnes, Moser, McNaughton and Moser (2005) Independent Codes for Spatial and Episodic Memory in Hippocampal Neuronal Ensembles
% https://doi.org/10.1126/science.1114037
% Leutgeb, Leutgeb, Moser and Moser (2007) Pattern separation in the dentate gyrus and CA3 of the hippocampus
% https://doi.org/10.1126/science.1135801

% DESCRIPTION IN REF
% Spatial firing rate distributions (‘place fields’) for each wellisolated neuron were constructed in the standard manner by 
% summing the total number of spikes that occurred in a given location bin (5 by 5 cm), dividing by the amount of time that
% the animal spent in that location, and smoothing with a Gaussian centered on each bin. The average rate in each bin x was estimated as:
% [equation in paper]
% where g is a smoothing kernel, h is a smoothing factor, n is the number of spikes, si the location of the i-th spike, y(t) the location 
% of the rat at time t, and [0, T) the period of the recording. A Gaussian kernel was used for g and h = 5 cm. Positions more than 5 cm away 
% from the tracked path were regarded as unvisited.         

% SIMPLE DESCRIPTION
% For every bin we calculate the distance to every spike and every position data point. These distances are weighted with a gaussian so that
% data close to the bin centre = high and data far away = low. The sum of the spike distances are then divided by the sum of the position
% distances to obtain firing rate. The standard deviation of the gaussian can be changed, with a result similar to more smoothing in a standard
% firing rate map. For speed the effect of the gaussian can be clipped at config.maxdist, meaning we have to perform computations only
% on position data close to each bin
            
            % prepare bin centres
            [X1,Y1] = meshgrid(movmean(xvec,2,'EndPoints','discard'),movmean(yvec,2,'EndPoints','discard')); 
            xcen = X1(:);
            ycen = Y1(:);
            
            % prepare kernel
            coeff1 = (sqrt(2*pi) .* rmset.ssigma);
            coeff2 = (1 ./ (sqrt(2*pi) .* rmset.ssigma));
            k = @(x) ( (exp(-0.5 * (x./rmset.ssigma).^2) ./ coeff1) ./ coeff2 ); % gaussian function

            % the dwellmap only needs to be computed if one wasn't provided
            if ~isempty(speedlift) && ~all(isnan(speedlift(:)))
                dwellmap = speedlift; 
                expected_size = [length(yvec)-1,length(xvec)-1];
                if ~all(size(dwellmap) == expected_size)
                    error(sprintf('Size of provided dwellmap (%s) does not match the expected size {%s)... exiting',mat2str(size(dwellmap)),mat2str(expected_size)));        
                end                
            else
                pos_ksde = mvksdensity([pox,poy],[xcen,ycen],'bandwidth',1,'Function','pdf','Kernel',k,'BoundaryCorrection','reflection');%,'Support',rmset.maplims([1 2; 3 4]));
                pos_ksde = pos_ksde .* (rmset.binsize.^2) .* numel(pox) .* (1/rmset.srate); % convert to probability and then time      
                dwellmap = reshape(pos_ksde,size(X1));
            end  
            dwellmap(dwellmap < rmset.mindwell | isinf(abs(dwellmap))) = NaN;            
            speedlift = dwellmap;            
            
            % Generate spikemap and ratemap    
            spk_ksde = mvksdensity([spx,spy],[xcen,ycen],'bandwidth',1,'Function','pdf','Kernel',k,'BoundaryCorrection','reflection');%,'Support',rmset.maplims([1 2; 3 4]));        
            spk_ksde = spk_ksde .* (rmset.binsize.^2) .* numel(spx); % convert to probability and then spikes           
            spikemap = reshape(spk_ksde,size(X1));
            ratemap = spikemap ./ dwellmap;

            % make sure there are no weird values in the ratemap
            spikemap(isnan(dwellmap) | isinf(abs(spikemap))) = NaN;
            ratemap(isnan(dwellmap) | isinf(abs(ratemap))) = NaN;
            
%% >>>>>>>>>>>>>>>>>>>> Leutgeb 'pixelwise', or looping method                          
        case {'leutgeb_pixelwise'}
%% REF
% Leutgeb, Leutgeb, Barnes, Moser, McNaughton and Moser (2005) Independent Codes for Spatial and Episodic Memory in Hippocampal Neuronal Ensembles
% https://doi.org/10.1126/science.1114037
% Leutgeb, Leutgeb, Moser and Moser (2007) Pattern separation in the dentate gyrus and CA3 of the hippocampus
% https://doi.org/10.1126/science.1135801

%% DESCRIPTION IN REF
% Spatial firing rate distributions (‘place fields’) for each wellisolated neuron were constructed in the standard manner by 
% summing the total number of spikes that occurred in a given location bin (5 by 5 cm), dividing by the amount of time that
% the animal spent in that location, and smoothing with a Gaussian centered on each bin. The average rate in each bin x was estimated as:
% [equation in paper]
% where g is a smoothing kernel, h is a smoothing factor, n is the number of spikes, si the location of the i-th spike, y(t) the location 
% of the rat at time t, and [0, T) the period of the recording. A Gaussian kernel was used for g and h = 5 cm. Positions more than 5 cm away 
% from the tracked path were regarded as unvisited.         

%% SIMPLE DESCRIPTION
% For every bin we calculate the distance to every spike and every position data point. These distances are weighted with a gaussian so that
% data close to the bin centre = high and data far away = low. The sum of the spike distances are then divided by the sum of the position
% distances to obtain firing rate. The standard deviation of the gaussian can be changed, with a result similar to more smoothing in a standard
% firing rate map. For speed the effect of the gaussian can be clipped at config.maxdist, meaning we have to perform computations only
% on position data close to each bin

            % prepare bin centres
            [X1,Y1] = meshgrid(movmean(xvec,2,'EndPoints','discard'),movmean(yvec,2,'EndPoints','discard')); 
            xcen = X1(:);
            ycen = Y1(:);
            spikemap = zeros(size(X1));
    
            % prepare smoothing coefficients
            coeff1 = (sqrt(2*pi) .* rmset.ssigma);
            coeff2 = (1 ./ (sqrt(2*pi) .* rmset.ssigma));
            k = @(x) ( (exp(-0.5 * (x./rmset.ssigma).^2) ./ coeff1) ./ coeff2 ); % gaussian function
            
            % the dwellmap only needs to be computed if one wasn't provided
            dwell = false;
            if ~isempty(speedlift) && ~all(isnan(speedlift(:)))
                dwellmap = rmset.speedlift; 
                if ~all(size(dwellmap) == size(X1))
                    error(sprintf('Size of provided dwellmap (%s) does not match the expected size {%s)... exiting',mat2str(size(dwellmap)),mat2str(size(X1))));        
                end                
            else
                dwellmap = zeros(size(X1));                
                dwell = true;
            end  
            
            % Run through every bin            
            for bb = 1:numel(xcen) % for every bin  
                % generate dwellmap if necessary            
                if dwell % if we need a new dwellmap
                    % logical index which we can use to cut position data to within a square with side length config.maxdist of the bin
                    rindx = pox>xcen(bb)-rmset.maxdist & pox<xcen(bb)+rmset.maxdist & poy>ycen(bb)-rmset.maxdist & poy<ycen(bb)+rmset.maxdist;  
                    if ~any(rindx) % if there are no position data within this distance of the bin, leave it and the spikemap empty
                        continue
                    end                     
                    dp = sqrt(sum(([pox(rindx) poy(rindx)]-[xcen(bb) ycen(bb)]).^2,2)); % calculate the distance to every position data point
                    
                    % if this bin is too far from the position data, leave it and the spikemap empty
                    if min(dp,[],'omitnan') > rmset.mindist
                        continue
                    end 
                    dp_norm = k(dp(:)); % gaussian weight the distance values
                    
                    % calculate time as the sum of these distances multiplied by sampling rate to scale it correctly with time, accumulate this value in our dwellmap
                    tbin = sum(dp_norm,'omitnan')*(1/rmset.srate); % time in voxel
                    dwellmap(bb) = tbin; % accumulate data
                end

                % Create the spikemap in a similar way
                % logical index which we can use to cut spike data to within a square with side length config.maxdist of the bin                
                rindx = spx>xcen(bb)-rmset.maxdist & spx<xcen(bb)+rmset.maxdist & spy>ycen(bb)-rmset.maxdist & spy<ycen(bb)+rmset.maxdist;  
                if ~any(rindx) % if there are no spike data within this distance of the bin, leave it empty
                    continue
                end                   
                ds = sqrt(sum(([spx(rindx) spy(rindx)]-[xcen(bb) ycen(bb)]).^2,2)); % calculate the distance to every spike data point
                ds_norm = k(ds(:)); % gaussian weight the distance values
                
                % spike count is the sum of these distances, accumulate this value in our spikemap                
                sbin = sum(ds_norm,'omitnan'); % the number of spikes which fall into the bin
                spikemap(bb) = sbin;        
            end
            dwellmap(dwellmap < rmset.mindwell | isinf(abs(dwellmap))) = NaN;            
            ratemap = spikemap ./ dwellmap;

            % make sure there are no weird values in the ratemap
            spikemap(isnan(dwellmap) | isinf(abs(spikemap))) = NaN;
            ratemap(isnan(dwellmap) | isinf(abs(ratemap))) = NaN;

%% >>>>>>>>>>>>>>>>>>>> Bootstrap method
        case {'bootstrap'}  
            % the dwellmap only needs to be computed if one wasn't provided
            if ~isempty(speedlift) && ~all(isnan(speedlift(:)))
                dwellmap = speedlift; 
                expected_size = [length(yvec)-1,length(xvec)-1];                
                if ~all(size(dwellmap) == expected_size)
                    error(sprintf('Size of provided dwellmap (%s) does not match the expected size {%s)... exiting',mat2str(size(dwellmap)),mat2str(expected_size)));        
                end                
            else
                dwellmap = histcounts2(poy,pox,yvec,xvec) .* (1/rmset.srate);
                if rmset.smethod==1 && rmset.ssigma>0
                    dwellmap = imgaussfilt(dwellmap,rmset.ssigma./rmset.binsize,'FilterSize',rmset.ssize,'FilterDomain','spatial','Padding',0);
                end
            end                
            dwellmap(dwellmap < rmset.mindwell | isinf(abs(dwellmap))) = NaN; 
            speedlift = dwellmap;
                
            iti = 100;
            n = numel(spx);
            rmaps = NaN(length(yvec)-1,length(xvec)-1,iti);
            smaps = NaN(length(yvec)-1,length(xvec)-1,iti);            
            if rmset.smethod==2 && rmset.ssigma>0 % if we want to calculate ratio then smooth result
                H = fspecial('gaussian',rmset.ssize,rmset.ssigma./rmset.binsize); % kernel for smoothing                
            end           
            
            for ii = 1:iti % for every iteration
                rindx = randi(n,[n 1]);

                % generate spikemap and dwellmap if necessary
                % coordinates are backwards so the resulting map has the correct orientation in IJ
                spikemap_iti = histcounts2(spy(rindx),spx(rindx),yvec,xvec);

                % generate and smooth ratemap
                if rmset.smethod==1 && rmset.ssigma>0 % if we want to smooth spikes and time then calculate ratio
                    spikemap_iti = imgaussfilt(spikemap_iti,rmset.ssigma./rmset.binsize,'FilterSize',rmset.ssize,'FilterDomain','spatial','Padding',0);
                    ratemap_iti = spikemap_iti ./ dwellmap;    
                elseif rmset.smethod==2 && rmset.ssigma>0 % if we want to calculate ratio then smooth result
                    ratemap_iti = nanconv(spikemap_iti./dwellmap,H,'nanout');                       
                elseif rmset.smethod==3 || rmset.ssigma==0 % if we want no smoothing at all 
                    ratemap_iti = spikemap_iti ./ dwellmap;        
                end                
                rmaps(:,:,ii) = ratemap_iti;
                smaps(:,:,ii) = spikemap_iti;
            end
            ratemap = mean(rmaps,3,'omitnan');
            spikemap = mean(smaps,3,'omitnan');
            rmset.confidence95 = 1.960 .* ( mean(rmaps,3,'omitnan') ./ sqrt(iti) ); % mean 95% confidence interval
            rmset.coeffv = std(rmaps,[],3,'omitnan') ./ mean(rmaps,3,'omitnan'); % coefficient of variation
                
%% >>>>>>>>>>>>>>>>>>>> Adaptive (kernel accelerated method)                                   
        case {'fyhn'}  
% SIMPLE DESCRIPTION            


% REF for original method
% Fyhn et al (2004) Spatial Representation in the Entorhinal Cortex
% https://doi.org/10.1126/science.1099901

% DESCRIPTION IN REF for original method
% To characterize firing fields, we estimated a spike density function by convolving the
% spike train with a smoothing kernel and sampled this ‘signal’ synchronously with the
% tracker position. The kernel was a 2 s wide Blackman window normalized to have unit
% gain at zero frequency, which approximates a Gaussian kernel with 330 ms standard
% deviation. The reason for smoothing the rate map in the temporal dimension has to do
% with the low frequency of the tracker signal compared to the spike frequency. Our
% algorithm required the analog spike train to be discretized in synchrony with the video
% tracker (i.e. intervals at 20 ms). This necessitated the use of an anti-aliasing lowpass
% filter with substantial attenuation at 25 Hz (40 ms periods). After determining the firing
% rate for each tracker position, we calculated the spatially averaged firing rate for each
% location in a 5 × 5 cm grid centered above the open field. The average rates were
% estimated by weighted means of the sampled spike density function. The weights were
% Euclidian distance multiplied by a 30 cm wide Blackman window. In order to avoid
% error from extrapolation, we considered positions more than 5 cm away from the
% tracked path as unvisited.  

            time_window = rmset.twindow; % seconds
            samp_interval = 1/rmset.srate;
            pot = rmset.pot; 
            spt = rmset.spt; 

            edg = [pot; pot(end)+samp_interval] - (samp_interval/2);
            spike_train = histcounts(spt,edg); % binned spike train
            
            if time_window>samp_interval
                black_window = blackman(ceil(rmset.srate*time_window)); % blackman window
                black_window = black_window ./ sum(black_window(:),'omitnan'); % normalize for unit gain at zero frequency
                insta_spikes = conv(spike_train,black_window,'same'); % smoothed spike count (will sum to numel(spx))
                insta_rate = insta_spikes ./ (1/rmset.srate); % instantaneous firing rate
            else
                insta_rate = spike_train ./ (1/rmset.srate); % instantaneous firing rate                
            end

            % prepare bin centres
            [X1,Y1] = meshgrid(movmean(xvec,2,'EndPoints','discard'),movmean(yvec,2,'EndPoints','discard')); 
            xcen = X1(:);
            ycen = Y1(:);
            ratemap = NaN([size(X1),length(rmset.ssigma)]);
            spikemap = NaN(size(X1));
            dwellmap = NaN(size(X1));

            % Run through every bin            
            for bb = 1:numel(xcen) % for every bin  
                % logical index which we can use to cut position data to within a square with side length rmset.binsize/2 of the bin
                rindx = pox>xcen(bb)-rmset.maxdist & pox<xcen(bb)+rmset.maxdist & poy>ycen(bb)-rmset.maxdist & poy<ycen(bb)+rmset.maxdist;                  
                if ~any(rindx) % if there are no position data within this distance of the bin, leave it empty
                    continue
                end                     
                dp = sqrt(sum(([pox(rindx) poy(rindx)]-[xcen(bb) ycen(bb)]).^2,2)); % calculate the distance to every position data point
                fr = insta_rate(rindx); % firing rate values for the samples within the bin

                % if this bin is too far from the position data, leave it empty
                if min(dp,[],'omitnan') > rmset.mindist
                    continue
                end 

                [i,j] = ind2sub(size(X1),bb);
                for ss = 1:length(rmset.ssigma) % for every smoothing value
                    dp_norm = normpdf(dp,0,rmset.ssigma(ss)); % gaussian weight the distance values

                    % calculate firing rate as the weighted sum of these rates (weighted by distance to bin)
                    ratemap(i,j,ss) = sum(fr(:).*dp_norm(:),'omitnan'); % firing rate in bin
                end
            end

        otherwise
                error(sprintf('Unknown mapping method {%s)... exiting',rmset.method));        


    end




















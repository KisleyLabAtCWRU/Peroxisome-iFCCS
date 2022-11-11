%% PeroxisomeCorrelationMain_LK20221006.m
% image analysis to perform image correlation on mCherry and GFP-labeled
% peroxisome proteins
% Analysis includes sections to:
% 1. Load data (read comments in this section about the input data format)
% 2. Image registration of GFP and mCherry channels based on control point selection,
% correlation, and image transformation using fiducial bead markers
% 3. Find location of peroxisomes using existing particle tracking code
% (requires Troika single particle tracking code maintained by Prof.
% Christy Landes' group available at:
% https://github.com/LandesLab/Troika-Single-particle-tracking; relevant
% citation is DOI: 10.1039/C3CP53968G
% 4. For each peroxisome, perform image fluorscence cross correlation
% spectroscopy at the 17 pixels surrounding the centroid location (read
% comments in this section regarding format of the data saved)
% 5. Save data (optional) 
% 6. Make figures of results (optional)
% Contact Prof. Lydia Kisley, Case Western Reserve University,
% lydia.kisley@case.edu for any questions


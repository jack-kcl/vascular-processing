# vascular-processing


Rayburst_radius_est.m
Estimates the radii thoughout the segmented vascular images using the Rayburst algorithm.

AnalyseSkeleton.m
Extracts the segment and junction information from skeletonised vascular images. Preserves loops.

import_skeleton.m
Imports the nodes, segments, knots into Matlab from xls.

analyse/

calculate_conductance.m
Computes the global fluid conductance from a given vascular network.

label_segments.m
Labels each segment as being radial, spiral or canal vessel.

postprocess.m
Remove spurious branches and topological artefacts for final analysis.

report_stats.m
Calculate and report geometric caharacteristics of the network, by segment.

report_stats_nodes.m
Caclulate and report geometric characterisitics of the network, by nodes.


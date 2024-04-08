# PTC_feature_extraction

use: 
bounds = mask2bounds(PTC_mask_normal);
[feats] = extract_all_features(bounds);

Complete Feature list: 
Voronoi: Area std.
Voronoi: Area average
Voronoi: Area 5% / 95%
Voronoi: Area disorder
Voronoi: Perimeter std.
Voronoi: Perimeter average
Voronoi: Perimeter 5% / 95%
Voronoi: Perimeter disorder
Voronoi: Chord std.
Voronoi: Chord average
Voronoi: Chord 5% / 95%
Voronoi: Chord disorder
Delaunay: Side length 5% / 95%
Delaunay: Side length std.
Delaunay: Side length average
Delaunay: Side length disorder
Delaunay: Triangle area 5% / 95%
Delaunay: Triangle area std.
Delaunay: Triangle area average
Delaunay: Triangle area disorder
MST: Edge length average
MST: Edge length std.
MST: Edge length 5% / 95%
MST: Edge length disorder
Arch: Area of polygons
Arch: Number of polygons
Arch: Density of polygons
Arch: Average distance to 3 nearest neighbors
Arch: Average distance to 5 nearest neighbors
Arch: Average distance to 7 nearest neighbors
Arch: Std. distance to 3 nearest neighbors
Arch: Std. distance to 5 nearest neighbors
Arch: Std. distance to 7 nearest neighbors
Arch: Disorder of distance to 3 nearest neighbors
Arch: Disorder of distance to 5 nearest neighbors
Arch: Disorder of distance to 7 nearest neighbors
Arch: Avg. nearest neighbors in a 10 pixel radius
Arch: Avg. nearest neighbors in a 20 pixel radius
Arch: Avg. nearest neighbors in a 30 pixel radius
Arch: Avg. nearest neighbors in a 40 pixel radius
Arch: Avg. nearest neighbors in a 50 pixel radius
Arch: Std. neighbors in 10 pixel radius
Arch: Std. neighbors in 20 pixel radius
Arch: Std. neighbors in 30 pixel radius
Arch: Std. neighbors in  40 pixel radius
Arch: Std. neighbors in 50 pixel radius
Arch: Disorder neighbors in 10 pixel radius
Arch: Disorder neighbors in 20 pixel radius
Arch: Disorder neighbors in 30 pixel radius
Arch: Disorder neighbors in 40 pixel radius
Arch: Disorder neighbors in 50 pixel radius
Shape: Mean area ratio
Shape: Mean distance ratio
Shape: Mean Std. of distance
Shape: Mean variance of distance
Shape: Mean distance ratio
Shape: Mean perimeter ratio
Shape: Mean smoothness
Shape: Mean invariant moment 1
Shape: Mean invariant moment 2
Shape: Mean invariant moment 3
Shape: Mean invariant moment 4
Shape: Mean invariant moment 5
Shape: Mean invariant moment 6
Shape: Mean invariant moment 7
Shape: Mean fractal dimension
Shape: Mean Fourier descriptor 1
Shape: Mean Fourier descriptor 2
Shape: Mean Fourier descriptor 3
Shape: Mean Fourier descriptor 4
Shape: Mean Fourier descriptor 5
Shape: Mean Fourier descriptor 6
Shape: Mean Fourier descriptor 7
Shape: Mean Fourier descriptor 8
Shape: Mean Fourier descriptor 9
Shape: Mean Fourier descriptor 10
Shape: Std. area ratio
Shape: Std. distance ratio
Shape: Std. Std. of distance
Shape: Std. variance of distance
Shape: Std. distance ratio
Shape: Std. perimeter ratio
Shape: Std. smoothness
Shape: Std. invariant moment 1
Shape: Std. invariant moment 2
Shape: Std. invariant moment 3
Shape: Std. invariant moment 4
Shape: Std. invariant moment 5
Shape: Std. invariant moment 6
Shape: Std. invariant moment 7
Shape: Std. fractal dimension
Shape: Std. Fourier descriptor 1
Shape: Std. Fourier descriptor 2
Shape: Std. Fourier descriptor 3
Shape: Std. Fourier descriptor 4
Shape: Std. Fourier descriptor 5
Shape: Std. Fourier descriptor 6
Shape: Std. Fourier descriptor 7
Shape: Std. Fourier descriptor 8
Shape: Std. Fourier descriptor 9
Shape: Std. Fourier descriptor 10
Shape: Median area ratio
Shape: Median distance ratio
Shape: Median Std. of distance
Shape: Median variance of distance
Shape: Median distance ratio
Shape: Median perimeter ratio
Shape: Median smoothness
Shape: Median invariant moment 1
Shape: Median invariant moment 2
Shape: Median invariant moment 3
Shape: Median invariant moment 4
Shape: Median invariant moment 5
Shape: Median invariant moment 6
Shape: Median invariant moment 7
Shape: Median fractal dimension
Shape: Median Fourier descriptor 1
Shape: Median Fourier descriptor 2
Shape: Median Fourier descriptor 3
Shape: Median Fourier descriptor 4
Shape: Median Fourier descriptor 5
Shape: Median Fourier descriptor 6
Shape: Median Fourier descriptor 7
Shape: Median Fourier descriptor 8
Shape: Median Fourier descriptor 9
Shape: Median Fourier descriptor 10
Shape: 5% / 95% area ratio
Shape: 5% / 95% distance ratio
Shape: 5% / 95% std. of distance
Shape: 5% / 95% variance of distance
Shape: 5% / 95% distance ratio
Shape: 5% / 95% perimeter ratio
Shape: 5% / 95% smoothness
Shape: 5% / 95% invariant moment 1
Shape: 5% / 95% invariant moment 2
Shape: 5% / 95% invariant moment 3
Shape: 5% / 95% invariant moment 4
Shape: 5% / 95% invariant moment 5
Shape: 5% / 95% invariant moment 6
Shape: 5% / 95% invariant moment 7
Shape: 5% / 95% fractal dimension
Shape: 5% / 95% Fourier descriptor 1
Shape: 5% / 95% Fourier descriptor 2
Shape: 5% / 95% Fourier descriptor 3
Shape: 5% / 95% Fourier descriptor 4
Shape: 5% / 95% Fourier descriptor 5
Shape: 5% / 95% Fourier descriptor 6
Shape: 5% / 95% Fourier descriptor 7
Shape: 5% / 95% Fourier descriptor 8
Shape: 5% / 95% Fourier descriptor 9
Shape: 5% / 95% Fourier descriptor 10
CPT: Mean tensor contrast energy
CPT: Std. tensor contrast energy
CPT: Range tensor contrast energy
CPT: Mean tensor contrast inverse moment
CPT: Std. tensor contrast inverse moment
CPT: Range tensor contrast inverse moment
CPT: Mean tensor contrast average
CPT: Std. tensor contrast average
CPT: Range tensor contrast average
CPT: Mean tensor contrast var
CPT: Std. tensor contrast var
CPT: Range tensor contrast var
CPT: Mean tensor contrast entropy
CPT: Std. tensor contrast entropy
CPT: Range tensor contrast entropy
CPT: Mean tensor intensity average
CPT: Std. tensor intensity average
CPT: Range tensor intensity average
CPT: Mean tensor intensity variance
CPT: Std. tensor intensity variance
CPT: Range tensor intensity variance
CPT: Mean tensor intensity entropy
CPT: Std. tensor intensity entropy
CPT: Range tensor intensity entropy
CPT: Mean tensor entropy
CPT: Std. tensor entropy
CPT: Range tensor entropy
CPT: Mean tensor energy
CPT: Std. tensor energy
CPT: Range tensor energy
CPT: Mean tensor correlation
CPT: Std. tensor correlation
CPT: Range tensor correlation
CPT: Mean tensor information measure 1
CPT: Std. tensor information measure 1
CPT: Range tensor information measure 1
CPT: Mean tensor information measure 2
CPT: Std. tensor information measure 2
CPT: Range tensor information measure 2
Sub-Graph: Number of nodes
Sub-Graph: Number of edges
Sub-Graph: Average degree
Sub-Graph: Average eccentricity
Sub-Graph: Diameter
Sub-Graph: Radius
Sub-Graph: Average eccentricity 90%
Sub-Graph: Diameter 90%
Sub-Graph: Radius 90%
Sub-Graph: Average path length
Sub-Graph: Clustering coefficient c
Sub-Graph: Clustering coefficient d
Sub-Graph: Clustering coefficient e
Sub-Graph: Number of connected components
Sub-Graph: Giant connected component ratio
Sub-Graph: Average connected component size
Sub-Graph: Number isolated nodes
Sub-Graph: Percentage isolated nodes
Sub-Graph: Number end nodes
Sub-Graph: Percentage end nodes
Sub-Graph: Number central nodes
Sub-Graph: Percentage central nodes
Sub-Graph: Mean edge length
Sub-Graph: Std. edge length
Sub-Graph: Skewness edge length
Sub-Graph: Kurtosis edge length
mean PTC size
mean PTC eccentricity
mean PTC major axis
mean PTC minor axis
mean PTC Aspect ratio
mean PTC extent
mean PTC orientation
std PTC size
std PTC eccentricity
std PTC major axis
std PTC minor axis
std PTC Aspect ratio
std PTC extent
std PTC orientation
The 39 co-occurring PTC tensor (CPT) features measured the disorder of neighborhoods of PTCs using the entropy of orientation of the major axes of PTCs within a local neighborhood. 
The 26 subgraph features described the connectivity and clustering of small PTC neighborhoods using PTC centroids. 

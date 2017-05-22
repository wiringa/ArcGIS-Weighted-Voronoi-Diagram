# ArcGIS Weighted Voronoi Diagram

Generate a weighted Voronoi diagram (Thiessen polygons) from a set of
points and attribute values, using an intermediate raster for
distance calculations. Instead of assigning raster cells to source
points based solely on distance, the distance is divided by some
numeric attribute provided with the point data.

## Parameters/Usage

All parameters are required.

* **Point source** is a point feature class containing the points
around which to build polygons. It must contain at least one
numeric field, the weight field. The feature class must be in a
geographic coordinate system.
* **Weight field** is a numeric field in the point source that the
distance between any two raster cells is divided by to determine
assignments.
* **Diagram output feature class** contains the resulting polygon
features. OIDs from the source are carried through to the output.
* **Buffer point source** is the number of decimal degrees to
buffer the extent of source by.
* **Raster cell size** is the number of decimal degrees in both
X and Y to use when creating the assignment raster.
* **Geodesic distance methods** must be one of `Haversine` or
`Vincenty`. Distances between cells are calculated between
the latitudes and longitudes of cell centers. The Haversine
method calculates distance on the spheroid (assuming 6371 km
as the average earth radius) while the Vincenty method
calculates distance on the ellipsoid (using WGS84 datum values
for semi-major and semi-minor axis lengths). Note that the
Haversine method is significantly faster than the Vincenty
method, roughly 6x-8x.

## Environment variables

The `arcpy.env.extent` environment variable is temporarily
overridden to allow for buffering of the point source extent.

## Python packages included

The [Haversine 0.4.5](https://pypi.python.org/pypi/haversine)
and [Vincenty 0.1.4](https://pypi.python.org/pypi/vincenty)
Python packages are included with the toolbox.

## Compatibility

Tested using ArcGIS for Desktop 10.4.1 using 64-bit background
geoprocessing.

## Known Issues and Limitations

* It requires the source data to be in a GCS (e.g., WGS84 or NAD83).
* It's slow for higher-resolution runs. Better use of NumPy may
help.
* Both the Haversine and Vincenty geodesic distance methods have
fixed parameters for earth radius and ellipsoid axis lengths.

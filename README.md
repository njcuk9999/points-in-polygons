# points-in-polygons
Python code to count and identify points (xi,yi) inside a set of polygons

    Takes a list of polygons and counts how many points are inside.
    Points are defined by xarr, yarr such that we look whether
    (xarr[i], yarr[i]) is contained in the polygons.


    :param polygons: list of polygon objects, a polygon collection: this is
                     the list or array of polygons, each polygon should be
                     a list of arrays that contain a set of vertices each
                     with a list of (x, y) coordinates:

                     polygons = [polygon1, polygon2, ... , polygonN]

                     where:

                     polygon1 = [vertexlist1, vertexlist2, ..., vertexlistN]

                     vertexlist1 = [(x0, y0), (x1, y1), ...., (xN, yN)]

                     i.e. a single polygon (a square) could be:
                         [[[(0, 0), (1, 0), (1, 1), (0, 1)]]]

    :param xarr: array of floats, x coordinates for points
                 (same length as yarr)
    :param yarr: array of floats, y coordinates for points
                 (same length as xarr)

    :param polygons_in: Array of bools Same shape as polygons
                        (except no vertices). This controls whether a polygon
                        is "inside" another polygon if include_holes is False
                        we assume any polygon with polygon_in = True is a hole
                        and thus the count should take these points as NOT
                        being in the polygons.

                        polygons_in = None

                        polygons_in = [polygon_in1, polygon_in2, ...,
                                       polygon_inN]

                        where:

                        polygon_in = [True, False, ..., True]

    :param include_holes: bool, whether we count polygons that lie inside
                          (i.e. have polygons_in[i][j] = True) as holes
                          and thus points inside these polygons are not
                          counted as inside the polygons


    :return count: array of counts for each polygon (in polygon)
                   one count for each polygon in polygons

    :return xpoints: array of x points where both xarr[i] and yarr[i] inside
                     the polygon
                     one array of points for each polygon in polygons

    :return ypoints: array of y points where both xarr[i] and yarr[i] inside
                     the polygon
                     one array of points for each polygon in polygons

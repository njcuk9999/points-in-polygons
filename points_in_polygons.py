import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib.path as mplPath


# =============================================================================
# Define functions
# =============================================================================
def mask_from_polygons(xarr, yarr, polygons, polygons_in):
  """
  Takes a list of polygons and returns a array of bools (mask) where [i]
  is True if xarr[i] and yarr[i] are inside all polygons.
  
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


  :return insideany: array of bools for, mask True where xarr[i] and yarr[i] are inside the polygon
  """
  xyarr = np.array(zip(xarr, yarr))
  falsearray = np.array([False] * len(xarr), dtype=bool)
  insideany = falsearray.copy()
  for k in range(len(polygons)):
      polygon = polygons[k]
      # +++++++++++++++++++++++++++++++++++++++++++++
      # deal with annoying contained polygons
      if polygons_in is None:
          polygon_in = None
      else:
          polygon_in = polygons_in[k]
      # +++++++++++++++++++++++++++++++++++++++++++++
      for j in range(len(polygon)):
          poly = polygon[j]
          # mask out the points outside the poly clip box
          pmax_x, pmax_y = np.max(poly, axis=0)
          pmin_x, pmin_y = np.min(poly, axis=0)
          mask1 = (xarr > pmin_x) & (xarr < pmax_x)
          mask2 = (yarr > pmin_y) & (yarr < pmax_y)
          mask = mask1 & mask2
          # if no points inside poly clip box then don't bother counting
          if len(mask[mask]) == 0:
              continue
          # -----------------------------------------------------------------
          # deal with annoying contained polygons
          if polygon_in is None:
              poly_in = False
          else:
              poly_in = polygon_in[j]
          # -----------------------------------------------------------------
          # Creates a mask for points (xs, ys) based on whether they are
          # inside a polygon poly from http://stackoverflow.com/a/23453678
          bbPath = mplPath.Path(poly)
          inside = falsearray.copy()
          inside[mask] = bbPath.contains_points(xyarr[mask])
          # -----------------------------------------------------------------
          # if polygon is inside another polygon (as defined by poly_in)
          # do not count it as inside
          insideany |= inside
          # -----------------------------------------------------------------
  return insideany

# count the number of points inside a list of "polygon collections"
def count_objects_inside(polygons, xarr, yarr, polygons_in=None,
                         include_holes=True):
    """
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
    """
    # zip the x and y coordinates
    xyarr = np.array(zip(xarr, yarr))
    # set up the counts
    counts = np.zeros(len(polygons), dtype=int)
    # set up a blank object array for filling with points that are inside
    oix, oiy = np.zeros((2, len(polygons)), dtype=object)
    # set up the false array to copy mask into
    # (quicker to do this once and copy each time rather than make each time)
    falsearray = np.array([False] * len(xarr), dtype=bool)
    # loop round each polygon in polygons
    for k in range(len(polygons)):
        # get this polygon set
        polygon = polygons[k]
        # deal with annoying contained polygons
        if polygons_in is None:
            polygon_in = None
        else:
            polygon_in = polygons_in[k]
        # copy inside any from previously created false array
        insideany = falsearray.copy()
        # loop round each set of vertices in this polygon
        for j in range(len(polygon)):
            # get the vertices for this polygon
            vertices = polygon[j]
            # mask out the points outside the poly clip box (this is quicker
            # than running the full set of points through Path.contains_points)
            pmax_x, pmax_y = np.max(vertices, axis=0)
            pmin_x, pmin_y = np.min(vertices, axis=0)
            mask1 = (xarr > pmin_x) & (xarr < pmax_x)
            mask2 = (yarr > pmin_y) & (yarr < pmax_y)
            mask = mask1 & mask2
            # if no points inside poly clip box then don't bother counting
            if len(mask[mask]) == 0:
                continue
            # deal with polygons contained within polygons
            # (must be known before runtime)
            if polygon_in is None:
                poly_in = False
            else:
                poly_in = polygon_in[j]
            # Creates a mask for points (xs, ys) based on whether they are
            # inside a polygon poly from http://stackoverflow.com/a/23453678
            bbPath = mplPath.Path(vertices)
            # inside = point_inside_polygon_multi3(xarr, yarr, poly)
            inside = falsearray.copy()
            inside[mask] = bbPath.contains_points(xyarr[mask])
            # if polygon is inside another polygon (as defined by poly_in)
            # do not count it as inside
            if poly_in and not include_holes:
                insideany &= ~inside
            else:
                insideany |= inside
        # store counts and points in created arrays
        counter = len(insideany[insideany])
        oix[k] = xarr[insideany]
        oiy[k] = yarr[insideany]
        counts[k] = counter
    # return counts, array of points for each polygon (x and y separately)
    return np.array(counts), oix, oiy


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    xmin, xmax, ymin, ymax, total = -2, 8.5, -2, 8.5, 25000
    # create some polygons
    # polygons = [polygon1, polygon2, ... polygonN]
    # polygon1 = [vertexset1, vertexset2, ... vertexsetN]
    # vertexset1 = [(x0, y0), (x1, y1), ..., (xN, yN)]
    # This should produce two polygons (a diamond with a diamond hole and a
    # rectangle with a rectangle hole)
    my_polygons = [[[(0, -1), (-1, 0), (0, 1), (1, 0)],
                    [(0, -0.7), (-0.7, 0), (0, 0.7), (0.7, 0.0)]],
                   [[(2, 2), (8, 2), (8, 8), (2, 8)],
                    [(4, 4), (4, 5), (5, 5), (5, 4)]]]
    # define whether polygons are inside any other polygon
    # polygons_in = [polygon_in1, polygon_in2, ..., polygon_inN]
    # polygon_in1 = [True, False, ..., True]
    my_polygons_in = [[False, True], [False, True]]
    # create random grid of points
    xs = np.random.uniform(xmin, xmax, total)
    ys = np.random.uniform(ymin, ymax, total)
    # uses the count function
    args = [my_polygons, xs, ys, my_polygons_in]
    count, xpoints, ypoints = count_objects_inside(*args, include_holes=False)
    # plots a graphical representation of the function
    plt.close()
    fig, frame = plt.subplots(ncols=1, nrows=1)
    # loop around my_polygons and plot them in blue (unless they are "inside"
    # then plot them in white)
    for ki, my_polygon in enumerate(my_polygons):
        for ji, my_vertices in enumerate(my_polygon):
            if my_polygons_in[ki][ji]:
                p = mpatch.Polygon(my_vertices, color='w', zorder=2)
            else:
                p = mpatch.Polygon(my_vertices, color='b', zorder=1)
            frame.add_artist(p)
    # add the original points to the plot in black
    frame.scatter(xs, ys, color='k', zorder=3)
    # loop round xpoints and plot them in red
    for ii in range(len(xpoints)):
        frame.scatter(xpoints[ii], ypoints[ii], color='r', zorder=4)
    # display the number of points as the title
    args = [np.sum(count), len(xs.flatten())]
    msg = 'There are {0}/{1} points inside the area occupied by the polygons'
    plt.title(msg.format(*args))
    # force the limits of the axis to the size of the points
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    # finally show and close our figure
    plt.show()
    plt.close()
# =============================================================================
# End of code
# =============================================================================

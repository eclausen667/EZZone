import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import Delaunay
import math
import pyproj
from shapely.geometry import MultiPolygon, Polygon, asShape, LineString, mapping, MultiPoint, Point, MultiLineString, GeometryCollection
from shapely.ops import cascaded_union, polygonize
from shapely.affinity import rotate
from rasterstats import zonal_stats
import geojson
from jenks import jenks
import brewer2mpl
import matplotlib.pyplot as plt
from descartes import PolygonPatch
from matplotlib.collections import PatchCollection

def extrapolate_nans(x, y, v):
    '''
    Extrapolate the NaNs or masked values in a grid INPLACE using nearest
    value.

    .. warning:: Replaces the NaN or masked values of the original array!

    Parameters:

    * x, y : 2D arrays
        Arrays with the x and y coordinates of the data points.
    * v : 2D array
        Array with the scalar value assigned to the data points.

    Returns:

    * v : 2D array
        The array with NaNs or masked values extrapolated.
    '''

    if np.ma.is_masked(v):
        nans = v.mask
    else:
        nans = np.isnan(v)
    notnans = np.logical_not(nans)
    try:
        v[nans] = griddata((x[notnans], y[notnans]), v[notnans],
        (x[nans], y[nans]), method='nearest').ravel()
        return v
    except ValueError:
        return

def alpha_shape(points, alpha):
    """
    Compute the alpha shape (concave hull) of a set
    of points.
 
    @param points: Iterable container of points.
    @param alpha: alpha value to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    """
    if len(points) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return MultiPoint(list(points)).convex_hull
 
    def add_edge(edges, edge_points, coords, i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add( (i, j) )
        edge_points.append(coords[ [i, j] ])
 
    coords = np.array([point.coords[0]
                       for point in points])
 
    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the
    # triangle
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
 
        # Lengths of sides of triangle
        a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
 
        # Semiperimeter of triangle
        s = (a + b + c)/2.0
 
        # Area of triangle by Heron's formula
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))
        circum_r = a*b*c/(4.0*area)
 
        # Here's the radius filter.
        #print circum_r
        if circum_r < 1.0/alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)
 
    m = MultiLineString(edge_points)
    triangles = list(polygonize(m))
    return cascaded_union(triangles), edge_points
    
def geojFCparser(obj, geomtype):
    """Parse a geometry from a GeoJSON FeatureCollection
    
    Parameters
    __________
        
        obj : dict
            the FeatureCollection object as type dict, or will be loaded to dict
            
        geomtype : str
            the desired geometry type to extract
            must be one of the following:("Point", "MultiPoint", "LineString", "MultiLineString", "Polygon", "MultiPolygon")
    
    Returns
    _______
        
        dict of geometry object if inside FeatureCollection, else string error message
    """        
    if geomtype not in ("Point", "MultiPoint", "LineString", "MultiLineString", "Polygon", "MultiPolygon"):
        return 'Invalid GeomType'
    #if our object is a raw string
    if type(obj) == str:
        obj = geojson.loads(obj)
    for i in obj['features']:
        for v in i.iteritems():
            for x in v:
                if isinstance(x,dict) and geomtype in x.values():
                    return x
    return 'GeomType not in FeatureCollection'

def boundJson(x,y, geometryType = 1, inProjection = 'epsg:4326' ,outProjection = 'epsg:3857'):
    """
    Create a GeoJSON boundary from csv input
    Args:
        filename(file): the filename to perform the operation on
        geometryType(int): type of geometric shape, 1 = rectangular, 2 = circular, 3 = irregular
        inProjection(int): EPSG projection number for csv data
        outProjection(int): EPSG projection number for output GeoJSON string (default 3857 for google maps)

    Returns:
        GeoJSON string of polygon boundary 
    """
    g = geometryType
    inP = inProjection
    outP = outProjection
    jFeat = []
    crs = {
    "type": "name",
    "properties": {
        "name": outP
        }}

    inP = pyproj.Proj(init=inP)
    outP = pyproj.Proj(init=outP)


    #get the point coordinates
    points = np.column_stack([x,y])

    if g == 1:
        xMax = [x.max(),y[np.where(x==x.max())].item(0)]
        xMax = list(pyproj.transform(inP,outP,xMax[0],xMax[1]))

        yMax = [x[np.where(y==y.max())].item(0),y.max()]
        yMax = list(pyproj.transform(inP,outP,yMax[0],yMax[1]))

        xMin = [x.min(),y[np.where(x==x.min())].item(0)]
        xMin = list(pyproj.transform(inP,outP,xMin[0],xMin[1]))

        yMin = [x[np.where(y==y.min())].item(0),y.min()]
        yMin = list(pyproj.transform(inP,outP,yMin[0],yMin[1]))

        geoJson = geojson.Polygon([[xMax,yMax,xMin,yMin,xMax]],crs=crs)

    if g == 2 or g == 3:
        points = MultiPoint(points)
        m = mapping(points.convex_hull)
        for i in m['coordinates']:
            for j in i:
                jFeat.append(pyproj.transform(inP,outP,j[0],j[1]))

        geoJson = geojson.Polygon([jFeat], crs= crs)

    #return geoJson
    return geojson.loads(geojson.dumps(geoJson))#, points
    
def getSpacing(boundary, rowSpacing):
    new = rowSpacing
    while ((math.sqrt((boundary.area)/1000)) > new):
        new = new + rowSpacing
    return new


#orders y elements least to greatest of a 2d numpy array and eliminates duplicate values, 
#returning the new array
def unique(a):
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1)
    return a[ui]

#returns the difference between two x points
def deltax(xy1, xy2):
    return (xy2[0] - xy1[0])

#returns integer of epsg coordinates:
def epsgInt(stringEpsg):
    li = []
    for l in stringEpsg:
        if l.isdigit():
            li.append(l)
    return int(''.join(li))

#reproject geojson object/string, shapely polygon/multip and return a shapely polygon, 
#can pass projection as string or pyproj object
def reprojectPoly(obj, inProjection = 'epsg:3857', outProjection = 'epsg:4326'):
    
    #if our object is a raw string
    if type(obj) == str:
        obj = geojson.loads(obj)
    #if our object is a geojson dict
    if type(obj) == dict:
        #is our geojson a feature class?
        try:
            poly = asShape(obj['geometry'])
        #or a regular geometry?
        except KeyError:
            poly = asShape(obj)
    else:
        poly = obj
        
    #pyproj transformation
    if type(inProjection) != str:
        inP = inProjection
    else:
        inP = pyproj.Proj(init=inProjection)
    if type(outProjection) != str:
        outP = outProjection
    else:
        outP = pyproj.Proj(init=outProjection)
       
    #transform coords and return polygon
    extCoords = [pyproj.transform(inP,outP,i[0],i[1]) for i in list(poly.exterior.coords)]
    try:
        intCoords = [[pyproj.transform(inP,outP,j[0],j[1]) for j in list(i.coords)] for i in poly.interiors]
        poly = Polygon(extCoords,intCoords)
    except AttributeError:
        poly = Polygon(extCoords)
    return poly
    
def reprojectMultiP(obj, inProjection = 'epsg:3857', outProjection = 'epsg:4326'):
    #array for new multipolygons
    multiPolys = []
    #iterate through list of multipolygons
    for polys in obj:
        #array for reprojected polygons
        newPolys = []
        #try to iterate through multipolygon, except when it is only a polygon
        try:
            for p in polys:
                newPolys.append(reprojectPoly(p, inProjection, outProjection))
            #append multipolygon to list of multipolygons
            multiPolys.append(MultiPolygon(newPolys))
        except TypeError:
            multiPolys.append(reprojectPoly(polys, inProjection, outProjection))
    return multiPolys

#reproject geojson/shapely point
def reprojectPoint(obj, inProjection = 'epsg:3857', outProjection = 'epsg:4326'):
    
    #if our object is a raw string
    if type(obj) == str:
        obj = geojson.loads(obj)
    #if our object is a geojson dict
    if type(obj) == dict:
        #is our geojson a feature class?
        try:
            point = asShape(obj['geometry'])
        #or a regular geometry?
        except KeyError:
            point = asShape(obj)
    else:
        point = obj
        
    #pyproj transformation
    if type(inProjection) != str:
        inP = inProjection
    else:
        inP = pyproj.Proj(init=inProjection)
    if type(outProjection) != str:
        outP = outProjection
    else:
        outP = pyproj.Proj(init=outProjection)
        
    newCoords =  [pyproj.transform(inP,outP,i[0],i[1]) for i in list(point.coords)]
    return Point(newCoords) 
    
#reproject geojson/shapely point
def reprojectLine(obj, inProjection = 'epsg:3857', outProjection = 'epsg:4326'):
    
    #if our object is a raw string
    if type(obj) == str:
        obj = geojson.loads(obj)
    #if our object is a geojson dict
    if type(obj) == dict:
        #is our geojson a feature class?
        try:
            line = asShape(obj['geometry'])
        #or a regular geometry?
        except KeyError:
            line = asShape(obj)
    else:
        line = obj
        
    #pyproj transformation
    if type(inProjection) != str:
        inP = inProjection
    else:
        inP = pyproj.Proj(init=inProjection)
    if type(outProjection) != str:
        outP = outProjection
    else:
        outP = pyproj.Proj(init=outProjection)
        
    newCoords =  [pyproj.transform(inP,outP,i[0],i[1]) for i in list(line.coords)]
    return LineString(newCoords)

#reproject arrays of separate x and y values    
def reprojectArray(x,y,inP,outP):
    #stack x, y coords
    points = np.column_stack([x,y])
    
    #pyproj transformation
    if type(inP) == str:
        inP = pyproj.Proj(init = inP)
    if type(outP) == str:
        outP = pyproj.Proj(init = outP)
    #2d array of coords reprojected
    return np.array([pyproj.transform(inP,outP,p[0],p[1]) for p in points])

    
#return a pyproj projection object given an array of coordinates or shapely polygon in lat/long dd
def utmzone(xy):
    #if input is shapely polygon
    if isinstance(xy, Polygon):
        xy = xy.exterior.coords[0]
    #calculate zone number, include 1 if southern hemisphere
    zone = (int(math.floor((xy[0] +180)/6) + 1), 1 if xy[1] < 0 else 0)
    #if in southern hemisphere
    if zone[1] == 1:
        return pyproj.Proj('+proj=utm +zone='+ str(zone[0]) + ', +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
    return pyproj.Proj('+proj=utm +zone='+ str(zone[0]) + ', +ellps=WGS84 +datum=WGS84 +units=m +no_defs')


#create a rectangular grid of polygons given a shapely polygon in UTM coords and a row spacing variable
def rectGrid (polygon, rowSpace = 2):
    '''
    ----------------------------
    POLYGON GRID CREATION
    ----------------------------
    '''
    #grab coordinates from shapely object
    polyCoords = list(polygon.exterior.coords)

    #eliminate duplicate closing polygon point
    xy = np.array(polyCoords)
    xy = unique(xy)

    #find max y and second to max y to find top line of polygon, min to find side line
    maxy = xy[3]
    nmaxy = xy[2]

    #find delta x between two top points, if positive, maxy is top left, if negative, nmaxy is top left
    #moving from left -------> right always, use lastPtSide and lastPtTop/firstPtTop to calculate side line length
    dxtop = deltax(maxy,nmaxy)


    #find first point and top and side points from dxtop
    if dxtop > 0:
        firstPtTop = maxy
        lastPtTop = nmaxy
        lastPtSide = xy[1]
    if dxtop < 0:
        firstPtTop = nmaxy
        lastPtTop = maxy
        lastPtSide = xy[0]

    #calculate top and side of field lines to use for number of row and column cubes
    topField = LineString([firstPtTop,lastPtTop])
    sideField = LineString([firstPtTop,lastPtSide])
    rowCube = topField.length/rowSpace
    colCube = sideField.length/rowSpace

    #array of completed polys, structured as so:
    '''
    12345
    12345
    12345
    12345
    12345
    '''


    #calculate delta x and delta y for angle calculations, these angles don't change
    dxT = lastPtTop[0] - firstPtTop[0]
    dyT = lastPtTop[1] - firstPtTop[1]
    topLineAngle = math.atan2(dyT,dxT)

    dxS = lastPtSide[0] - firstPtTop[0]
    dyS = lastPtSide[1] - firstPtTop[1]
    sideLineAngle = math.atan2(dyS,dxS)

    angle1 = sideLineAngle
    angle2 = topLineAngle
    angle3 = math.pi + angle1

    firstPolyPt = firstPtTop
    polys = []

    #iterate through rows and columns and create polygons, column by column
    for j in range(1,int(math.ceil(rowCube)+1)):
            for i in range(1,int(math.ceil(colCube)+1)):
            #special case for first entry
                if i == 1:
                    #first polygon is special, needed for next column of polygons
                    firstPoly=[]
                    firstPoly.append(firstPolyPt)

                    #second point
                    x2 = firstPolyPt[0] + math.cos(angle1) * rowSpace
                    y2 = firstPolyPt[1] + math.sin(angle1) * rowSpace
                    secondPt = [x2,y2]
                    firstPoly.append(secondPt)

                    #third point
                    x3 = x2 + math.cos(angle2) * rowSpace
                    y3 = y2 + math.sin(angle2) * rowSpace
                    thirdPt = [x3,y3]
                    firstPoly.append(thirdPt)

                    #fourth point
                    x4 = x3 + math.cos(angle3) * rowSpace
                    y4 = y3 + math.sin(angle3) * rowSpace
                    fourthPt = [x4,y4]
                    firstPoly.append(fourthPt)

                    #create a shapely polygon and append to polys list, prepare to repeat
                    polys.append(Polygon(firstPoly))
                    oldPoly = firstPoly

                    #first polygon point for next column is fourth point of polygon
                    firstPolyPt = fourthPt

                else:
                    #first polygon point for next row is second point of previous polygon
                    firstPt  = oldPoly[1]

                    #second point
                    x2 = firstPt[0] + math.cos(angle1) * rowSpace
                    y2 = firstPt[1] + math.sin(angle1) * rowSpace
                    secondPt = [x2,y2]

                    #third point
                    x3 = x2 + math.cos(angle2) * rowSpace
                    y3 = y2 + math.sin(angle2) * rowSpace
                    thirdPt = [x3,y3]

                    #fourth point
                    x4 = x3 + math.cos(angle3) * rowSpace
                    y4 = y3 + math.sin(angle3) * rowSpace
                    fourthPt = [x4,y4]

                    nextPoly=[firstPt,secondPt,thirdPt,fourthPt]
                    polys.append(Polygon(nextPoly))
                    oldPoly = nextPoly
                    nextPoly = []
    return polys
    
    #create a rectangular grid of polygons given a shapely polygon in UTM coords and a row spacing variable
def irregGrid (boundary, rowLine, rowSpace = 2):
    '''
    ----------------------------
    POLYGON GRID CREATION
    ----------------------------
    '''
    #grab coordinates from shapely object
    bbox = boundary.bounds
    center = boundary.centroid
    
    pxmin = bbox[0]
    pymin = bbox[1]
    pxmax = bbox[2]
    pymax = bbox[3]
    
    linec = list(rowLine.coords)    
    
    firstPtTop = (pxmin,pymax)
    lastPtTop = (pxmax, pymax)
    lastPtSide = (pxmin, pymin)
    
    
    #calculate top and side of field lines to use for number of row and column cubes
    topField = LineString([firstPtTop,lastPtTop])
    sideField = LineString([firstPtTop,lastPtSide])
    rowCube = topField.length/rowSpace
    colCube = sideField.length/rowSpace
    
    
    #array of completed polys, structured as so:
    '''
    12345
    12345
    12345
    12345
    12345
    '''


    #calculate delta x and delta y for angle calculations, these angles don't change
    dxT = lastPtTop[0] - firstPtTop[0]
    dyT = lastPtTop[1] - firstPtTop[1]
    topLineAngle = math.atan2(dyT,dxT)

    dxS = lastPtSide[0] - firstPtTop[0]
    dyS = lastPtSide[1] - firstPtTop[1]
    sideLineAngle = math.atan2(dyS,dxS)

    angle1 = sideLineAngle
    angle2 = topLineAngle
    angle3 = math.pi + angle1

    firstPolyPt = firstPtTop
    polys = []

    firstPtLine = linec[0]
    lastPtLine  = linec[1]
    
    lineAngle = math.atan2((lastPtLine[1] - firstPtLine[1]),(lastPtLine[0] - firstPtLine[0]))
    
    
    #iterate through rows and columns and create polygons, column by column
    for j in range(1,int(math.ceil(rowCube)+1)):
            for i in range(1,int(math.ceil(colCube)+1)):
            #special case for first entry
                if i == 1:
                    #first polygon is special, needed for next column of polygons
                    firstPoly=[]
                    firstPoly.append(firstPolyPt)

                    #second point
                    x2 = firstPolyPt[0] + math.cos(angle1) * rowSpace
                    y2 = firstPolyPt[1] + math.sin(angle1) * rowSpace
                    secondPt = [x2,y2]
                    firstPoly.append(secondPt)

                    #third point
                    x3 = x2 + math.cos(angle2) * rowSpace
                    y3 = y2 + math.sin(angle2) * rowSpace
                    thirdPt = [x3,y3]
                    firstPoly.append(thirdPt)

                    #fourth point
                    x4 = x3 + math.cos(angle3) * rowSpace
                    y4 = y3 + math.sin(angle3) * rowSpace
                    fourthPt = [x4,y4]
                    firstPoly.append(fourthPt)

                    #create a shapely polygon, check if within boundary, if so rotate and add to list
                    poly = Polygon(firstPoly)
                    if poly.intersects(boundary):
                        poly = rotate(poly, lineAngle, origin = center, use_radians = True)
                        polys.append(poly) 
                    oldPoly = firstPoly

                    #first polygon point for next column is fourth point of polygon
                    firstPolyPt = fourthPt

                else:
                    #first polygon point for next row is second point of previous polygon
                    firstPt  = oldPoly[1]

                    #second point
                    x2 = firstPt[0] + math.cos(angle1) * rowSpace
                    y2 = firstPt[1] + math.sin(angle1) * rowSpace
                    secondPt = [x2,y2]

                    #third point
                    x3 = x2 + math.cos(angle2) * rowSpace
                    y3 = y2 + math.sin(angle2) * rowSpace
                    thirdPt = [x3,y3]

                    #fourth point
                    x4 = x3 + math.cos(angle3) * rowSpace
                    y4 = y3 + math.sin(angle3) * rowSpace
                    fourthPt = [x4,y4]
                    
                    #create shapely poly, add and rotate if in boundary
                    nextPoly=[firstPt,secondPt,thirdPt,fourthPt]
                    poly = Polygon(nextPoly)
                    if poly.intersects(boundary):
                        poly = rotate(poly, lineAngle, origin = center, use_radians = True)
                        polys.append(poly)

                    #prepare for next round
                    oldPoly = nextPoly
                    nextPoly = []
    return polys
    
# helper function to calculate point from relative polar coordinates (degrees)
def polarPoint(originPoint, angle,  distance):
    return [originPoint.x + math.sin(math.radians(angle)) * distance, originPoint.y + math.cos(math.radians(angle)) * distance]
    

def circGrid(spacing,length,center, boundary):
    
    # circle radius
    radius = spacing 
    # end of radius
    radiusEnd = int(math.ceil(length/radius))
    #width of sector in degrees
    sectorWidth = 4.0    
    polys = []     #array for storing features to plot
    
    for x in xrange(0,int(360.0/sectorWidth)):
        for r in xrange(1,radiusEnd + 1):
            if r==1:
                segmentVertices = []
                #first point is center
                centerPoint = polarPoint(center, 0,0)
                segmentVertices.append(centerPoint)
            
                #second point
                firstVertex = polarPoint(center, x*sectorWidth,r * radius)
                segmentVertices.append(firstVertex)
            
                #third point
                secondVertex = polarPoint(center, x * sectorWidth+sectorWidth, r * radius)
                segmentVertices.append(secondVertex)
            
                #center point ends polygon
                segmentVertices.append(centerPoint)
                oldVertices = segmentVertices            
                
                #add to polys if in field boundary
                poly = Polygon(segmentVertices)
                if poly.intersects(boundary):
                    polys.append(poly)                
            else:
                newVertices = []
                firstVertex = oldVertices[1]            
                newVertices.append(firstVertex)
                
                secondVertex = polarPoint(center, x*sectorWidth,r * radius)
                newVertices.append(secondVertex)
                
                thirdVertex = polarPoint(center, x * sectorWidth+sectorWidth, r * radius)
                newVertices.append(thirdVertex)
                
                fourthVertex = oldVertices[2]
                newVertices.append(fourthVertex) 
                
                #add to polys if in field boundary
                poly = Polygon(newVertices)
                if poly.intersects(boundary):
                    polys.append(poly)
                      
                oldVertices = newVertices
               
    return polys
    
#interpolate point dataset to gridded surface covering boundary of input polygon, input poly should be in UTM
def interpolateRaster(x,y,z,inP,outP,poly):
    
    #reproject x and y points into UTM
    points = reprojectArray(x,y,inP,outP)
    
    #get dimensions for grid
    bbox = poly.bounds
    xmin = bbox[0]
    ymin = bbox[1]
    xmax = bbox[2]
    ymax = bbox[3]
    
    #create x by y grid from data point min and maxes
    gridX, gridY = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    
    #interpolate to grid
    gridZ = griddata(points, z, (gridX, gridY), method = 'linear')
    
    #extrapolate beyond boundaries
    extrapolate_nans(gridX,gridY,gridZ)
    
    #return interpolated grid plus bounds for use in raster transformation
    return gridZ, [xmin,xmax,ymin,ymax]

#get stats that correspond to shapely polygon coverage of a numpy array ie raster
def getPolyStats(array,polys,bounds,stat = 'mean'):

    #get individual bounds variables
    xmin,xmax,ymin,ymax = bounds

    #get transform for raster stats
    nrows,ncols = np.shape(array)
    xres = (xmax-xmin)/float(ncols)
    yres = (ymax-ymin)/float(nrows)
    geotransform=(xmin,xres,0,ymin,0,yres)
    
    stats=[]
    for poly in polys:
        try:
            st = zonal_stats(poly, array.T, stats = stat, transform = geotransform)
            stats.append(st)
        except Exception:
            continue
    stats = np.array([i.get('mean') for j in stats for i in j])
    stats = stats.astype('float')
    boolm = np.isfinite(stats)
    polys = np.array(polys)[boolm]
    stats = stats[boolm]
    return stats, polys
    
    '''
    #calculate raster stats
    stats = zonal_stats(polys,array.T, stats=stat,transform = geotransform )

    #take list of dicts from stats and convert to numpy array of stat variable
    return np.array([0 if str(i.get('mean')) == 'nan' else i.get('mean') for i in stats])
    '''

def classify(value, breaks):
  for i in range(1, len(breaks)):
    if value < breaks[i]:
      return i
  return len(breaks) - 1
  

def plotMultip(polys):
    plt.close('all')
    fig = plt.figure()    
    ax = fig.add_subplot(111)
    if len(polys) < 7:
        colors = ['red','blue','green','orange','yellow','purple']
        for i,x in enumerate(polys):
            try:
                polyp = [PolygonPatch(p) for p in x]
                ax.add_collection(PatchCollection(polyp, facecolors = colors[i]))
            except TypeError:
                polyp = PolygonPatch(x, fc = colors[i])
                ax.add_patch(polyp) 
    for i,x in enumerate(polys):
            try:
                polyp = [PolygonPatch(p) for p in x]
                ax.add_collection(PatchCollection(polyp))
            except TypeError:
                polyp = PolygonPatch(x)
                ax.add_patch(polyp)    
    ax.autoscale()    
    plt.show()
    
def simplifyPolys(polys, spacing):
    #input list of polygons/multipolygons and spacing attribute used to produce polygons      
    #area
    area = math.floor(spacing**2)    
    #list of indexes for input polygons
    ilist = [i for i,x in enumerate(polys)]
    #iterate through array of multipolygons
    for idx,multip in enumerate(polys):
        #array for polygons to be removed from multip
        badp = []
        #list of indexes for input polygons - multip
        iilist = [i for i in ilist if (ilist[i] != idx)]
        #assume multip is a Multipolygon object
        try:   
            #iterate through multip    
            for ix,p in enumerate(multip):
                #if area less than desired
                if math.floor(p.area) <= area:         
                    #iterate through polys - multip                  
                    for ixx, multip2 in enumerate([polys[i] for i in iilist]):
                        #assume multip2 is multipolygon
                        try:
                            #for polygon in multip2
                            for p2 in multip2:
                                #if it intersects p, then append the index of p, reunion polygons and replace multip2 in polys
                                if p.intersects(p2):
                                    badp.append(ix)
                                    multip2 = cascaded_union([multip2,p])
                                    polys[iilist[ixx]] = multip2
                                    break
                            else:
                                # executed if the loop ended normally (no break)
                                continue  
                            # executed if 'continue' was skipped (break)
                            break  
                        #else it is a polygon
                        except TypeError:
                            if p.intersects(multip2):
                                badp.append(ix)
                                multip2 = cascaded_union([multip2,p])
                                polys[iilist[ixx]] = multip2
                                break
            #recreate list of polygons in multip, delete ones that were unioned with other multip's, and recreate multip
            multip = [p for p in multip]
            for index in sorted(badp, reverse=True):
                del multip[index]
            multip = cascaded_union(multip)
            polys[idx] = multip
        #except if multip is Polygon object
        except TypeError:
            p = multip
            #if area less than desired
            if math.floor(p.area) <= area:         
                #iterate through polys - multip                  
                for ixx, multip2 in enumerate([polys[i] for i in iilist]):
                    #assume multip2 is multipolygon
                    try:
                        #for polygon in multip2
                        for p2 in multip2:
                            #if it intersects p, then append the index of p, reunion polygons and replace multip2 in polys
                            if p.intersects(p2):
                                polys.remove(p)
                                multip2 = cascaded_union([multip2,p])
                                polys[iilist[ixx]] = multip2
                                break
                        else:
                            # executed if the loop ended normally (no break)
                            continue  
                        # executed if 'continue' was skipped (break)
                        break  
                    #else it is a polygon
                    except TypeError:
                        if p.intersects(multip2):
                            polys.remove(p)
                            multip2 = cascaded_union([multip2,p])
                            polys[iilist[ixx]] = multip2
                            break
    return [cascaded_union(p).simplify(.01) for p in polys if type(p) != GeometryCollection]

def gjsonJenks(polyStats, polys, inP, classes, spacing = None):
    #get the break points
    classes = jenks(polyStats, classes)

    #do the actual classification
    classified = np.array([classify(i,classes) for i in polyStats])

    #max value of zones
    maxz = max(classified)

    #nested list of zone indices
    zoneIndices = [[idx  for idx,val in enumerate(classified) if zone + 1 == val] for zone in range(maxz)]

    #nested list of polygons corresponding to each zone number
    polySort = [[polys[index] for index in zone] for zone in zoneIndices]

    #merge geometries, generate list of zones, create geojson feature collection from list
    polyComb = [cascaded_union(polyz).simplify(.01) for polyz in polySort]
    
    #if simplifying is needed        
    if spacing is not None:
        #simplify polys
        polyComb = simplifyPolys(polyComb, spacing)
    
    #reproject multipolygons to epsg:4326
    polyComb = reprojectMultiP(polyComb, inProjection = inP) 
        
    #create features and dump geojson
    features = [geojson.Feature(geometry=mapping(polyComb[zone]), id= zone, properties={"zone": zone+1}) for zone in range(len(polyComb))]
    
    crs = {
        "type": "name",
        "properties": {
            "name": "EPSG:4326"
            }}
    
    featureColl = geojson.FeatureCollection(features, crs = crs)
    return geojson.loads(geojson.dumps(featureColl))

def jsonColorCodes(numClasses):
    #special case for two zones
    if numClasses == 2:
        cmap = ['#FF0000', '#00FF00']
        return {key + 1: {'fillColor': value} for (key,value) in enumerate(cmap)}
    else:
        #get color map of blues based on number of classes    
        cmap = brewer2mpl.get_map('Blues', 'Sequential', numClasses)
        return {key + 1: {'fillColor': value} for (key,value) in enumerate(cmap.hex_colors)}

def asshape(geojson):
    return asShape(geojson)
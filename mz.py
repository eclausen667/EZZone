import numpy as np
from scipy.interpolate import griddata
from scipy.signal import medfilt
import math
import pyproj
from shapely.geometry import MultiPolygon, Polygon, asShape, LineString, mapping, MultiPoint
from shapely.ops import cascaded_union
from rasterstats import zonal_stats
import geojson
import random
from jenks import jenks

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

def boundJson(x,y, geometryType = 1, inProjection = 'epsg:4326' ,outProjection = 'epsg:3857'):
    """
    Create a GeoJSON boundary or point centroid from csv input
    Args:
        filename(file): the filename to perform the operation on
        geometryType(int): type of geometric shape, 1 = rectangular, 2 = circular, 3 = irregular
        inProjection(int): EPSG projection number for csv data
        outProjection(int): EPSG projection number for output GeoJSON string (default 3857 for google maps)

    Returns:
        GeoJSON string of polygon boundary for the csv data for geometry type 1 or 3,
        GeoJSON string of point centroid for the csv data for geometry type 2
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

    return geoJson #geojson.dumps(geoJson)


#orders y elements least to greatest of a 2d numpy array and eliminates duplicate values, returning the new array
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

#reproject geojson or shapely polygon and return a shapely polygon, specify geojson = True for geojson conversion, specify inPProj = True or outPProj = True to pass a pyproj.Proj object
def reproject(obj, inProjection = 'epsg:3857', outProjection = 'epsg:4326',gjson = False, inPProj = False, outPProj = False,multiP = False, gjsonString = False):
    if gjson == True:
        #convert to shapely object, easier to operate on coordinates this way
        poly = asShape(obj['geometry'])
    if gjsonString == True:
        obj = geojson.loads(obj)
        poly = asShape(obj['geometry'])
    else:
        poly = obj
    #pyproj transformation
    if inPProj == True:
        inP = inProjection
    else:
        inP = pyproj.Proj(init=inProjection)
    if outPProj == True:
        outP = outProjection
    else:
        outP = pyproj.Proj(init=outProjection)
    #if our input is a multipolygon
    if multiP == True:
        #array for new multipolygons
        multiPolys = []
        #iterate through list of multipolygons
        for polys in poly:
            #array for reprojected polygons
            newPolys = []
            #try to iterate through multipolygon, except when it is only a polygon
            try:
                for p in polys:
                    newCoords = [pyproj.transform(inP,outP,i[0],i[1]) for i in list(p.exterior.coords)]
                    poly = Polygon(newCoords)
                    newPolys.append(poly)
                #append multipolygon to list of multipolygons
                multiPolys.append(MultiPolygon(newPolys))
            except TypeError:
                newCoords = [pyproj.transform(inP,outP,i[0],i[1]) for i in list(polys.exterior.coords)]
                poly = Polygon(newCoords)
                multiPolys.append(poly)
        return multiPolys
    #transform coordinates and add to new array
    else:
        newCoords = [pyproj.transform(inP,outP,i[0],i[1]) for i in list(poly.exterior.coords)]
        poly = Polygon(newCoords)
        return poly

def utmzone(xy, shapelyPoly = False):
    if shapelyPoly == True:
        xy = list(xy.exterior.coords)[1]
    zone = (int(math.floor((xy[0] +180)/6) + 1), 1 if xy[1] < 0 else 0)
    if zone[1] == 1:
        return pyproj.Proj('+proj=utm +zone='+ str(zone[0]) + ', +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
    return pyproj.Proj('+proj=utm +zone='+ str(zone[0]) + ', +ellps=WGS84 +datum=WGS84 +units=m +no_defs')


#create a rectangular grid of polygons given a shapely polygon and a row spacing variable
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

#interpolate point dataset to gridded surface covering boundary of input polygon, input poly should be in UTM
def interpolateRaster(x,y,z,poly):

    #filter data values
    z= medfilt(z)

    #get min and max points of x and y polygon coords
    polyX = list(poly.exterior.coords.xy[0])
    polyY = list(poly.exterior.coords.xy[1])
    xmin=min(polyX)
    xmax=max(polyX)
    ymin=min(polyY)
    ymax=max(polyY)

    #create x by y grid from data point min and maxes
    gridX, gridY = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]

    #random subset of xyz, get sample, make array, transpose array
    rand = np.array(random.sample(zip(x,y,z),200)).T
    x = rand[0]
    y = rand[1]
    z = rand[2]
    #rand interpolate to grid
    gridZ = griddata((x,y), z, (gridX, gridY), method = 'linear')
    #extrapolate beyond boundaries
    extrapolate_nans(gridX,gridY,gridZ)
    return gridZ, [xmin,xmax,ymin,ymax]

#get stats that correspond to polygon coverage of a numpy array ie raster
def getPolyStats(array,polys,bounds,stat = 'mean'):

    #get individual bounds variables
    xmin,xmax,ymin,ymax = bounds

    #get transform for raster stats
    nrows,ncols = np.shape(array)
    xres = (xmax-xmin)/float(ncols)
    yres = (ymax-ymin)/float(nrows)
    geotransform=(xmin,xres,0,ymin,0,yres)

    #calculate raster stats
    stats = zonal_stats(polys,array.T, stats=stat,transform = geotransform)

    #take list of dicts from stats and convert to numpy array of stat variable
    return np.array([0 if np.isnan(i.get('mean')) else i.get('mean') for i in stats])

def classify(value, breaks):
  for i in range(1, len(breaks)):
    if value < breaks[i]:
      return i
  return len(breaks) - 1

def gjsonJenks(polyStats, polys, inProjection):
    #get the break points
    classes = jenks(polyStats, 5)

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
    #reproject multipolygons to epsg:4326
    polyComb = reproject(polyComb, inProjection = inProjection, inPProj = True, multiP = True)

    features = [geojson.Feature(geometry=mapping(polyComb[zone]), id= zone, properties={"zone": zone+1}) for zone in range(len(polyComb))]
    crs = {
        "type": "name",
        "properties": {
            "name": "EPSG:4326"
            }}
    featureColl = geojson.FeatureCollection(features, crs = crs)
    return featureColl


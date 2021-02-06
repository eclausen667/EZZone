'''
local app was changed to include a run.py and __init__.py file.
'''
from flask import Flask
# from app import app # added this for local app
from flask import render_template, request, send_file, Response
from werkzeug.utils import secure_filename
import mz
import numpy as np
import os
import tempfile

app = Flask(__name__)

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(APP_ROOT, 'static/uploads')
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['ALLOWED_EXTENSIONS'] = set(['txt', 'csv'])
#global variables
inP = 0
fieldType = 0
x = 0
y = 0
z = 0
rowSpacing = 0
pivotLength = 0
pivotSpacing = 0
classes = 0
fieldBound = 0
smooth = 0


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/upload', methods=['POST', 'GET'])
def upload():
    #global variables
    global inP, fieldType, x, y, z, classes, rowSpacing, pivotSpacing, pivotLength, fieldBound, smooth
    # get uploaded file, if ok then continue, else render badfile template
    fileup = request.files['file']
    if fileup and allowed_file(fileup.filename):
        filename = secure_filename(fileup.filename)
        filename = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        fileup.save(filename)
        # unpack values from file, else render badfile template
        try:
            x, y, z = np.loadtxt(filename, skiprows=1,
                                 delimiter=",", unpack=True)
        except ValueError:
            return render_template('badFile.html')
        # grab rest of values from form, modify to correct type
        classes = int(request.form['classes'])
        inP = request.form['coord']
        inP = 'epsg:' + str(inP)
        fieldType = request.form['field']
        fieldType = int(fieldType)
        # render templates according to field type
        if fieldType == 1:
            # check for row spacing value
            try:
                rowSpacing = float(request.form['row'])
                smooth = int(request.form['smooth'])
            except ValueError:
                return render_template('badParam.html')
            # check for valid epsg
            try:
                json = mz.boundJson(
                    x=x, y=y, geometryType=fieldType, inProjection=inP)
            except RuntimeError:
                return render_template('badEpsg.html')
            #render in browser
            return render_template('rect.html', string=json)
        if fieldType == 2:
            # check for valid
            try:
                pivotSpacing = float(request.form['spacing'])
            except ValueError:
                return render_template('badParam.html')
            try:
                pivotLength = float(request.form['length'])
            except ValueError:
                return render_template('badParam.html')
            try:
                json = mz.boundJson(
                    x=x, y=y, geometryType=fieldType, inProjection=inP)
                fieldBound = mz.asshape(json)
            except RuntimeError:
                return render_template('badEpsg.html')
            return render_template('vri.html', string=json)
        if fieldType == 3:
            # check for row spacing value
            try:
                rowSpacing = float(request.form['row'])
                smooth = int(request.form['smooth'])
            except ValueError:
                return render_template('badParam.html')
            # check for valid epsg
            try:
                json = mz.boundJson(
                    x=x, y=y, geometryType=fieldType, inProjection=inP)
            except RuntimeError:
                return render_template('badEpsg.html')
            #render in browser
            return render_template('irreg.html', string=json)
    else:
        return render_template('badFile.html')


@app.route('/createRect', methods=['POST'])
def createRect():
    # get the customized geojson, reproject to 4326, get utmzone, reproject to utm
    content = str(request.form['jsonval'])
    bound = mz.reprojectPoly(obj=content)
    utmPyProj = mz.utmzone(bound)
    bound = mz.reprojectPoly(
        obj=bound, inProjection='epsg:4326', outProjection=utmPyProj)
    # create grid, interpolate, get grid stats, create geojson
    rSpace = mz.getSpacing(bound, rowSpacing)
    polys = mz.rectGrid(bound, rSpace)
    array, bounds = mz.interpolateRaster(x, y, z, inP, utmPyProj, bound)
    polyStats, polys2 = mz.getPolyStats(array, polys, bounds)
    gjsonObj = mz.gjsonJenks(
        polyStats, polys2, utmPyProj, classes, (rSpace*smooth))
    # create colormap
    cmap = mz.jsonColorCodes(classes)
    #render in browser
    return render_template('displayZones.html', string=gjsonObj, cmap=cmap)


@app.route('/createVri', methods=['POST'])
def createVri():
    # get the customized geojson, reproject to 4326, get utmzone, reproject to utm
    content = str(request.form['jsonval'])
    bound = mz.geojFCparser(content, 'Polygon')
    point = mz.geojFCparser(content, 'Point')
    bound = mz.reprojectPoly(bound)
    utmPyProj = mz.utmzone(bound)
    bound = mz.reprojectPoly(
        bound, inProjection='epsg:4326', outProjection=utmPyProj)
    point = mz.reprojectPoint(point)
    point = mz.reprojectPoint(
        obj=point, inProjection='epsg:4326', outProjection=utmPyProj)
    # create grid, interpolate, get grid stats, create geojson
    polys = mz.circGrid(pivotSpacing, pivotLength, point, bound)
    array, bounds = mz.interpolateRaster(x, y, z, inP, utmPyProj, bound)
    polyStats, polys2 = mz.getPolyStats(array, polys, bounds)
    gjsonObj = mz.gjsonJenks(polyStats, polys2, utmPyProj, classes)
    # create colormap
    cmap = mz.jsonColorCodes(classes)
    #render in browser
    return render_template('displayZones.html', string=gjsonObj, cmap=cmap)


@app.route('/createIrreg', methods=['POST'])
def createIrreg():
    # get the customized geojson, reproject to 4326, get utmzone, reproject to utm
    content = str(request.form['jsonval'])
    bound = mz.geojFCparser(content, 'Polygon')
    line = mz.geojFCparser(content, 'LineString')
    bound = mz.reprojectPoly(bound)
    line = mz.reprojectLine(line)
    utmPyProj = mz.utmzone(bound)
    bound = mz.reprojectPoly(
        bound, inProjection='epsg:4326', outProjection=utmPyProj)
    line = mz.reprojectLine(
        line, inProjection='epsg:4326', outProjection=utmPyProj)
    # create grid, interpolate, get grid stats, create geojson
    rSpacing = mz.getSpacing(bound, rowSpacing)
    polys = mz.irregGrid(bound, line, rSpacing)
    array, bounds = mz.interpolateRaster(x, y, z, inP, utmPyProj, bound)
    polyStats, polys2 = mz.getPolyStats(array, polys, bounds)
    gjsonObj = mz.gjsonJenks(
        polyStats, polys2, utmPyProj, classes, (rSpacing*smooth))
    # create colormap
    cmap = mz.jsonColorCodes(classes)
    #render in browser
    return render_template('displayZones.html', string=gjsonObj, cmap=cmap)


@app.route('/sendFile', methods=['POST'])
def sendFile():
    content = str(request.form['jsonval'])
    return Response(content,
                    mimetype='application/json',
                    headers={'Content-Disposition': 'attachment;filename=zones.geojson'})


app.secret_key = '12ad432gfd'
if __name__ == '__main__':
    # main()
    port = int(os.environ.get('PORT', 5000))
    app.run(host='127.0.0.1', port=port)

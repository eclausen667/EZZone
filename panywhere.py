from flask import render_template, request, Flask
from werkzeug import secure_filename
import mz
import numpy as np
import os

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

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST','GET'])
def upload():
    global inP, fieldType, x, y, z
    fileup = request.files['file']
    if fileup and allowed_file(fileup.filename):
        filename = secure_filename(fileup.filename)
        filename = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        fileup.save(filename)
        global inP, fieldType, x, y, z
        x,y,z = np.loadtxt(filename,skiprows=1,delimiter=",", unpack=True)
        inP = request.form['coord']
        inP = 'epsg:' + str(inP)
        fieldType = request.form['field']
        fieldType = int(fieldType)
        json = mz.boundJson(x = x, y = y,geometryType = fieldType, inProjection =inP)
        return render_template('rect.html',string=json)

@app.route('/create', methods=['POST'])
def create():
    content = request.form['jsonval']
    poly = mz.reproject(obj = content, gjsonString = True)
    utmPyProj = mz.utmzone(poly, shapelyPoly = True)
    poly = mz.reproject(obj = poly, inProjection = 'epsg:4326', outProjection = utmPyProj, outPProj = True)
    polys = mz.rectGrid(poly, rowSpace = 2)
    array, bounds = mz.interpolateRaster(x,y,z,poly)
    polyStats = mz.getPolyStats(array,polys,bounds)
    gjsonObj = mz.gjsonJenks(polyStats,polys,utmPyProj)
    return render_template('rectdone.html', string = gjsonObj)
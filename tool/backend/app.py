#!/usr/bin/env python

##################################################################################
# app.py
from flask import Blueprint
main = Blueprint('main', __name__)

import urllib
import json
from engine import ICS5200Engine

import logging
logging.basicConfig(filename = "~/ics5200.log", level=logging.INFO)
logger = logging.getLogger(__name__)

from flask import Flask, request
from flask_cors import CORS, cross_origin

from decimal import Decimal

from cheminfo import ChemInfo
from moleculehelper import *
from pythonhelper import *

################################################################################################################################################################
# Helper class
class CustomEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, Decimal):
            return float(o)
        return super(CustomEncoder, self).default(o)

################################################################################################################################################################
# Ligands
@main.route("/getTestLigandsDS/", methods=["GET"])
def getTestLigandsDS():
    return json.dumps(ics5200Engine.getTestLigandsDS(), cls=CustomEncoder)

@main.route("/getSmilesSVG/<input>/<inputType>", methods=["GET"])
def getSmilesSVG(input, inputType):    
    if inputType.lower() == "smiles":
        smiles = urllib.unquote(input)
    else:
        smiles = ics5200Engine.getSmiles(input)
    if len(smiles) > 0:
        PythonHelper.writeToJupyterConsole(">smilesToSVG: " + smiles)
        logger.debug("smilesToSVG: " + smiles)
        return ChemInfo.smilesToSVG(smiles)
    else:
        return ""

@main.route("/getSmiles/<molRegNo>", methods=["GET"])
def getSmiles(molRegNo):    
    return json.dumps(ics5200Engine.getSmiles(molRegNo))
        
@main.route("/getLigandBindings/<smiles>", methods=["GET"])
def getLigandBindings(smiles):
    return json.dumps(ics5200Engine.getLigandBindings(smiles))

@main.route("/doLigandExperiment/<molRegNo>/<fingerprint>/<similarity>/<threshold>", methods=["GET"])
def doLigandExperiment(molRegNo, fingerprint, similarity, threshold):
    return json.dumps(ics5200Engine.doLigandExperiment(molRegNo, LigandHelper, fingerprint, similarity, float(threshold)), cls=CustomEncoder)

@main.route("/doLigandExperimentFromSmiles/<smiles>/<fingerprint>/<similarity>/<threshold>", methods=["GET"])
def doLigandExperimentFromSmiles(smiles, fingerprint, similarity, threshold):
    return json.dumps(ics5200Engine.doLigandExperimentFromSmiles(urllib.unquote(smiles), LigandHelper, fingerprint, similarity, float(threshold)), cls=CustomEncoder)

@main.route("/isLigandInChEMBL/<smiles>", methods=["GET"])
def isLigandInChEMBL(smiles):
    return json.dumps(ics5200Engine.isLigandInChEMBL(urllib.unquote(smiles)))

################################################################################################################################################################
# Proteins
@main.route("/isProteinInChEMBL/<sequence>", methods=["GET"])
def isProteinInChEMBL(sequence):
    return json.dumps(ics5200Engine.isProteinInChEMBL(sequence), cls=CustomEncoder)

@main.route("/doProteinExperiment/<compId>", methods=["GET"])
def doProteinExperiment(compId):
    return json.dumps(ics5200Engine.doProteinExperiment(compId), cls=CustomEncoder)

################################################################################################################################################################
# Creating Flask App
def createApp(sparkContext, datasetPath):
    global ics5200Engine

    ics5200Engine = ICS5200Engine(sparkContext, datasetPath)

    app = Flask(__name__)
    CORS(app)
    app.register_blueprint(main)
    return app


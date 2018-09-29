#!/usr/bin/env python

##################################################################################
# server.py
import time, sys, cherrypy, os
from paste.translogger import TransLogger
from app import createApp
from pyspark import SparkContext, SparkConf

def initSparkContext():
    # load spark context
    conf = SparkConf().setAppName("ics5200-server")
    # IMPORTANT: pass aditional Python modules to each worker
    sc = SparkContext(conf=conf, pyFiles=['moleculehelper.py', 'pythonhelper.py'])
    sc.setLogLevel("INFO")

    return sc

def runServer(app):

    # Enable WSGI access logging via Paste
    appLogged = TransLogger(app)

    # Mount the WSGI callable object (app) on the root directory
    cherrypy.tree.graft(appLogged, '/')

    # Set the configuration of the web server
    cherrypy.config.update({
        'engine.autoreload.on': True,
        'log.screen': True,
        'server.socket_port': 5432,
        'server.socket_host': 'hadoop1'
    })

    # Start the CherryPy WSGI web server
    cherrypy.engine.start()
    cherrypy.engine.block()


if __name__ == "__main__":
    # Init spark context and load libraries
    sc = initSparkContext()
    dataset_path = '' # os.path.join('datasets', 'ml-latest')
    app = createApp(sc, dataset_path)

    # start web server
    runServer(app)


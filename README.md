SpectraViewer
=============

[![Build Status](https://travis-ci.org/chrismit/SpectraViewer.svg?branch=master)](https://travis-ci.org/chrismit/SpectraViewer)

A cross-platform viewer to open up mass spec peaklist files.

To use, simply drag and drop a file onto the program.  Multiple tabs are provided to allow for any number of files to be open concurrently.
Dragging multiple files into a window will open all them at once for a combined view.
Holding Control+dropping a file will append it to the current view (this way, you can even mix & match file types).

The current goal is functionality -- elegance can come later.

Installation:
You need the following installed:
PyQt4 -- install instructions are here: http://pyqt.sourceforge.net/Docs/PyQt4/
pyqtgraph

The rest can be done with:
pip install -r requirements.txt

Supported:
mgf files - on first load the mgf file will be indexed for random access.  This is a one
	time cost, but it is worth it.

X!Tandem xml files

Proteome Discovered .msf files - These load very fast, anyone used to PD's abysmal (and unthreaded) speed will appreciate this

MZML/PepXML

Mascot DAT Files

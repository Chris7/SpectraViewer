SpectraViewer
=============

A cross-platform viewer to open up mass spec peaklist files.

To use, simply drag and drop a file onto the program.  Multiple tabs are provided to allow for any number of files to be open concurrently.

The current goal is functionality -- elegance can come later.

Supported:
mgf files - on first load the mgf file will be indexed for random access.  This is a one
	time cost, but it is very worth it.
X!Tandem xml files - X!Tandem has a strange way of storing the m/z and peptide matches, 
	it appears to treat everything as a singly charged ion (so you will see a spectra with
	2+ charge having the peptides matched as if they are singly charged.  I have no idea 
	if this is buggy behavior at the moment). 
GFF files - this is a feature that is currently, and probably will remain only used by me.  
	I translate my peptide results to a GFF file which I view in GBrowse, this allows me 
	to see peptides covering exon junctions, etc. nicely.  If you add a GFF attribute of 
	spectraId1 to your GFF files (such as A1.xxxx.xxx.x.dta) and spectraId2 (A1.xxx.xx.x)
	this feature will work for you.  Set the search path to the mgf files via settings -
	add search path, a config file will remember what dirs to search in the future).


Soon to be added:
dta files
mzXML files
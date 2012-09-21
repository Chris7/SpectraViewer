SpectraViewer
=============

A cross-platform viewer to open up mass spec peaklist files.

To use, simply drag and drop a file onto the program.  Multiple tabs are provided to allow for any number of files to be open concurrently.
Dragging multiple files into a window will open all them at once for a combined view (though they need to be the same file type at the moment.)
Holding Control+dropping a file will append it to the current view (this way, you can even mix & match file types).

The current goal is functionality -- elegance can come later.

Supported:
mgf files - on first load the mgf file will be indexed for random access.  This is a one
	time cost, but it is worth it.
	
X!Tandem xml files

Proteome Discovered .msf files - These load very fast, anyone used to PD's abysmal (and unthreaded) speed will appreciate this

GFF files - this is a feature that is currently, and probably will remain only used by me.  
	I translate my peptide results to a GFF file which I view in GBrowse, this allows me 
	to see peptides covering exon junctions, etc. nicely.  If you add a GFF attribute of 
	spectraId1 to your GFF files (such as A1.xxxx.xxx.x.dta) and spectraId2 (A1.xxx.xx.x)
	this feature will work for you.  Set the search path to the mgf files via settings -
	add search path, a config file will remember what dirs to search in the future).


Soon to be added:

dta files
mzXML files
some simple higher level analysis tools like box plots & venn diagrams
selected spectra export
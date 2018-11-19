# IDL-read-write-Nifti
### reading/writing Nifit files using interactive data language (IDL)

### The reader/writer can handle the following formats:###
gz; nii, and hdr+img pairs.

#### Some remarks 
With the three files above you should be able to read/write Nifit data. The software has been mainly tested on a Siemens Trio, but should also work for other scanner types if a conversion to Nifti standard was applied. However, the code may not handle all kinds of scanner specific variations (or format conversions thereafter). in addition, you may want to change some default parameters; e.g., in our institution we applied a vertically flip to our data, but most of this can be controlled by using keywords.

Meanwhile, I have switched over to python and I have not used the reader/writer for several years. The last modifications were applied in 2014.

### plot_dscam_v3.py  

This script was authored by W. Xu and adapted from a script originally written by R. Arnold.  

Language:  
Python  

Inputs:  
.gtf file - This is a genome annotation file.  
sj.out.tab file(s) - This is a splice junction file, produced from STAR aligner.  

Purpose:  
To draw support levels for each exon.  
  - height of bar represents read depth for splices at that junction.  
  - green is canonical donor site.  
  - yellow is canonical acceptor site.  
  - blue represents position of exon as defined in the annotation file.  

Output:  
.svg - A vector graphic



### batch_resize_svg

This is a bash script for resizing all .svg files in a folder.
  - Execute this script from the folder which contains all the svgs you want to resize
  - Change the -w and -h accordingly  

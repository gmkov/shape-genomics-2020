# From pictures to wing data- Fiji

This tutorial allows extracting wing shape and size variables from images taken in standard conditions

## Pre-requisites

* Download Fiji (https://imagej.net/Fiji/Downloads)
* Download images associated to this paper, two options:
    * from [Earthcape](https://heliconius.ecdb.io/), by searching for the Unit IDs (individual ids) on the menu Sample Images > Filter List, then paste list of IDs as such “CAM016112_d.JPG” to search. 
    * from zenodo directly with the links provided in the data section. you can use ‘wget -i’ in the command line to download directly.

## Script overview

1) set a scale so that the measurements it takes represent real sizes: 
make a line of 5 mm on the ruler → Analyze → Set Scale → change “Known distance” to 5 mm

2) crop the picture up and grab the "rectangle"/wing you want to study
make a rectangle around the wing → Image → Crop

3) edit the image to enhance the wing edges 
Image → Adjust → change as required
Process → Enhance contrast (saturated=0)

4) make the image binary (background white, wing black) 
Process → Binary → Options → Set to 25 iterations, Count 1, pad edges when eroding
Process → Binary → Make binary 
Process → Binary → Fill holes

5) carry out particle analyses (measure the biggest object -wing- size, shape, length etc., create an outline of this object) 
Analyze → Analyze Particles → Set Size to 200-Infinity, Show Outlines, tick: display results, summarize, exclude on edges

6) store results in a csv, export an outline tiff/jpeg of the object(wing) that has been measured (for post-fiji checks)

## Automatization with Macros:

These steps are written as a macro batch script, Fiji will run the same set of commands on all pictures within a folder (you determine input/output folders). If wanting to reproduce or reuse these scripts first it is recommended to test and understand how macros work.

### 1. Learn how to write a Macro

This step is useful, to know how to script a macro and then apply it to a batch of photos, is good to “see” how fiji codes

Open fiji and load an image via File → Open

Go to **Plugins → Macros →  Record**
This will record what you do. Open an image (one from input 1) with Fiji

Repeat your favourite steps as above and see how the code for the Macro is written

### 2. Test your macro script

Open Fiji
Click **Process > Batch > Macro**
Select Input directory (where the images are)
Select output directory
Output format JPEG

Then paste your script (edited from the one below according to the script you wrote in step 1)

Click Test (to test it!) - it'll do one image
When you think things are looking good and the test worked: first COPY AND PASTE your edited version of the script onto a text file, then click PROCESS (and itll do all the images within input1)
Voila!

_These are the scripts used for most images associated to the paper, slight alterations depending on background and light conditions_

```
# grey batch
dir="/Users/XXX/size.batch.haplotagging/output1/";
name = getTitle; 
mainTitle=getTitle();
makeLine(1749, 3250, 2004, 3248);
run("Set Scale...", "distance=255.0078 known=5 pixel=1 unit=mm global");
makeRectangle(2596, 4, 2588, 1880);
run("Crop");
run("Duplicate...", "title=originalimage");
selectWindow(mainTitle);
run("Options...", "iterations=5 count=1 pad");
run("Enhance Contrast...", "saturated=3");
setMinAndMax(0, 210);
run("RGB Stack");
run("Convert Stack to Images")
selectWindow("Blue");
setThreshold(0, 115);
rename(name);
run("Convert to Mask", "method=Default background=Light");
run("Fill Holes", "stack");
run("Fill Holes", "stack");
run("Set Measurements...", "area fit shape feret's limit display add redirect=None decimal=3");
run("Analyze Particles...", "size=200-Infinity show=Ellipses display exclude summarize record in_situ add");
run("Overlay Options...", "stroke=white width=100 fill=none set");
roiManager("Set Color", "white");
roiManager("Set Line Width", 7);
selectWindow("originalimage");
roiManager("Show All without labels"); 
run("Flatten");
saveAs("Results",  dir+ "results.csv");



```


Images with green backgrounds

```
dir="/Users/XXXX/output8/";
name = getTitle; 
mainTitle=getTitle();
makeLine(1749, 3250, 2004, 3248);
run("Set Scale...", "distance=190.0026 known=5 pixel=1 unit=mm global");
makeRectangle(2706, 42, 2436, 1578);
run("Crop");
run("Duplicate...", "title=originalimage");
selectWindow(mainTitle);
run("Options...", "iterations=5 count=1 pad");
run("Brightness/Contrast...");
setMinAndMax(23, 179);
run("Subtract Background...", "rolling=50 light sliding disable");
run("Gaussian Blur...", "sigma=2");
run("RGB Stack");
run("Convert Stack to Images")
selectWindow("Green");
setThreshold(0, 115);
rename(name);
run("Convert to Mask", "method=Default background=Light");
run("Fill Holes", "stack");
run("Fill Holes", "stack");
run("Set Measurements...", "area fit shape feret's limit display add redirect=None decimal=3");
run("Analyze Particles...", "size=200-Infinity show=Ellipses display exclude summarize record in_situ add");
run("Overlay Options...", "stroke=white width=100 fill=none set");
roiManager("Set Color", "white");
roiManager("Set Line Width", 7);
selectWindow("originalimage");
roiManager("Show All without labels"); 
run("Flatten");
saveAs("Results",  dir+ "results.csv");
```


**Tips: **things to edit depending on the wings you have are the saturation levels (in enhance contrast line), the iterations and counts and anything you can image to do to help fiji find a wing (and not the background) when making the image binary. 

Fiji will run these steps for each picture in 1-5s (without any user interaction), so in a couple of minutes you can get thousands of wing measurements 

Outputs can be checked visually (the white outline is the particle that got analysed by the script): check pdf version 
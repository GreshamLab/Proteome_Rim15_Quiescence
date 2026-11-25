//Select Directories
input = getDirectory("Choose Source Directory");
output = getDirectory("Choose Destination Directory");

//===============================

run("Bio-Formats Macro Extensions"); 									// need this to open images with the special importer that does not have the popup window

list = getFileList(input);
for (i = 0; i < list.length; i++)
        action(input, output, list[i]);

function action(input, output, filename) {
        Ext.openImagePlus(input+list[i]); 								// to open images without the bioformats importer
        run("Z Project...", "[Max Intensity]"); 				// does the max projection
        run("Grays"); 													// does gray LUT
        saveAs("Tiff", output + "Max_" + filename);					// save as TIFF into the "output" folder and adds "_max-pro" to the end of the file name
        run("Close All");
}

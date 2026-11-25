//Select Directories
input = getDirectory("Choose Source Directory");
output = getDirectory("Choose Destination Directory");

//===============================

run("Bio-Formats Macro Extensions"); // need this to open images with the special importer that does not have the popup window

list = getFileList(input);
for (i = 0; i < list.length; i++)
        action(input, output, list[i]);

function action(input, output, filename) {
        Ext.openImagePlus(input+list[i]); //to open images without the bioformats importer
        makeRectangle(100, 100, 800, 800);
        run("Crop");
        run("Grays");
        saveAs("Tiff", output + filename + "_crop");
        close();
}

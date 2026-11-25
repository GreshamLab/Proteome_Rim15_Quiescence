//Select Directories
dir1 = getDirectory("Choose Source Directory ");
dir2 = getDirectory("Choose Destination Directory Results ");
//===============================


run("Bio-Formats Macro Extensions");

//Looping function through the files in a directory

list = getFileList(dir1);
for (i=0; i<list.length; i++) 
{
	if (endsWith(list[i], ".tif")) 
	{ 
		Ext.openImagePlus(dir1+list[i]);
		SingleParticleTracking();
	}
}

//all files set to run
//need to check for when files are finished - rename them as needed
for (i=0; i<list.length; i++) 
{
	if (endsWith(list[i], ".tif")) 
	{
		
		origFile = dir1 + "Traj_" + list[i] + ".csv";
		newFile = dir2 + "Traj_" + replace(list[i], ".tif", ".csv");
		while(!File.exists(origFile))
		{
			wait(100);
		}
		wait(100);
		s=File.rename(origFile, newFile);	/// IF FILE ALREADY EXISTS IT WON'T OVERWRITE!  TRY TO FIX ! 
		
		origFile = dir1 + "Traj_" + list[i] + ".txt";
		newFile = dir2 + "Traj_" + replace(list[i], ".tif", ".txt");
		// s=File.rename(origFile, newFile);	/// IF FILE ALREADY EXISTS IT WON'T OVERWRITE!  TRY TO FIX ! 
		
	}
}

function SingleParticleTracking() 
{
	//waitForUser("W1");	
	run("Particle Tracker 2D/3D", "radius=3 cutoff=0 per/abs=0.5 link=1 displacement=5 dynamics=Brownian");
	//waitForUser("W2");
	run("Close All");
}

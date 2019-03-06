/*
Cheetah for EIGER
 Written by Takanori Nakane 

Don't forget to set HDF5_PLUGIN_PATH!
*/

//  Distantly based on cheetah-sacla, cheetah-eiger from Takanory Nakane and cheetah-euxfel
//   Created by Anton Barty on 20/1/14.
//   Copyright (c) 2014 Anton Barty. All rights reserved.

#include <iostream>
#include <stdlib.h>
#include <getopt.h>    /* for getopt_long; POSIX standard getopt is in unistd.h */
#include <string>
#include <stdio.h>
#include <math.h>

#ifdef __WITH_MPI__
#include <mpi.h>
#endif /* __WITH_MPI__ */

#include <hdf5.h>
#include <hdf5_hl.h>

#include "cheetah.h"

// This is for parsing getopt_long()
struct tCheetahEigerparams {
    const char *masterFile;
    const char *iniFile;
    int frameFirst;
    int frameStep;
    int frameLast;
} CheetahEigerparams;

void parse_config(int, char *[], tCheetahEigerparams*);
int cheetah_process_file(tCheetahEigerparams*);

int main(int argc, char * argv[]) {

  int ret;

  #ifdef __WITH_MPI__
  int my_rank, num_procs;
  char hname[1024];
  int lenhname;

  int irequired, iprovided;
  irequired = MPI_THREAD_MULTIPLE;
  MPI_Init_thread(&argc, &argv, irequired, &iprovided);

  /* Identify this process */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Get_processor_name(hname, &lenhname);
    
  /* Find out how many total processes are active */
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  /* Until this point, all programs have been doing exactly the same.
     Here, we check the rank to distinguish the roles of the programs */
  if (my_rank == 0) {
    printf("(MPI_Init_thread) required: %d, provided: %d\n", irequired, iprovided);
    printf("We have %i processes.\n", num_procs);
  #endif /* __WITH_MPI__ */ 

        printf("Cheetah for EIGER\n");
	printf(" by Takanori Nakane\n");
	printf(" based on cheetah-sacla by Anton Barty\n");
	printf("\nIf this program fails, make sure HDF5_PLUGIN_PATH points to lz4 plugin!\n");

   #ifdef __WITH_MPI__
   }
   #endif /* __WITH_MPI__ */

	// Parse configurations
	parse_config(argc, argv, &CheetahEigerparams);

	std::string filename(CheetahEigerparams.masterFile), cheetahini(CheetahEigerparams.iniFile);
	
   #ifdef __WITH_MPI__
   if (my_rank == 0) {
   #endif /* __WITH_MPI__ */
   
    printf("Program name: %s\n",argv[0]);
    printf("Input data file: %s\n", filename.c_str());
    printf("Cheetah .ini file: %s\n", cheetahini.c_str());
	
	printf("first: %d\n", CheetahEigerparams.frameFirst);
	printf("step: %d\n", CheetahEigerparams.frameStep);
	printf("last: %d\n", CheetahEigerparams.frameLast);

    #ifdef __WITH_MPI__
    }
    // define specific first frame for each rank
    CheetahEigerparams.frameFirst += my_rank*CheetahEigerparams.frameStep;
    CheetahEigerparams.frameStep *= num_procs;
    #endif /* __WITH_MPI__ */

	// Process file
    ret = cheetah_process_file(&CheetahEigerparams);

    #ifdef __WITH_MPI__
    /* Tear down the communication infrastructure */
    MPI_Finalize();
    #endif /* __WITH_MPI__ */

    return 0;
}

int cheetah_process_file(tCheetahEigerparams *global) {

	// Parse metadata
	printf("\n** Parsing metadata in the HDF5 file **\n");

	int xpixels = -1, ypixels = -1, beamx = -1, beamy = -1, ival = -1, nimages = -1;
	double pixelsize = -1, wavelength = -1, distance = -1;
	double count_time = -1, frame_time = -1, osc_width = -1;
	char detector_sn[256] = {}, description[256] = {};

	hid_t hdf;	

	hdf = H5Fopen(global->masterFile, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (hdf < 0) {
		fprintf(stderr, "Failed to open file %s\n", global->masterFile);
		return -1;
	}
	printf("Successfully opened %s\n", global->masterFile);

	H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/nimages", &ival);
    nimages = ival;
	printf("Found nimages=%d in the input\n", ival);
	H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/ntrigger", &ival);
	printf("Found ntrigger=%d in the input\n", ival);
	nimages *= ival;
	printf("Asssuming total number of %d images in the input\n", nimages);

	if (global->frameLast<0)
		global->frameLast = nimages-1;

	H5Eset_auto(0, NULL, NULL); // Comment out this line for debugging.

	fprintf(stderr, "Metadata in HDF5:\n");
	H5LTread_dataset_string(hdf, "/entry/instrument/detector/description", description);
	fprintf(stderr, " /entry/instrument/detector/description = %s\n", description);
	H5LTread_dataset_string(hdf, "/entry/instrument/detector/detector_number", detector_sn);
	fprintf(stderr, " /entry/instrument/detector/detector_number = %s\n", detector_sn);
	H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/x_pixels_in_detector", &xpixels);
	H5LTread_dataset_int(hdf, "/entry/instrument/detector/detectorSpecific/y_pixels_in_detector", &ypixels);
	fprintf(stderr, " /entry/instrument/detector/detectorSpecific/{x,y}_pixels_in_detector = (%d, %d) (px)\n",
			xpixels, ypixels);
	H5LTread_dataset_int(hdf, "/entry/instrument/detector/beam_center_x", &beamx);
	H5LTread_dataset_int(hdf, "/entry/instrument/detector/beam_center_y", &beamy);
	fprintf(stderr, " /entry/instrument/detector/beam_center_{x,y} = (%d, %d) (px)\n", beamx, beamy);
	H5LTread_dataset_double(hdf, "/entry/instrument/detector/count_time", &count_time); // in m
	fprintf(stderr, " /entry/instrument/detector/count_time = %f (sec)\n", count_time);
	H5LTread_dataset_double(hdf, "/entry/instrument/detector/frame_time", &frame_time); // in 
	fprintf(stderr, " /entry/instrument/detector/frame_time = %f (sec)\n", frame_time);
	H5LTread_dataset_double(hdf, "/entry/instrument/detector/x_pixel_size", &pixelsize); // in 
	H5LTread_dataset_double(hdf, "/entry/instrument/detector/detector_distance", &distance);
	fprintf(stderr, " /entry/instrument/detector/detector_distance = %f (m)\n", distance);
	fprintf(stderr, " /entry/instrument/detector/x_pixel_size = %f (m)\n", pixelsize);

	// There are many ways to store wavelength!
	H5LTread_dataset_double(hdf, "/entry/instrument/beam/wavelength", &wavelength);
	if (wavelength > 0) {
		fprintf(stderr, " /entry/instrument/beam/wavelength = %f (A)\n", wavelength);
	} else {
		fprintf(stderr, "  /entry/instrument/beam/wavelength not present. Trying another place.\n");
		H5LTread_dataset_double(hdf, "/entry/instrument/monochromator/wavelength", &wavelength);
		if (wavelength > 0) {
			fprintf(stderr, " /entry/instrument/monochromator/wavelength = %f (m)\n", wavelength);
		} else {
			fprintf(stderr, "  /entry/instrument/monochromator/wavelength not present. Trying another place.\n");
			H5LTread_dataset_double(hdf, "/entry/instrument/beam/incident_wavelength", &wavelength);
			if (wavelength > 0) {
				fprintf(stderr, " /entry/instrument/beam/incident_wavelength = %f (m)\n", wavelength);
			}
		}
	}
	if (wavelength < 0) {
		fprintf(stderr, " wavelength not defined!");
	}

	// Same for oscillation range (but not relevant here)
	H5LTread_dataset_double(hdf, "/entry/sample/goniometer/omega_range_average", &osc_width);
	if (osc_width > 0) {
		fprintf(stderr, " /entry/sample/goniometer/omega_range_average = %f (deg)", osc_width);
	} else {
		fprintf(stderr, " oscillation width not defined. \"Start_angle\" is set to 0!\n");
		osc_width = 0;
	}

	
	hid_t entry, group;
	entry = H5Gopen2(hdf, "/entry", H5P_DEFAULT);
	if (entry < 0) {
		fprintf(stderr, "/entry does not exist!\n");
		return -1;
	}
	
	// Check if /entry/data present
	group = H5Gopen2(entry, "data", H5P_DEFAULT);  
	if (group < 0) {
		group = entry; // FIXME: leak!
	}

	int start_block_number = 1;
	if (H5LTfind_dataset(group, "data_000000")) {
		fprintf(stderr, "This dataset starts from data_000000.\n");
		start_block_number = 0;
	} else {
		fprintf(stderr, "This dataset starts from data_000001.\n");
	}

	char data_name[20] = {};
	hid_t data, dataspace;
	int number_per_block = 0;
	
	// Open the first data block to get the number of frames in a block
	snprintf(data_name, 20, "data_%06d", start_block_number); 
	data = H5Dopen2(group, data_name, H5P_DEFAULT);
	dataspace = H5Dget_space(data);
	if (data < 0) {
		fprintf(stderr, "failed to open /entry/%s\n", data_name);
		return -1;
	}
	if (H5Sget_simple_extent_ndims(dataspace) != 3) {
		fprintf(stderr, "Dimension of /entry/%s is not 3!\n", data_name);
		return -1;    
	}
	hsize_t dims[3];
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	number_per_block = dims[0];
	fprintf(stderr, "The number of images per data block is %d.\n", number_per_block);

	H5Sclose(dataspace);
	H5Dclose(data);
	
	//	Initialise Cheetah
	printf("** Setting up Cheetah **\n");
	static uint32_t ntriggers = 0;
	long frameNumber = 0;
    long runNumber = 0;
	static cGlobal cheetahGlobal;
	signed short *buf = (signed short*)malloc(sizeof(signed short) * xpixels * ypixels);

	static time_t startT = 0;
	time(&startT);
    strcpy(cheetahGlobal.configFile, global->iniFile);
	cheetahInit(&cheetahGlobal);

    for (frameNumber = global->frameFirst; frameNumber <= global->frameLast; ) {
        printf("Processing frame %ld\n", frameNumber);
		
		//  Cheetah: Calculate time beteeen processing of data frames
		time_t	tnow;
		double	dtime, datarate;
		time(&tnow);
		
		dtime = difftime(tnow, cheetahGlobal.tlast);
		if(dtime > 1.) {
			datarate = (frameNumber - cheetahGlobal.lastTimingFrame)/dtime;
			cheetahGlobal.lastTimingFrame = frameNumber;
			time(&cheetahGlobal.tlast);
			cheetahGlobal.datarate = datarate;
		}
			
		// Cheetah: Create a new eventData structure in which to place all information
		cEventData	*eventData;
		eventData = cheetahNewEvent(&cheetahGlobal);
		ntriggers++;
            
        // Cheetah: Populate event structure with meta-data
		eventData->frameNumber = frameNumber;
		eventData->fiducial = frameNumber;
		eventData->runNumber = runNumber;
		eventData->nPeaks = 0;
		eventData->pumpLaserCode = 0;
		eventData->pumpLaserDelay = 0;
		eventData->photonEnergyeV = 12390 / wavelength;
		eventData->wavelengthA = wavelength;
		eventData->pGlobal = &cheetahGlobal;
		
        // Cheetah: Copy image data into
        
		// 1. open the required data
		int block_number = start_block_number + frameNumber / number_per_block;
		int frame_in_block = frameNumber % number_per_block;
		fprintf(stderr, "frame %d is in block %d, frame %d (all 0-indexed).\n", 
				frameNumber, block_number, frame_in_block);
	
		snprintf(data_name, 20, "data_%06d", block_number); 
		data = H5Dopen2(group, data_name, H5P_DEFAULT);
		dataspace = H5Dget_space(data);
		if (data < 0) {
			fprintf(stderr, "failed to open /entry/%s\n", data_name);
			return -1;
		}
		if (H5Sget_simple_extent_ndims(dataspace) != 3) {
			fprintf(stderr, "Dimension of /entry/%s is not 3!\n", data_name);
			return -1;    
		}
		
		// 2. get the frame
		H5Sget_simple_extent_dims(dataspace, dims, NULL);
		hsize_t offset_in[3] = {frame_in_block, 0, 0};
		hsize_t offset_out[3] = {0, 0, 0};
		hsize_t count[3] = {1, ypixels, xpixels};
		hid_t memspace = H5Screate_simple(3, dims, NULL);   
		if (memspace < 0) {
			fprintf(stderr, "failed to create memspace\n");
			return -1;
		}
		int ret;
		ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset_in, NULL, 
								  count, NULL);
		if (ret < 0) {
			fprintf(stderr, "select_hyperslab for file failed\n");
			return -1;
		}
		ret = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, 
								  count, NULL);
		if (ret < 0) {
			fprintf(stderr, "select_hyperslab for memory failed\n");
			return -1;
		}
		H5Dread(data, H5T_NATIVE_SHORT, memspace, dataspace, H5P_DEFAULT, buf);
		H5Sclose(dataspace);
		H5Sclose(memspace);
		H5Dclose(data);

		// 3. copy to the event data
		long    detID = 0;
		long    pix_nn = cheetahGlobal.detector[detID].pix_nn;

//		int underflow = 0, overflow = 0;
		for(long ii=0; ii<pix_nn; ii++) {
			long tmp = buf[ii];
//			if (ii < 10000) printf("%d ", tmp);
			if (tmp < 0) {
//				underflow++; //printf("%ld ", tmp);
				tmp = 0; 
			} else if (tmp >= SHRT_MAX) { // SHRT_MAX == 32767 seems panel gap. 
//				overflow++;
				tmp = 0; 
			}
			eventData->detector[detID].data_raw16[ii] = (uint16_t)tmp;
		}
//		printf("underflow %d overflow %d\n", underflow, overflow);
		cheetahProcessEventMultithreaded(&cheetahGlobal, eventData);

		frameNumber += global->frameStep;
    }

	// Finish

	cheetahExit(&cheetahGlobal);

	H5Gclose(group);
	H5Fclose(hdf);
	free(buf);
	
	time_t endT;
	time(&endT);
	double dif = difftime(endT,startT);
	std::cout << "time taken: " << dif << " seconds\n";	
	std::cout << "Clean exit!\n";
    return 0;
}

/*
 *  Print some useful information
 */
void print_help(void){
    std::cout << "Cheetah interface for Eiger\n";
    std::cout << "authors & developers, March 2019-\n";
    std::cout << std::endl;
    std::cout << "usage: cheetah-eiger -i <INIFILE> proteinXXX_master.h5 \n";
    std::cout << std::endl;
    std::cout << "\t--inifile=<file>     Specifies cheetah.ini file to use\n";
    //std::cout << "\t--experiment=<name>  String specifying the experiment name (used for lableling and setting the file layout)\n";
    std::cout << "\t--first=<n>          Start processing from the <first> frame (inclusive, frame numbering starts with 0)\n";
    std::cout << "\t--step=<n>           Process only each <step>th frame\n";
	std::cout << "\t--last=<n>           End processing with the <last> frame (inclusive, set to -1 to process all frames)\n";
    std::cout << std::endl;
    std::cout << "End of help\n";
}

/*
 *	Configuration parser (getopt_long)
 */
void parse_config(int argc, char *argv[], tCheetahEigerparams *global) {

	if(argc<3) {
		printf("At least two arguments required.\n");
		print_help();
		exit(1);
	}
	
	global->masterFile = NULL;
	global->iniFile = NULL;
	global->frameFirst = 0;
	global->frameStep = 1;
	global->frameLast = -1;

	// legacy cheetah-eiger: if we have only two arguments we are done
	if(argc<=3) {
		global->masterFile = argv[1];
		global->iniFile = argv[2];
		return;
	}

	// Add getopt-long options
    // three legitimate values: no_argument, required_argument and optional_argument
	const struct option longOpts[] = {
		/*{ "master",  required_argument, NULL, 'm' },*/
		{ "inifile", required_argument, NULL, 'i' },
        { "first",   required_argument, NULL, 'f' },
        { "step",    required_argument, NULL, 's' },
        { "last",    required_argument, NULL, 'l' },
		{ "verbose", no_argument, NULL, 'v' },
		{ "help",    no_argument, NULL, 'h' },
		{ NULL,      no_argument, NULL,  0  }
	};
	const char optString[] = "i:f:s:l:vh?";
	
	int opt;
	int longIndex;
	while( (opt=getopt_long(argc, argv, optString, longOpts, &longIndex )) != -1 ) {
		switch( opt ) {
		case 'v':
			//global->verbose++;
			break;
		case 'h':   /* fall-through is intentional */
		case '?':
			print_help();
			exit(1);
			break;
		case 'i':
			global->iniFile = optarg;
			std::cout << "cheetah.ini file set to " << global->iniFile << std::endl;
			break;
        case 'f':
			global->frameFirst = atol(optarg);
			std::cout << "the first frame number to process set to " << global->frameFirst << std::endl;
			break;
		case 's':
			global->frameStep = atol(optarg);
			std::cout << "the step for frames processing set to " << global->frameStep << std::endl;
			break;
        case 'l':
			global->frameLast = atol(optarg);
			std::cout << "the last frame number to process set to " << global->frameLast << std::endl;
			break;
		case 0: /* long option without a short arg */
			// TODO
			break;
		default:
			/* You won't actually get here. */
			break;
		}
	}
	
	// This is where unprocessed arguments end up
	//std::cout << "optind: " << optind << std::endl;
	if( argc-optind == 1 ) {
		global->masterFile = argv[optind];
	} else {
		std::cout << "Number of unprocessed arguments: " << argc-optind << std::endl;
		std::cout << "There should be exactly one for the master data file" << std::endl;
		exit(1);
	}
}

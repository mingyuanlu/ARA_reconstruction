# RIvR
----------
Radio Interferometric Neutrino Reconstruction. This is a software that utilizes waveform interferometry and the Radiospline library to attempt the reconstruction of neutrino events from the Askaryan Radio Array (ARA) experiment.

# Requirements
----------
1. Radiospline. See https://github.com/WIPACrepo/radiospline
2. CFITSIO
3. Healpix
4. ROOT, required components: MathMore, Gui
5. clFFT
6. AraRoot. See http://www.hep.ucl.ac.uk/uhen/ara/araroot/
7. AraSim. See http://www.physics.ohio-state.edu/~connolly/AraSim/arainstr.html

# Building
----------
To make use of the software tools:
1. In your .bashrc or equivalent file, include line

		$ scl enable devtoolset-2 $SHELL

2. Set up a working AraRoot copy. Set up a working AraSim copy in the same directory, otherwise AraSim events analysis will not run properly. An area containing these (AraRoot analysis + AraSim) can be found at the WIPAC server at /data/user/mlu27/analysis/RIvR_skeleton.

3. In this common work area, create libAraSim.a with

		$ ar rcs libAraSim.a *.o

4. Modify CMakeLists.txt according to the enclosed CMakeLists_example.txt. Point the include paths to where the softwares tools are installed. In particular, if you set up the common work area for AraRoot and AraSim as in #1, ${ARASIM_INCLUDE_DIR} should be the same as ${CURRENT_WORK_DIR}. Also, make sure the line

		$ find_package(ROOT REQUIRED COMPONENTS MathMore Gui)

	precedes the two instances of the command

		$ ROOT_GENERATE_DICTIONARY(...)

5. In CMakeLists.txt, add lines

		$ add_executable(analysis analysis.cxx)
		$ target_link_libraries(analysis ${ARAEVENT_LIBRARIES} ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES} ${ExtraLibs} ${HealpixLibs})

 This enables you to run the analysis.

6. Do

		$ sh INSTALL.sh 0

# Usage
----------
The "analysis" executable is generated in #5 above.
To reconstruct AraSim events, do:

		./build/analysis ${OUTPUT_DIR} recoSetupFile_default.txt 0 ${PATH_TO_SIMULATION}/AraOut.root

To reconstruct real data, do:

		./build/analysis ${OUTPUT_DIR} recoSetupFile_default.txt 0 ${PATH_TO_REAL_DATA}/eventXXXX.root ${PATH_TO_PEDESTAL_FILE}/pedestalValues.run00XXXX.dat

In the above examples, '0' is the run number given to the analysis runs. The analysis outputs will be a ROOT file named recoOut.${recoSetupFile}.run${RunNumber}.root, plus interferometric skymaps (in FITS format) named recoSkymap.${recoSetupFile}.run${RunNumber}(.constantN/2ndRay).fits. These output files will be located in ${OUTPUT_DIR}.

# Configuring the reconstruction
----------
All configuration parameters for the reconstruction are defined and explained in recoSettings.h. recoSetupFile_default.txt sets all parameters to the default values. To reconstruct with different parameters, modify recoSetupFile and pass it as an argument to analysis.

# OpenCL kernels
----------
All OpenCL kernels used are implemented in kernel_3D_analysis.c

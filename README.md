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

1. Set up a working AraRoot copy. Set up a working AraSim copy.
2. In the AraSim work area, create libAraSim.a with

		$ ar rcs libAraSim.a *.o

3. Go to the AraRoot work area. Modify CMakeLists.txt according to the enclosed CMakeLists_example.txt. Point the include paths to where the softwares tools are installed. In particular, point the ARASIM_INCLUDE_DIR to the AraSim work area in #2.
4. In CMakeLists.txt, add lines

		$ add_executable(analysis analysis.cxx)
		$ target_link_libraries(analysis ${ARAEVENT_LIBRARIES} ${ROOT_LIBRARIES} ${ZLIB_LIBRARIES} ${ExtraLibs} ${HealpixLibs})
 
 This enables you to run the analysis.

5. Do
 
		$ sh INSTALL.sh 0
		
# Usage
----------
The "analysis" executable is generated in #5 above. 
To reconstruct AraSim events, do:

		./build/analysis recoSetupFile_default.txt 0 ${PATH_TO_SIMULATION}/AraOut.root

To reconstruct real data, do:

		./build/analysis recoSetupFile_default.txt 0 ${PATH_TO_REAL_DATA}/eventXXXX.root ${PATH_TO_PEDESTAL_FILE}/pedestalValues.run00XXXX.dat
		
# Configuring the reconstruction
----------
All configuration parameters for the reconstruction are defined and explained in recoSettings.h. recoSetupFile_default.txt sets all parameters to the default values. To reconstruct with different parameters, modify recoSetupFile and pass it as an argument to analysis.

# OpenCL kernels
----------
All OpenCL kernels used are implemented in kernel_3D_analysis.c



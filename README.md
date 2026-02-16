# MRISimVecField: MRI Simulation and Reconstruction Framework for Magnetic Vector Fields

This simulator is a MATLAB simulation package capable of full consideration of angular magnetic field deviations. An interface to magnetic field simulation software like CST is fairly easy via import of txt-files. A more detailed description of this repository is published in:
	- F.Bschorr et al. MRI simulation and reconstruction framework for arbitrarily-oriented encoding and transmit/receive magnetic vector fields. Zeitschrift für Medizinische Physik. https://doi.org/10.1016/j.zemedi.2026.02.002


#### General:
- Written in MATLAB
- Freely available under BSD License
- Simulator for Arbitrary Magnetic Vector Fields

- For visualization purposes arrayShow is used (Version 0.35): https://github.com/tsumpf/arrShow
- The following MATLAB Toolboxes ares used:
	- Signal Processing Toolbox
	- Image Processing Toolbox
	- Communications Toolbox
	- Parallel Computing Toolbox V7.6

#### Installation
- Download / clone this repository
- Open SimMain.m in MATLAB
- Run the script
- The result (reco_rho_straight_rrr_img) should be similar to that saved in "ExampleImg_reconstructed.mat"


More documentation can be found in "doc_V1.pdf"

#### Citation
If you intend to use any part of the code in the repository, please cite

#### Acknowledgments
NOLIMS was developed within the ExCaVI group of Ulm University and Ulm University Medical Center.

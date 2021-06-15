# developed tools for PTMs identification and characterization improvement. 


Post-translational modifications play a very important role in multiple diseases. In order to improve its identification and quantification, 
several tools have been developed in this work: FDRFilterer, DM0Solver, TrunkSolver, PDmTablemaker, SiteSolver and Sticker. These modules are
independent but we recommend their use in the previously initiated order. Before executing these programs, the row files of the mass spectrometer
have been passed through Comet-PTM (https://github.com/CNIC-Proteomics/Comet-PTM) and some programs from a work in development (DMCalibrator, 
DMModeller, PeakInspector, PeakSelector, PeakAssignator and PeakFDRer https://github.com/CNIC-Proteomics/SHIFTS-4). It should be noted that before
and after SiteSolver, PeakAssigantor must be run again.

From the output of PeakFDRer, the FDRfilterer will filter the FDR at the global, local, and peak levels. DM0Solver and TrunkSolver recover unmodified 
peptides taking into account what is indicated by the user in the configuration file. PDMTableMaker will help SiteSolver to find the correct residue where
the modification is located, with the support of some lists created by the user. In addition, PDmTableMaker together with Sticker that will label the 
modifications according to lists will lead to a faster and easier data inspection. 

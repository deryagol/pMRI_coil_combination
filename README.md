# pMRI_coil_combination
Matlab codes and Reconstruction Examples of the paper:
Derya Gol Gungor, Lee C Potter, "A Subspace-based Coil Combination Method for Phased-array Magnetic Resonance Imaging", submitted to Magnetic Resonance in Medicine in 2015.

demo.m is a demo file for the proposed coil combination approach 
ptl_convmtx2.m calculates 2D partial convolution matrix - accuracy is not guaranteed!
crop.m, sos.m, sos_kspace.m can be obtained form Michael Lustig's SPIRiT package. 
Spine dataset used in the decome can be obtained from PULSAR Toolbox. Webpage addresses are provided in demo.m. 
experiments folder includes reconstruction results applied to 2013 ISMRM Recon Challenge cardiac datasets. *_unmod.gif refers to BCC reconstructed images. 

Please acknowledge this code and cite the paper appropriately if used/distributed/modified. 
Contact with Derya at deryagol@gmail.com for any questions and comments. 

(c) Derya Gol Gungor, 2015

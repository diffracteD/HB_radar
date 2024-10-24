Base Program Name HB_Radar  
Basic Principle: Finds hydrogen bonds (H-bonds) from X-ray Crystal structures and calculates Potential Energy Surface using Density Functional Theory.  
>Working protocol:  
>Harvest .pdb files from RCSB-PDB   
>Adds hydrogen to eah pdb and outputs .atom file  
>Searches each file and harvest coordinates involved in H-bond formation (updated version can isolate nearby coordinates of H-bond, given certain radius).  
>Uses NWChem program as a wrapper to calculate QM single point energy (SPE).   
>Writes Fortran scripts independently to run NWChem as a sub-process and returns the value to specific array once a SPE calculation finishes.  
>**Update.v.3: Program can iteratively self-write and try different Fortran scripts for NWChem to try differnt grid size and exchange-correlation function, till a calculation converges.  
>Outputs the H-bond coordinates and corresponding SPE in single file.  
>... repeates for every file supplied  
>Once SPE calculations are over, it plots all the information w.r.t. H-bond's geometric parameters.  

Smaple output and studies performed: https://onlinelibrary.wiley.com/doi/10.1002/prot.25271  


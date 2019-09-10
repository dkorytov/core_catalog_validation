# core_catalog_validation
Validation for the Core Catalog

This validation suite performs a number of sanity checks on core catalog outputs. These include confirmation of uniqueness of identifiers, distributions of positions and infall quantities, consistency checks between adjacent timesteps and against merger tree files, and appropriate assignation of the central flag.   

# Usage: 
python check_cores.py time_step1 time_step2
 
Note that you may need to update the paths to genericio.py and to the core catalog outputs in the check_cores.py file. The default points to the location of the most up-to-date core catalogs on datastar as of 09/10/19. 


# Example output: python check_cores.py 88 86

Minimum infall mass has positive value: 21.0 particles

List of infall steps:
[43 44 45 46 48 49 50 52 53 54 56 57 59 60 62 63 65 67 68 70 72 74 76 77 79
 81 84 86 88]
 
Minimum infall step has positive value: 43

('Position', 'x')
	 min/max: -0.0487365722656/256.076782227
	 non-finite vals: 0, 0.0%
	 
('Position', 'y')
	 min/max: -0.0309143066406/256.050720215
	 non-finite vals: 0, 0.0%
	 
('Position', 'z')
	 min/max: -0.0411987304688/256.08770752
	 non-finite vals: 0, 0.0%
	 
Maximum number of cores in a halo is 44

Number of halos is 475950

All halos have a central core

Core tag is unique

Infall tree node index is unique

Warning, merger tree files should only match properly at step 499 (at higher redshifts N_mergertree > N_centralcores due to fragmentation). As you are looking at a higher redshift, expect a ~5% level of missing cores with a range of masses (including 0 mass fragments).

Number of halos in merger tree file = 496789

Number of central cores in core file = 475950

Percentage of cores present in earlier timestep missing in subsequent timestep = 1.93906515383

Total number of issues found = 0



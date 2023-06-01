# Charm
*Ch*romosome re*ar*rangement *m*odeler.
There ara python the based tool to simulate of Hi-C-map with the predetermened chromosomal rearrangments. This tools allow simulate such structural variations as CNVs, inversions, translocation and the extra-chromosome creation.

## Requirements
1) Python >= 3.7 with numpy module
2) Jucier Tools https://github.com/aidenlab/juicer (for dumping the contacts from the .hic-files and creation the new .hic-files)

## Quick Start
### Unix
charm.sh [-i ini_file] [-S stage] 
* [ini_file]: the path to ini-file contains all needed paramaters to simulate SVs: path to work directory, hic-file, uniq SV id, model paramaters and others, the ini-file description see in test.ini
* [stage]: optitional, must be one from "pre+","SVs+","sim+","lift+","wt+","hic" (default "pre+")
  - "pre+" is default parameter, from hic-file with wild_type to hic-file with SVs.
  - Use "SVs+" when database with contact statistics are done, but the file with rearrangment description has not yet been created 
  - Use "sim+" when the rearrangment description and database are done, but the contacts of mutant genome are not yet been simulated
  - Use "lift+" when the contacts of mutant genome are simulated ("in_mut.[simulation_id]" files created), but this contacts are not yet lifted on the reference genome
  - Use "wt+" when the simulation is fully processed ("in_mut.[simulation_id]" created), but wild-type replicas are not done
  - Use "hic" when the all previously stages are successfully complited, but .hic-file is not created
* advanced settings for [stage]:
  - Use "pre" to creation the database with contact statistics only
  - Use "SVs" to creation the rearrangment description files only; in this mode Charm can process the file containing any number of independent SVs 
  - Use "sim" to only simulate contacts within *mutanted* genome based on handly defined database
  - Use "lift" to lifover contacts from handly defined mutant genome
  - Use "wt" to simulated wild-type replicas
### others OS
python3 scripts/charm_manager.py [-S stage] -i [ini_file]

## Test dataset
1) dowload testdataset folder
2) open **TEST_EXAMPLE.ini**
3) write to ini-file the path to you work directory in the [global] section, on the "work_dir" key
4) write to ini-file the path to juicer tools .jar-file in the [global] section, on the "path_to_juicertools" key  
4.5 optitional, if path to java not defined, write to ini-file the path to java directory in the [global] section, on the "path_to_java_dir" key 
5) write to ini-file the path to test.chr.sizes file in the [global] section, on the "chrom_sizes" key 
6) write to ini-file the path to test.hic file in the [prerpoceesing] section, on the "path_to_hic" key
7) write to ini-file the path to test.svs_list.txt in the [SVs] section, on the "path_to_svs_list" key
8) don't change other paramaters
9) run Charm with the "-S pre+" parameter

## First run
1) open **TEST_EXAMPLE.ini**
2) write to ini-file the path to you work directory in the [global] section, on the "work_dir" key
3) write to ini-file the unique simulation ID in the [global] section, on the "simulation_id" key
4) write to ini-file the path to juicer tools .jar-file in the [global] section, on the "path_to_juicertools" key  
4.5 optitional, if path to java not defined, write to ini-file the path to java directory in the [global] section, on the "path_to_java_dir" key 
5) write to ini-file the path to **you_chromosome_sizes_file** file in the [global] section, on the "chrom_sizes" key; the file format see below 
6) write to ini-file the path to **you_hic_map.hic** file in the [prerpoceesing] section, on the "path_to_hic" key
7) write to ini-file the path to **you_SVs_list_file** in the [SVs] section, on the "path_to_svs_list" key; the file format see below
8) change other parameters as you wish
9) run Charm

###
The whole genome Hi-C datasets of chromosomal rearrangments can be founded [here](https://genedev.bionet.nsc.ru/ftp/by_Project/Charm/wgHi-C/) 


### The chromosome sizes file
This file contains chromosome sizes (see https://github.com/aidenlab/juicer/wiki/Pre). The chromosome names and chromosome sizes must be correspondent to the chromosome sizes and chromosome names in .hic-file. 
File format (see the example "test.chr.sizes")
```
<chromosome name> <chromosome size bp>
```

### The SVs description 
To simulated SVs, Charm requers the description of rearrangment given in the correct format
File format (see ethe xample "test.svs_list.txt" in the testdataset folder)
```
<reference genome id> <mutant genome uniq id> <chromosome> <coordinate of rearrangmnet locus start> <coordinate of rearrangmnet locus end> <the indicator> <new chromosome> <coordinate of new position of locus> <copy number of locus on OLD position> <copy number of locus on NEW position>
```
The indicator variants:
  - Use "->" for the plain SVs, this indicator designes the start and the end of SVs description
  - Use "!>" for the start of description of complicated SVs
  - Use ">>" for the continuation of description SVs
  - Use ">!" for the end of description of complicated SVs, the all lines between "!>" and ">!" are procesed by Charm as one SV.
 
Examples:
  1) the general translocation; the moving of locus chr1:1,000,000-2,000,000 to chr1:7,000,000 position: 
  ```
  hg19 my_translocation chr1 1000000 2000000 -> chr1 7000000 0 1
  ```
  2) the tandem duplication of locus chr1:1,000,000-2,000,000
  ```
  hg19 my_duplication1 chr1 1000000 2000000 -> chr1 2000000 1 1
  ```
  or
  ```
  hg19 my_duplication2 chr1 1000000 2000000 -> chr1 2000000 0 2
  ```
  3) the deletion of locus chr1:1,000,000-2,000,000 
  ```
  hg19 my_deletion chr1 1000000 2000000 -> chr1 1000000 0 0
  ```
  4) the cnv, the repeating locus chr1:1,000,000-2,000,000  10 times
  ```
  hg19 my_cnv10 chr1 1000000 2000000 -> chr1 1000000 0 10
  ```
  5) the inversion of locus chr1:1,000,000-2,000,000 
  ```
  hg19 my_inversion chr1 1000000 2000000 -> chr1 2000000 0 *-*1
  ```
  6) the translocation of locus chr1:1,000,000-2,000,000 to chromosome end with the locus saving on old position and cnv x10 on new position
  ```
  hg19 my_translocation2 chr1 1000000 2000000 -> chr1 + 1 10
  ```
  7) the complicated rearrangment: the translocation and inversion within one genome
  ```
  hg19 my_complicated2 chr1 1000000 2000000 !> chr1 7000000 0 1
  hg19 my_complicated2 chr1 6500000 7500000 >! chr1 6500000 0 -1
  ```
  8) the very complicated rearrangment: the translocation with cnv and creation of new chromosome
  ```
  hg19 my_complicated3 chr1 1000000 2000000 !> chrNew + 1 3
  hg19 my_complicated3 chr1 3000000 4000000 >> chrNew + 1 -3
  hg19 my_complicated3 chr1 5000000 7500000 >! chrNew + 1 -3
  ```

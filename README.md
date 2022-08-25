# Multispecies_occupancy_lynx_wildcat_wolf
 Code and data for Dyck et al., 2022 - Dracula’s ménagerie: A multispecies occupancy analysis of lynx, wildcat, and wolf in the Romanian Carpathians

<hr>

### GENERAL INFORMATION

**Title:** Data and code for - Dracula’s ménagerie: A multispecies occupancy analysis of lynx, wildcat, and wolf in the Romanian Carpathians

**Author Information:**  
    Principal Investigator Contact Information  
		Name: Marissa A. Dyck  
		Institution: Ohio University  
		Address: Athens, OH USA  
		Email: marissadyck17@gmail.com  

**Date of data collection:** 2018 - 2020 (seasonally)

**Geographic location of data collection:** Southern Carpathians, Romania

**Information about funding sources that supported the collection of the data:** Field work was funded by the OAK Foundation grant number OCAY-11-136 and by the European Commission through the Operational Programme 'Environment', grant number SMIS 102086. VDP was partially supported by a grant from the Romanian National Authority for Scientific Research (PN-III-P1-1.1-TE-2019-0835). Travel for MAD to Romania was provided by the Ohio University College of Arts and Sciences.


### SHARING/ACCESS INFORMATION

**Licenses/restrictions placed on the data:** None

**Link to publications that cite:**  <a href = "https://onlinelibrary.wiley.com/doi/10.1002/ece3.8921" target = "_blank">[PDF WILEY]</a> 

**Recommended citation for this dataset:** Dyck, Marissa (2022), Data from: Dracula’s ménagerie: A multispecies occupancy analysis of lynx, wildcat, and wolf in the Romanian Carpathians, Dryad, Dataset, https://doi.org/10.5061/dryad.47d7wm3gp 

**Recommended citation for this manuscript:** Dyck, Marissa A., et al. "Dracula’s ménagerie: A multispecies occupancy analysis of lynx, wildcat, and wolf in the Romanian Carpathians." *Ecology and Evolution* 12.5 (2022): e8921.



### DATA & FILE OVERVIEW

**File List:**  

* <span style = "color: #7B0F17;">**multi_occ_LWW.Rproj**</span> R project to run code for Dyck et al., 2022
* <span style = "color: #7B0F17;">**multi_occ_LWW_script.R**</span> R script with code for multispecies occupancy analysis for Dyck et al., 2022
* <span style = "color: #7B0F17;">**om_predict_fix.R**</span> source script provided by Ken Kellner to fix package bug related to using the predict function
* <span style = "color: #7B0F17;">**species_matrix**</span> contains presence/absence data for all 3 species
* <span style = "color: #7B0F17;">**cams_data**</span> contains landcover and environmental data associated with the location of each camera trap
* <span style = "color: #7B0F17;">**Trap_effort**</span> contains data for camera trap activity in each session. If one camera was working every day of a session the value will = 14; any days where both cameras were inactive subtracts 1 from this value.

*'winter' files contain data from 2018/12/17 to 2019/3/31 and 'autumn' contain data from 2019/10/9 to 2020/1/15* 


### METHODOLOGICAL INFORMATION

**Description of methods used for collection/generation of data:** Methodological details are provided in the paper <a href = "https://onlinelibrary.wiley.com/doi/10.1002/ece3.8921" target = "_blank">[PDF WILEY]</a> 

**Methods for processing the data:** We used multi species occupancy models to process the data <a href = "https://onlinelibrary.wiley.com/doi/10.1002/ece3.8921" target = "_blank">[See Rota et al., 2016]</a>. Full details can be found in the paper.


**Instrument- or software-specific information needed to interpret the data:** We used program R version 3.5.1 with package unmarked. 
* Download R <a href = "https://cran.r-project.org/bin/windows/" target = "_blank">[Windows link]</a> <a href = "https://cran.r-project.org/bin/macosx/" target = "_blank">[Mac link]</a>
* Downlad R Studio <a href = "https://www.rstudio.com/products/rstudio/" target = "_blank">[link]</a>

**People involved with sample collection, processing, analysis and/or submission:** The data was collected by an experienced team of wildlife rangers from  <a href = "https://www.carpathia.org/" target = "_blank">Foundation Conservation Carpathia</a> (FCC) with the help of volunteers. 

### DATA-SPECIFIC INFORMATION FOR: [cams_data_autumn]

* **Number of variables:** 26
*  **Number of cases/rows:** 76

**Variable List:**

* <span style = "color: #E8A12C;">**TrapCode**</span>, unique identifier for each camera trap station
* <span style = "color: #E8A12C;">**X/Y**</span>, GPS coordinates for each camera trap station
* <span style = "color: #E8A12C;">**Z**</span>, altitude for each camera trap station 
* <span style = "color: #E8A12C;">**Impact**</span>, description of anthropogenic impact at camera trap station
* <span style = "color: #E8A12C;">**distnatlro**</span>, distance to nearest national road (equivalent of highways), meters
* <span style = "color: #E8A12C;">**distsettle**</span>, distance to nearest human settlement, meters
* <span style = "color: #E8A12C;">**diststream**</span>, distance to nearest stream, meters
* <span style = "color: #E8A12C;">**denslocalr**</span>, the density of local roads measured using a moving window approach, km/km2
* <span style = "color: #E8A12C;">**distlocalr**</span>, distance to nearest local road, meters
* <span style = "color: #E8A12C;">**TRI5**</span>, Terrain Roughness Index within a 50m buffer, <a href = "http://download.osgeo.org/qgis/doc/reference-docs/Terrain_Ruggedness_Index.pdf" target = "_blank">[see Riley et al., 1999]</a>.
* <span style = "color: #E8A12C;">**CLC2018**</span>, the landcover type in 2018 as defined by <a href = "https://land.copernicus.eu/eagle/files/eagle-related-projects/pt_clc-conversion-to-fao-lccs3_dec2010" target = "_blank">CORINE Land Cover Classification system</a>.
* <span style = "color: #E8A12C;">**CLC112-512**</span>, proportions of each CORINE land cover type within each 2.7 x 2.7km cell, see <a href = "https://land.copernicus.eu/eagle/files/eagle-related-projects/pt_clc-conversion-to-fao-lccs3_dec2010" target = "_blank">[HERE]</a> for full list of all land cover types.


### DATA-SPECIFIC INFORMATION FOR: [cams_data_winter]

* **Number of variables:** 26
* **Number of cases/rows:** 64

**Variable List:**

* <span style = "color: #E8A12C;">**TrapCode**</span>, unique identifier for each camera trap station
* <span style = "color: #E8A12C;">**X/Y**</span>, GPS coordinates for each camera trap station
* <span style = "color: #E8A12C;">**Z**</span>, altitude for each camera trap station 
* <span style = "color: #E8A12C;">**Impact**</span>, description of anthropogenic impact at camera trap station
* <span style = "color: #E8A12C;">**distnatlro**</span>, distance to nearest national road (equivalent of highways), meters
* <span style = "color: #E8A12C;">**distsettle**</span>, distance to nearest human settlement, meters
* <span style = "color: #E8A12C;">**diststream**</span>, distance to nearest stream, meters
* <span style = "color: #E8A12C;">**denslocalr**</span>, the density of local roads measured using a moving window approach, km/km2
* <span style = "color: #E8A12C;">**distlocalr**</span>, distance to nearest local road, meters
* <span style = "color: #E8A12C;">**TRI5**</span>, Terrain Roughness Index within a 50m buffer, <a href = "http://download.osgeo.org/qgis/doc/reference-docs/Terrain_Ruggedness_Index.pdf" target = "_blank">[see Riley et al., 1999]</a>.
* <span style = "color: #E8A12C;">**CLC2018**</span>, the landcover type in 2018 as defined by <a href = "https://land.copernicus.eu/eagle/files/eagle-related-projects/pt_clc-conversion-to-fao-lccs3_dec2010" target = "_blank">CORINE Land Cover Classification system</a>.
* <span style = "color: #E8A12C;">**CLC112-512**</span>, proportions of each CORINE land cover type within each 2.7 x 2.7km cell, see <a href = "https://land.copernicus.eu/eagle/files/eagle-related-projects/pt_clc-conversion-to-fao-lccs3_dec2010" target = "_blank">[HERE]</a> for full list of all land cover types.

### DATA-SPECIFIC INFORMATION FOR: [species_matrix_autumn]

* **Number of variables:** 22
* **Number of cases/rows:** 76

**Variable List:**

* <span style = "color: #E8A12C;">**TrapCode**</span>, unique identifier for each camera trap station
* <span style = "color: #E8A12C;">**Lynx_1-Lynx_7**</span>, presence/absences data for Eurasian lynx (*Lynx lynx*) for each session
* <span style = "color: #E8A12C;">**Wildcat_1-Wildcat_7**</span>, presence/absence data for European wildcat (*Felis silvestris*) for each session
* <span style = "color: #E8A12C;">**Wolf_1-Wolf_7**</span>, presence/absence data for grey wolf (*Canis lupus*) for each session

### DATA-SPECIFIC INFORMATION FOR: [species_matrix_winter]

* **Number of variables:** 25
* **Number of cases/rows:** 64

**Variable List:**

* <span style = "color: #E8A12C;">**TrapCode**</span>, unique identifier for each camera trap station
* <span style = "color: #E8A12C;">**Lynx_1-Lynx_8**</span>, presence/absences data for Eurasian lynx (*Lynx lynx*) for each session
* <span style = "color: #E8A12C;">**Wildcat_1-Wildcat_8**</span>, presence/absence data for European wildcat (*Felis silvestris*) for each session
* <span style = "color: #E8A12C;">**Wolf_1-Wolf_8**</span>, presence/absence data for grey wolf (*Canis lupus*) for each session

### DATA-SPECIFIC INFORMATION FOR: [Trap_effort_autumn]

* **Number of variables:** 8
* **Number of cases/rows:** 76

**Variable List:** 

**Variable List:** 

* <span style = "color: #E8A12C;">**TrapCode**</span>, unique identifier for each camera trap station
* <span style = "color: #E8A12C;">**Session_1-Session_7**</span>, camera trap activity for each session. *If one camera was working every day of a session the value will = 14; any days where both cameras were inactive subtracts 1 from this value*.

### DATA-SPECIFIC INFORMATION FOR: [Trap_effort_winter]

* **Number of variables:** 9
* **Number of cases/rows:** 64

**Variable List:** 

* <span style = "color: #E8A12C;">**TrapCode**</span>, unique identifier for each camera trap station
* <span style = "color: #E8A12C;">**Session_1-Session_8**</span>, camera trap activity for each session. *If one camera was working every day of a session the value will = 14; any days where both cameras were inactive subtracts 1 from this value*.

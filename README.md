The following files allow one to reproduce analyses in the manuscript entitled "Warming-induced divergence in native rhizosphere fungi amplifies invasion-driven soil biotic homogenization".

DATA & FILE OVERVIEW

***In Data folder***

The experimental data are stored in Figshare [![DOI](https://zenodo.org/badge/DOI/10.6084/m9.figshare.28139549.svg)](https://doi.org/10.6084/m9.figshare.28139549.v2).
Before the manuscript is officially published, experimental and analytical data must remain confidential. 
If needed, please contact the first or corresponding author in advance to obtain the relevant experimental data. 
All data will be made available upon acceptance of the manuscript.

*List of experimental data files (.xlsx)*

    * 1. Field_ASVs_raw_data.xlsx  
    * 2. Field_data_group.xlsx  
    * 3. ASV_tax_information.xlsx (field survey)
    * 4. cohesion_data.xlsx
    * 5. Greenhouse_ASVs_row_data.xlsx  
    * 6. Greenhouse_data_group.xlsx  
    * 7. FungalTraits.xlsx  
    
*List of phylogenetic tree data files (.newick)*  

    * 1. IQ_tree_plant_2025.newick  (Phylogenetic tree of studying plant species)

***In Code folder***

The names of R-scripts correspond to the statistical analysis and visualization of the corresponding figures or tables in this manuscript.

*List of R-scripts*

    * 1. Figure 1.tiff (Draw in ppt)
    * 2. Figure 2 & S2a & S2c.R  
    * 3. Figure 3a & 3c.R  
    * 4. Figure 3b & 3d.R  
    * 5. Figure 4.R  
    * 6. Figure S1.R
    * 7. Figure S2b & S2d.R  
    * 8. Figure S3
    * 9. Figure S4.R  
    * 10. Table S2 & Table S3 (field survey part).R  
    * 11. Table S2 & Table S3 (greenhouse exp. part).R  
    * 12. Table S4.R
    
**Data-specific onformation for:** ***Field_ASVs_row_data.xlsx***

    * Abundance table of raw sequencing data of rhizosphere fungi from field survey (not rarefied to minimum sample size)

**Data-specific onformation for:** ***Field_data_group.xlsx***

    Variable list	 Description
    * Sample_ID: Sample id of  plant rhizosphere soil 
    * Latitude: Latitude of sampling point
    * Longitude: Longitude of sampling point
    * Year: Sampling year
    * Site: Name of sampling site
    * Chinese_name: Chinese name of study species
    * Species: Latin name of study species
    * Genus: Genus name of research species
    * Family: Family name of research species
    * Origin: Geographical origin of plants (native vs. exotic)
    * Soil_ph: Soil pH
    * Wcont: Soil water content
    * Soil_N: Soil total nitrogen content
    * Tave: Annual average temperature (℃)
    * Prec: Annual precipitation (mm)
    * Chol: Leaf chlorophyll
    * SLA: Specific leaf area (cm2 g-1)
    * LDMC: Leaf dry matter content (g g-1)
    * SRL: Specific root length (cm2 g-1)
    * FRR: Fine-to-total root mass (g g-1)
    * RMF: Root mass fraction (g g-1)
    * Rela_generalist: Relative abundance of generalist taxa
    * Rela_specialist: Relative abundance of specialist taxa
    * Site_pool: Total number of ASV per site per year
    * Field_SR: Fungal richness per soil samples
    
**Data-specific onformation for:** ***cohesion_data.xlsx***
    
    Variable list	 Description
    * Sample_ID: Sample id of  plant rhizosphere soil 
    * cohesion.pos: positive cohesion of rhizosphere fungal community networks
    * cohesion.neg: negative cohesion of rhizosphere fungal community networks
    
**Data-specific onformation for:** ***FungalTraits.xlsx***

    * The FungalTraits database (Põlme, S., Abarenkov, K., Henrik Nilsson, R., Lindahl, B.D., Clemmensen, K.E., Kauserud, H., et al. 2021. "FungalTraits: a user-friendly traits database of fungi and fungus-like stramenopiles." Fungal Diversity 105: 1-16.)

**Data-specific onformation for:** ***ASV_tax_information.xlsx***

    * Taxonomic information of fungal ASVs in rhizosphere samples during field investigation.
  
**Data-specific onformation for:** ***Greenhouse_ASVs_row_data.xlsx***

    * Abundance table of raw sequencing data of rhizosphere fungi from greenhouse experiment (not rarefied to minimum sample size)

**Data-specific onformation for:** ***Greenhouse_data_group.xlsx***

    Variable list	 Description
    * Sample_ID: Sample id of  plant rhizosphere soil 
    * Chinese_name: Chinese name of study species
    * Species: Latin name of study species
    * Genus: Genus name of research species
    * Family: Family name of research species
    * Repeats: Repeat number of soil samples
    * Origin: Geographical origin of plants (native vs. exotic)
    * Hmax: Plant height (cm)
    * Chol: Leaf chlorophyll
    * LA: Individual leaf area (cm2)
    * SLA: Specific leaf area (cm2 g-1)
    * LDMC: Leaf dry matter content (g g-1)
    * SRL: Specific root length (cm2 g-1)
    * FRR: Fine-to-total root mass (g g-1)
    * RMF: Root mass fraction (g g-1)

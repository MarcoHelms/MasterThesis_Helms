# MasterThesis_Helms
A collection of RScripts and properties of data sets used for calculations as part of my Master's Thesis




2. Creating a data frame for Linear mixed effects modeling with the following properties
# 15 columns as followed:
# Name: site                          Fill: Garden ID
# Name: gardentype_allotment          Fill: 0 or 1 (1 for allotment)
# Name: urbanization_500              Fill: Proportionate area of impervious ground in a 500 m radius around garden (ranged between 0-1)
# Name: Management_Intensity_ranged   Fill: Value of Management intensity index
# Name: plant_diversity               Fill: Values of Plant diversity (Shannon Index)
# Name: Habitat_diversity             Fill: Values of Habitat type diversity (Shannon Index)
# Name: openground                    Fill: Proportionate area of open ground within a garden (ranged between 0-1)
# Name: herblayer                     Fill: Proportionate area of herb layer within a garden (ranged between 0-1)
# Name: shrublayer                    Fill: Proportionate area of shrub layer within a garden (ranged between 0-1)
# Name: treelayer                     Fill: Proportionate area of tree layer within a garden (ranged between 0-1)
# Name: highly_managed                Fill: Proportionate area of highly managed area within a garden (ranged between 0-1)
# Name: artificial_areas              Fill: Proportionate area of artificial area within a garden (ranged between 0-1)
# Name: detritus                      Fill: Proportionate area of detritus within a garden (ranged between 0-1)
# Name: water                         Fill: Proportionate area of water within a garden (ranged between 0-1)
# Name: Taxa                          Fill: Name of one taxonomic group (e.g. Apidae)
# Richness                            Fill: Values of taxonomic richness (Number of species)
# Shannon                             Fill: Values of taxonomic diversity (Shannon Index)
# functional_richness                 Fill: Values of functional dispersion (as calculated in the previous script "fdis")
# functional_evenness                 Fill: Values of functional evenness (as calculated in the previous script "FEve.Ricotta")
# LCBD_tax_turn                       Fill: Values of taxonomic LCBD turnover (as calculated in the previous script)
# LCBD_tax_nest                       Fill: Values of taxonomic LCBD nestedness (as calculated in the previous script)
# LCBD_tax_all                        Fill: Values of taxonomic LCBD total (as calculated in the previous script)
# LCBD_fun_turn                       Fill: Values of functional LCBD turnover (as calculated in the previous script)
# LCBD_fun_nest                       Fill: Values of functional LCBD nestedness (as calculated in the previous script)
# LCBD_fun_all                        Fill: Values of functional LCBD total (as calculated in the previous script)
# 85 Rows filled with the values of the examined gardens
#
# copy the data frame eight times (arranged vertically) corresponding to the eight taxonomic groups 
# change the name in the column "Taxa" for each copy of the original data frame

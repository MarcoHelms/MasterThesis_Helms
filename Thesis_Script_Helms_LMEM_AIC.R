# Creating a data frame for Linear mixed effects modeling with the following properties
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

require(lmerTest)
require(cAIC4)
require(MASS)

# 1. Upload the dataframe ------------------------------------------------------
modeling_data <- read.table("better_garden_modeling_data_tax_fun.txt", sep="\t", row.names = 1)

descriptors=c("gardentype_allotment", "Management_Intensity_ranged", "urbanization_500", 
              "plant_diversity", "Habitat_diversity", "openground", "herblayer", "shrublayer", "treelayer",
              "highly_managed", "artificial_areas", "detritus", "water", "Taxa")

response <- "LCBD_fun_nest"

modeling_data <- modeling_data[ ,c(descriptors, response)]
colnames(modeling_data)
# 2. Check if a random effect is possible ---------------------------------
model1 <- lmer(formula = LCBD_fun_nest ~ 1  + (1 | Taxa), 
               data    = modeling_data)

cAIC(model1)
ranova(model1)

# if the results are significant proceed to the next step
# if the results are not significant proceed to step "3b" forther below

# 3. Create a first model ------------------------------------------------------
# Select the response variable for the model (e.g. LCBD_fun_nest)
# Include all predictors as fixed effects in the first model + random effects in the brackets
# run the model
model2 <- lmer(formula = LCBD_fun_nest ~ 1 + 
                 gardentype_allotment + Management_Intensity_ranged + urbanization_500 + 
                 plant_diversity + Habitat_diversity + openground +
                 herblayer + shrublayer + treelayer + highly_managed + 
                 artificial_areas + detritus + water + (1 |Taxa), 
               data    = modeling_data)

# Check results of the model and significance for each descriptor
drop1(model2)

# Check goodness of fit by using the Conditional Akaike Information Criterion
cAIC(model2)


# 4. Refine the model ----------------------------------------------------------
# Exclude the least significant predictor from the model and run the model again
model3 <- lmer(formula = LCBD_fun_nest ~ 1 + 
                 gardentype_allotment + Management_Intensity_ranged + urbanization_500 + 
                 plant_diversity + openground +
                 herblayer + shrublayer + treelayer + highly_managed + 
                 artificial_areas + detritus + water + (1 |Taxa), 
               data    = modeling_data)

# Check results of the model and significance for each descriptor
drop1(model3)

# Check goodness of fit by using the Conditional Akaike Information Criterion
# repeat this step as long as the value of the Conditional Akaike Information Criterion decreases
cAIC(model3)  

# The model with the lowest value of the Conditional Akaike Information Criterion
# is the best fitting model


# 5. Complete the final model---------------------------------------------------
# Selection of random slopes
model1_min_step <- stepcAIC(model3,
                            slopeCandidates = descriptors[-length(descriptors)],
                            numberOfPermissibleSlopes	= 1,
                            steps=1000,
                            direction="forward", data=modeling_data, trace=TRUE)

# Copy the model calculated by the previous step from the Console
model4 <- lmer(formula = LCBD_fun_nest ~ 1 + 
                 gardentype_allotment + Management_Intensity_ranged + urbanization_500 
               + plant_diversity + openground + herblayer 
               + shrublayer + treelayer + highly_managed + artificial_areas 
               + detritus + water + (1 + gardentype_allotment | Taxa) , data=modeling_data) 

#Check final properties of the model
drop1(model4)
ranova(model4)
cAIC(model4) 

#Save the model
saveRDS(model4, "Model_Garden_LCBD_fun_nest.rds")



################################################################################
# 3b. LM modeling if no random effects are possible --------------------------

# Create a first model
# Select the response variable for the model (e.g. LCBD_fun_nest)
# Include all predictors as fixed effects in the first model
# run the model

model2 <- lm(formula = LCBD_fun_nest ~ 1 + 
               gardentype_allotment + Management_Intensity_ranged + urbanization_500 + 
               plant_diversity + Habitat_diversity + openground +
               herblayer + shrublayer + treelayer + highly_managed + 
               artificial_areas + detritus + water, 
             data = na.omit(modeling_data))
summary(model2)

# Calculate a stepwise regression model
step.model <- stepAIC(model2, direction = "both", 
                      trace = FALSE)

# Check the properties of the final model
summary(step.model)

# Save the final model
saveRDS(step.model, "Model_Garden_LCBD_fun_nest.rds")


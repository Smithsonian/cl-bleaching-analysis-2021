R Analysis README file: Emergent properties in the responses of tropical corals to recurrent climate extremes.
Note: Analysis executed in RStudio.

1. The data folder contains all the data files used for this analysis. The data is in csv format and RData format. 

2. Each analysis folder corresponds to a specific bleaching statistical analysis (where bin.score is the response variable):
	* -DHW-region analysis- folder contains the analysis of the GLM model of the form DHW * region.
	* -DHW-year analysis- folder contains the analysis of the GLM model of the form DHW (when analyzing each individual year) and DHW * year by combining all years. 
	* -prior year heat stress- folder contains the analysis of the GLM model of the form DHW(2020) * DHW(previous year) - DHW(previous year) where "previous year" refers to 2016, 2017, and max-DHW(2016,2017). 
	* -prior year heat stress - others- folder contains the analysis of the GLM model of the form DHW(year) * DHW(previous year) - DHW(previous year) for the years 2002, 2016, and 2017.
	
3. Each RData file found in the "datasets" folder, inside each analysis folder, is taken from the "data" folder.

4. All routines follow the same structure. Execute in the order proposed by the file names:
	* 0-initial-setup.R: Loads all relevant packages and files. Creates coordinates matrices for the spatial analysis and any other dataframe that might be necessary. 
	* 1-model-selection.R: Creates a list of spatial structure candidates based on the values given to listw.candidates(). 
						   Uses the function select_weighting_matrix() that uses the ME() function to generate a table that contains information about each candidate including chosen eigenvectors, pvalue, and Moran's I.
						   These tables are stored in the "selected" folder inside each analysis folder.
	* (2,3,4...)-analysis-XXXX.R: Based on the selected candidates (procedure described in the Supplementary Materials) the model response, residuals, and spatial autocorrelation are examined.
	* analysis-functions.R: Defines all functions that are useful for all analysis scripts (for each year). 
	* confusion-matrix-calculation.R: Calculates the confusion matrix and accuracy of a model of interested. Results reported in the paper. 
	* me_modification.R: Modifies the ME() function of the spatialreg package. This can be implemented to the analysis by using trace(ME, edit = TRUE) in RStudio. Replace everything by this new function.

5. To replicate the results in the paper follow this order inside each analysis folder:
              * Run script 0-initial-setup.R
              * Use trace(ME, edit = T) to replace everything for the code found in me_modification.R. This step is fundamental. 
              * Run script 1-model-selection.R to generate the "candidates tables". To load the tables previously calculated (i.e., the tables we obtained when we did the analysis for the paper), load the last chunk of code of the script. The candidates_xxxx lines should always be executed as they are necessary for the next scripts.
              * Run remaining scripts in the order numbered.

							
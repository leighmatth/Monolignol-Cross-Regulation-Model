# Monolignol-Cross-Regulation-Model

This is the code associated with the article: Matthews ML, Wang JP, Sederoff R, Chiang VL, Williams CM, Modeling cross-regulatory influences on monolignol transcripts and proteins under single and combinatorial gene knockdowns in <i>Populus trichocarpa</i>


There are two parts to this code:

1) Identify Edges - These files are used to solve for the indirect influences among the monolignol transcripts and proteins. The main outputs of this code are the B and mue variables, which are used to define the transcript-protein model. The resulting B and mue that were used in the manuscript are stored in ModelDetails.mat. To replicate these results, uncomment the random number generator seed in lignin_modelSML.m. lignin_modelSML.m is the main script. All other functions are called from this script. 

	The code in this folder is built on the Sparse Maximum Likelihood Algorithm code from Cai X, Bazerque JA, Giannakis GB (2013) Inference of Gene Regulatory Networks with Sparse Structural Equation Models Exploiting Genetic Perturbations. PLOS Computational Biology 9(5): e1003068. https://doi.org/10.1371/journal.pcbi.1003068. 

	Changes were made to their base code for our model structure as described in the supporting material of our article (Matthews et al.). The original code from Cai et al., can be found in the supporting material of their paper.
		

2) Monolignol-Transcript-Protein Model - Use identified model (ModelDetails.mat) to estimate response of un-targeted transcripts and proteins for specified knockdowns

	TranscriptProteinModelPredictions.m is the main script, the other two functions are called from this file. ModelPredictions.m returns the predicted transcript and protein abundances from both the new and old models when an vector or matrix containing the desired knockdown abundance of the targeted transcripts are input (Xtarg) and a binary vector or matrix that indicates which transcripts are being targeted (Xmask). KnockdownVisulations.m produces barplots comparing the predicted results to the actual results for a specific gene and targeted knockdown, like those found in Matthews et al., when emulating the transgenic experiments.
	
	
The 2 matlab data files included are:

1) TranscriptProteinAbundances.mat

	This file contains the following variables:
	
	i. Yallexps: Table containing the normalized 20 transcript and 20 protein abundances for each of the 225 transgenic and wildtype experiments
	
	ii. Yallexps_imputed: Table where the missing values of Yallexps are imputed as described in Matthews et al.
	
	iii. Ymeanlines_imputed: Table of the average abundances over the replicates for each line of each transgenic construct in Yallexps_imputed.
	
	iv. Xallexps_mask: Binary table denoting which transcripts were targeted (1) corresponding to Yallexps and Yallexps_imputed.
	
	v. Xmeanlines_mask: Binary table denoting which transcripts were targeted (1) corresponding to Ymeanlines_imputed.
	
	vi. Ywt: Table containing the transcript and protein abundances from only the wildtypes.
	
	vii. Ywt_avg: Table containing the average wildtype abundance for each transcript and protein.
	

2) ModelDetails.mat

	This file contains the following variables:
	
	i. B: 40 x 40 matrix containing the inferred cross-regulatory influences between the monolignol transcripts and proteins
	
	ii. mue: 40 x 1 vector containing the constant terms for the model
	
	iii. Boldmodel: 40 x 40 matrix containing the one-to-one transcript to protein relationships used previously (Wang et al., 2018, Improving wood properties for wood utilization through multi-omics integration in lignin biosynthesis).
	
	iv. GeneLabels: 1 x 40 cell containing the list of the transcript and protein names in the order they appear in B, Boldmodel, mue.

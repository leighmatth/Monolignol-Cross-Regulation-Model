clear
clc

load ../TranscriptProteinAbundances.mat

%% To emulate experiments in paper
Ytable=Yallexps_imputed;
Xtable=Xallexps_mask;
[N,M]=size(Ytable);
Xmask=Xtable{:,:}';
Xtarg=Xmask.*Yallexps_imputed{:,:}';
Xfulltrans=Yallexps_imputed{:,1:M/2}';

Ypred=Model_Predictions(Xmask,Xtarg,Xfulltrans);

GeneLabels=Ytable.Properties.VariableNames;
Experiments=Ytable.Properties.RowNames;
[Gc, Const]=findgroups(cellfun(@(x) x(6:8),Experiments,'UniformOutput',false)); %Group experiments into their different constructs (targeted knockdowns)

potential_genes_of_interest=cellfun(@(x) x(2:end),GeneLabels(1:20),'UniformOutput',0); % List of potential genes of interest for visualization
potential_knockdowns_of_interest=Const; % List of potential constructs/targeted knockdowns for visualization

%% Create bar plots showing new and old model estimates (as seen in paper)
gene_of_interest='HCT1'; % See potential_genes_of_interest for other options
knockdown_of_interest='i29'; % See potential_knockdowns_of_interest for other options
legend_flag=1; %1 for legend on plots, 0 for no legend

Knockdown_Visualizations(gene_of_interest,knockdown_of_interest,Ytable,Xtable,Ypred,legend_flag)


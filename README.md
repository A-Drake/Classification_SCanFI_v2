Associated code for the manuscript: "Classifying boreal landbirds to facilitate multi-species approaches 
to national-level conservation and forest degradation analyses"

by Anna Drake, Elly Knight, Lisa Venier, Mannfred Boehm, Luc Guindon, David Correia, Diana Stralberg

Forest landbirds have high potential as indicators of changes in forest state including forest degradation. 
Multi-species responses may better-inform degradation assessments so long as species are grouped 
appropriately. We used Canada-wide landbird survey data and remotely sensed land cover information to 
demonstrate a workflow for classifying 117 landbirds according to local (within 200 m) and regional (within
2 km) habitat preference across 16 ecological regions. Using species distribution models, we quantified 
preference as habitat use relative to availability along three major axes of forest variability: percent 
cover, proportion conifer, and forest age. We describe this pattern as a selection curve bounded by the 
domain of a given covariate. Selection curves were generated for every species-ecoregion in our dataset. 

Describing these curves involved three steps: 
(1) quantifying predicted species use of covariate values (n=32) and the availability of these values (Script A-B), 
(2) calculating selection ratios from proportional use:availability (Script B-C) , and 
(3) fitting the resulting ratios to a selection curve and associated confidence intervals (Script D-E).

Jaccard Similarity analysis and full vs reduced model performance correspond to Script F & G

See manuscript for input data sources.

Analyses were completed in R. Those labelled "Compute Canada" were scripted to run using HPC. The associated 
shell script is provided at the bottom of these R. scripts

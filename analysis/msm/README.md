### data
        data acquired by following a notebook as above for all Isoleucines with Charmm36m/TIP3P and Amber (ff99SB*-nmr1/TIP3P) force fields

### Building_MSM_Ile_sidechain.ipynb
	Example jupyter notebook that can build a Markov Model for the side chain motion from "scratch"
	  (i.e. requires the trajectory and topology as input)
	
	If you are trying to build MSMs for a new system and want to use an example notebook as starting point, 
	this is the notebook you want to look at.

### Markov2Ct.py
	Helper functions used by notebook above

### custom_plots.py
	Plotting functions used by notebook above

### make_si_figs_msms.ipynb
	notebook allowing to reproduce the figuresin the si of the paper (based on data saved in data)

### make_paper_fig_comparison_msm_romance.ipynb
	allows to reproduce the figure shown in the article, as well as a simplified version of the figure (e.g. for presentations)

### make_paper_fig_comparison_msm_explicit_models.ipynb
	allows to reproduce the figure shown in the SI, comparing to explicit models

### populations_rotamers.ipynb
	allows to reproduce the figure in main text

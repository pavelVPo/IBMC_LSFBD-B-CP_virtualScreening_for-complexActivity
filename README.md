# IBMC_LSFBD-B-CP_virtualScreening_for-complexActivity

**The following text and code is a draft**

This repository contains the code and by extension some data to conduct the virtual screening for the compounds having cytotoxicity against the cell lines,
which could be used to model cervical cancer, affecting Shh pathway, affecting MMPs expression in an unsual manner; and
associated with the antiviral action; which could be of interest to our collegues from B&amp;CPP (https://en.ibmc.msk.ru/departments?view=article&id=41:laboratory-of-protein-biochemistry-and-pathology&catid=10:data).

## Preliminary search for the compounds having chances to have the desired activity profile using TDV (Tool for Diversity Visualization, https://way2drug.com/tdv/)

By occassion the catalougue of the fluorescent dyes (https://www.lumiprobe.com/catalog/reactive-fluorophores-quenchers) was screened using TDV to find smth having the same scaffolds as known drugs or interesting known bioactives. Three reagents were identified on the level of generic scaffolds: one is modified solubilizable cholesterol, which could be functionalized further using click-chemistry; two other reagents have the same generic scaffold as shikoin, a natural compound, which is cytotoxic against several cancer cell lines.

Interestingly, modified cholesterol fits the desired activity profile quite well, probably, except the cytotoxicity:

- It is well known that many of cholesterol's derivatives are anticancer drugs, and strictly speaking, cytotoxic properties of the exogenic cholesterol itself (not in complex with other molecules) is hard to define due to its insolubility. However, inclusion of the cholesterol in various nanoformulations, testifies in favour of its general safety in this context (** which is not directly related to the topic of cholesterol levels in blood, atherosclerosis, etc. **). Still, how solubilized cholesterol alone exactly affects the viability of the cell lines, which are models of the cervical cancer; is an open question according to the preliminary info search on this topic.

- Cholesterol is known factor affecting Shh pathway (https://www.pnas.org/doi/10.1073/pnas.0600124103).

- MMP expression probably could be affected by the cholesterol (https://www.researchgate.net/publication/320643835_Lipid_metabolism_fattens_up_hedgehog_signaling, https://pubmed.ncbi.nlm.nih.gov/17643435/).

- Cholesterol content to the large extent defines the membrane's biophysical properties, thus, it could influence the way how cell interacts with the viral particles.

## PASS-based virtual screening

- Virtual screening is a computational procedure, which allows one to select the most promising compounds for the biological evaluation among many available compounds.

- PASS (Prediction Activity Spectra for Substances, https://way2drug.com/PassOnline/, https://academic.oup.com/bioinformatics/article/16/8/747/190470) is the tool having self-explanatory title.

### For what to screen given the complex nature of the desired activity profile?

Cytotoxic compounds are numerous and could act on the multiple targets in the cell. Shh-pathway contains Smo-protein, which is relatively new pharmacological target associated somehow with all the desired biological phenomena. Thus, at first, it will be rational to find novel Smo-inhibitor, which is cytotoxic against the one or more of the cell lines used to model the cervical cancer: C33A, Ca-Ski, SiHa.

### Training data

ChEMBL v36 (https://www.ebi.ac.uk/chembl/) is used as the source of the training data (https://link.springer.com/article/10.1186/s13321-025-00963-z). Records containing IC50 values are considered.

### Data to screen

For this study it was proposed to screen the ChemRar's library, namely the file containing 300k diverse structures was screened (https://mol.chemrar.ru/diversity-libraries , '300K Разнообразная библиотека (алгоритм кластеризации Bemis-Murcko)' ).

### Propensity stage of the virtual screening

Using PASS-software it is possible to build classifiers in the framework of (Q)SAR / (Q)SPR (https://pubmed.ncbi.nlm.nih.gov/26754147/) methodology, which is a good thing for the folks using PASS, since almost every activity / property of the chemical compounds is determined by its chemical structure (surprisingly, exceptions exist). Primordial PASS was designed to propose the most promising directions for the biological evaluation of chemical compounds based on the extremely diverse data on structures and activities of newly synthesized and studied chemical compounds ([the appropriate reference will be added](http://www.russchembull.ru/rus/objects/papcat-5167.pdf)). This allows to use as much data as possible, but, at the same time, means that endpoints used could be not the established, standardized ones. So, this idea could be extended to the case, when PASS is trained to predict not an activity, but rather a propensity of the chemical compound to be evaluated for the certain activity (https://en.wikipedia.org/wiki/Propensity_score_matching , https://www.tandfonline.com/doi/abs/10.1080/1062936X.2019.1665580). Why to do that, if it is possible to immediately predict the activity instead of propensity? Because it was shown that it allows to increase the overall success of the virtual screening (https://www.tandfonline.com/doi/abs/10.1080/1062936X.2019.1665580 , https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2018.00133/full).

Thus, this part of the study is basically about predicting whether compound is likely to be evaluated for certain biological activity or not.

And 'activity' here means the ability to sufficiently inhibit some biological target.

Targets, considered at the moment are:

| tid    | pref_name          | chembl_id     |
| 103218 | Ca-Ski             | CHEMBL1075403 |
| 106482 | C-33-A             | CHEMBL2366313 |
| 80472  | SiHa               | CHEMBL612542  |
| 101400 | Protein smoothened | CHEMBL5971    |
|--------|--------------------|---------------|


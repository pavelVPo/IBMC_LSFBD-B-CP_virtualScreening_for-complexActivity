# IBMC_LSFBD-B-CP_virtualScreening_for-complexActivity

This repository contains the code and by extension some data to conduct the virtual screening for the compounds having cytotoxicity against the cell lines,
which could be used to model cervical cancer, affecting Shh pathway, affecting MMPs expression in an unsual manner; and
associated with the antiviral action; which could be of interest to our collegues from B&amp;CPP (https://en.ibmc.msk.ru/departments?view=article&id=41:laboratory-of-protein-biochemistry-and-pathology&catid=10:data).

## Preliminary search for the compounds having chances to have the desired activity profile using TDV (Tool for Diversity Visualization, https://way2drug.com/tdv/)

By occassion the catalougue of the fluorescent dyes (https://www.lumiprobe.com/catalog/reactive-fluorophores-quenchers) was screened using TDV to find smth having the same scaffolds as known drugs or interesting known bioactives. Three reagents were identified on the level of generic scaffolds: one is modified solubilizable cholesterol, which could be functionalized further using click-chemistry; two other reagents have the same generic scaffold as shikoin, a natural compound, which is cytotoxic against several cancer cell lines.

Interestingly, modified cholesterol fits the desired activity profile quite well, probably, except the cytotoxicity:

- It is well known that many of cholesterol's derivatives are anticancer drugs, and strictly speaking, cytotoxic properties of the cholesterol itself is hard to define due to its insolubility. However, inclusion of the cholesterol in various nanoformulations, testifies in favour of its general safety in this context (** which is not directly related to the topic of cholesterol levels in blood, atherosclerosis, etc. **). Still, how solubilized cholesterol alone exactly affects the viability of the cell lines, which are models of the cervical cancer; is an open question according to the preliminary info search on this topic.

- Cholesterol is known factor affecting Shh pathway (https://www.pnas.org/doi/10.1073/pnas.0600124103).

- MMP expression probably could be affected by the cholesterol (https://www.researchgate.net/publication/320643835_Lipid_metabolism_fattens_up_hedgehog_signaling).

- Cholesterol content to the large extent defines the membrane's biophysical properties, thus, it could influence the way how cell interacts with the viral particles.

## PASS-based virtual screening

- Virtual screening is a computational procedure, which allows one to select the most promising compounds for the biological evaluation among many available compounds.

- PASS (Prediction Activity Spectra for Substances, https://way2drug.com/PassOnline/, https://academic.oup.com/bioinformatics/article/16/8/747/190470) is the tool having self-explanatory title.

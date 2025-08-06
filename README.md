# Local, Regional, and Coastwide effects of interactions between storms and relative position in tidal frame on Louisiana coastal marshes

#### Authors: Donald R. Schoolmaster Jr., Camille L. Stagg, Gregg A. Snedden, and Brady R. Couvillion

#### Point of contact: Donald R. Schoolmaster Jr. (schoolmasterd@usgs.gov)
#### Repository Type:  _R_ script supporting publication
#### Year of Origin:   2025 (original publication)
#### Year of Version:  2025
#### Digital Object Identifier (DOI):TBD
#### USGS Information Product Data System (IPDS) no.: IP-180480

***

_Suggested Citation:_

This draft manuscript is distributed solely for purposes of scientific peer review. Its content is deliberative and predecisional, so it must not be disclosed or released by reviewers. Because the manuscript has not yet been approved for publication by the U.S. Geological Survey (USGS), it does not represent any official USGS finding or policy.

_Authors' [ORCID](https://orcid.org) nos.:_ \
DRS: 0000-0003-0910-4458\
CLS: 0000-0002-1125-7253\
GAS: 0000-0001-7821-3709\
BRC: 0000-0001-5323-1687
***
***


## Biological Subject Area and Programming Background

This code assumes the user is familiar with R, JAGS and uses the rjags() package.

## Software Version Details
The code runs using R version 4.4.2.

## File Details
The main file for this analysis is "compile_and_sample.R", it calls "land_change_analysis_v2.R" file. The output is saved as RMD files to the /Output/Results folder. These results are called by "create_figures_v2.R".

## Code Run Time
The code in the file compile_and_sample.R can many hours to run and some of the resulting files are over 100MB (too large to be included in the repository). The other scripts take less than a minute to run.

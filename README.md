This repository holds the data and source code for the following manuscript:

[Minikel EV, Nelson MR. **Human genetic evidence enriched for side effects of approved drugs.** _medRxiv._ 2023 Dec 13;2023.12.12.23299869.](https://doi.org/10.1101/2023.12.12.23299869)

Here, you can:

+ Run the source code to reproduce the figures from the input datasets. Just say `Rscript `[`src/drug_safety_analysis.R`](/src/drug_safety_analysis.R), noting the dependencies at the top of the script. It completes in about 16 minutes on a 2021 MacBook Pro. The script reproduces Figures 1-4 and S1-S3, Tables S1-S23, and stats_for_text.txt, all of which you can find in [display_items](/display_items).
+ If you're curious, you can also browse the source code for other scripts that prepared this releasable analytical dataset, in [src](/src). These scripts require some inputs that are either too large for GitHub, and/or not approved for public release, thus, you will not be able to successfully run them after cloning the repository; they are provided simply for reference in case you want to see what we did.
+ Browse the input datasets in [data](/data) and the matrix of all possible drug-side effect pairs in [intermediate](/intermediate).

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.


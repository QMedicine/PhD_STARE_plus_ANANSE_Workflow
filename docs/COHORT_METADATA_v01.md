# Assay-specific cohort metadata

The in-house Hi-C cohort was supplied by the user as female C57BL/6 murine HSCs from aged >24-month and young approximately 3-month mice. These descriptors must not be assigned automatically from the Old/Young labels to every assay.

The archived SRA/sample-information tables identify the selected public RNA/ATAC HSC samples as 10-week and 20-month mice. The primary publication reports female C57BL/6 mice and confirms the donor ages: [Itokawa et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9110722/), DOI [10.1038/s41467-022-30440-2](https://doi.org/10.1038/s41467-022-30440-2).

The matched sorting definition in the archived RNA metadata is IL-7− lineage− Sca-1+ c-kit+ CD34− CD150+ CD48− Flt3−. Hi-C and public RNA/ATAC are independent cohorts; matching replicate numbers across assays do not imply matched mice. The reconstruction is condition/state-level multimodal integration, with the age ranges recorded per assay.

Generic datasets supply their own cell_type, condition, cell_state and species fields. Run descriptions and cohort metadata should be provided explicitly; no biological age or sex is inferred from filenames or condition words in the executable workflow.

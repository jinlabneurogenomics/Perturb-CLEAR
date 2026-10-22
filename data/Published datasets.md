# Published datasets

Original published source data are not redistributed with this repository. We preserved the published sample and animal identifiers during data processing. Please follow the guidance below to get access to the associated data.

## Sorensen et al. 2026

Sorensen neurons are the primary paired RNA and morphology cohort. Start from the [publication](https://doi.org/10.1038/s41586-026-10424-8) and the Allen Institute [file manifest](https://raw.githubusercontent.com/AllenInstitute/exc_vis_manuscript/main/exc_patchseq_file_manifest.csv). The manifest provides Brain Image Library links for morphology reconstructions and NeMO links for transcriptomic files. [Supplementary Table 3](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-10424-8/MediaObjects/41586_2026_10424_MOESM5_ESM.xlsx) provides cell metadata, and [Supplementary Table 5](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-10424-8/MediaObjects/41586_2026_10424_MOESM7_ESM.xlsx) provides the measurement dictionary. The morphology workflow used the upright SWC reconstructions.

## Scala et al. 2021

Scala neurons are used as a sensitivity cohort. See the [publication](https://doi.org/10.1038/s41586-020-2907-3), the public [cell metadata](https://raw.githubusercontent.com/berenslab/mini-atlas/master/data/m1_patchseq_meta_data.csv), and the [Brain Image Library excitatory-neuron directory](https://download.brainimagelibrary.org/3a/88/3a88a7687ab66069/excitatory/).

## Gouwens et al. 2019

Gouwens neurons are used only for morphology normalization. See the [publication](https://doi.org/10.1038/s41593-019-0417-0). Reconstruction files were obtained through the [Allen Brain Atlas API](https://help.brain-map.org/display/api/Allen+Brain+Atlas+API), and [Supplementary Dataset 3](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41593-019-0417-0/MediaObjects/41593_2019_417_MOESM5_ESM.xlsx) provides the associated morphology metadata. These neurons are not used as a validation cohort for paired RNA-morphology mapping.
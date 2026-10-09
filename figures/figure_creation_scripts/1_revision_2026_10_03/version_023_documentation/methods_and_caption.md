# Figure 1 revision, version 023

## Methods text

The new panels were generated directly from `Broad_8K_Antibiotic_Activity_Dataset_Filtered_.csv`. The earlier Figure 2d property table was not used to define compounds or activity groups. Molecular structures were parsed with RDKit 2025.09.1 and identified by canonical isomeric SMILES. Stereochemistry and disconnected fragments were retained. No additional protonation, tautomer, or salt standardization was performed. The supplied file contains 7,997 records representing 7,985 unique structures.

Activity groups use the supplied binary labels for wild type (WT), Δ40, Δ3, and Δ6. A compound enters a strain group when its corresponding binary activity label is 1. Repeated records are combined by canonical isomeric structure, with positive calls combined by logical OR within each strain. Missing continuous inhibition measurements do not exclude compounds with an explicit binary hit label.

The groups contain 8 WT hits, 5 Δ40 hits, 39 Δ3 hits, and 124 Δ6 hits. This comparison uses all hits for each strain, so the same compound can appear in more than one group. The WT overlaps are 4 compounds for Δ40, 8 for Δ3, and 8 for Δ6. There remain 116 Δ6 hits inactive in WT, including the three compounds with missing continuous WT inhibition measurements. The plotted Δ6 group contains all 124 Δ6 hits. There is no inactive comparison group. Across the four strains, the plotted compounds represent 125 unique molecular structures.

Mean percentage growth inhibition was taken from the corresponding continuous measurement columns. Each strain's violin used all available measurements for that strain, with repeated measurements averaged within a canonical structure. The WT distribution contained 7,671 measured structures and each mutant distribution contained 7,985 measured structures. Missing continuous values were not imputed from binary activity labels. The displayed n = 7,994 is the screen annotation. Density sample sizes are reported separately above.

Gaussian kernel density estimates used a common bandwidth of two percentage points. Violin half-widths use a shared base-10 logarithmic density transformation to make low-density tails visible. For probability density $p(x)$ per percentage point, the half-width is $w(x) = 0.405 \max[\log_{10}(p(x)/p_0), 0] / \log_{10}(p_{\max}/p_0)$, where $p_0 = 10^{-5}$ and $p_{\max}$ is the largest density over all four strains. Density values below the common floor have zero displayed width. Each violin is mirrored around its strain centre, and widths are not rescaled separately by strain. The inhibition axis remains linear. Measurements below zero and above 100% were retained. Black bars indicate interquartile ranges and white points indicate medians. The dashed line marks 80% inhibition. Annotations report eight WT hits, five ΔOMP40 hits, 39 ΔE-pump3 hits, and 124 ΔE-pump6 hits. Hit percentages remain available in the analysis table using the actual measured structure counts, but are not displayed on panel f.

Descriptors were calculated with RDKit 2025.09.1. IPC was calculated with `GraphDescriptors.Ipc(mol, avg=False)` and transformed as $\log_{10}(\mathrm{IPC})$. All IPC values were finite and positive, so no pseudocount was used. IPC is an information-content descriptor calculated from coefficients of the molecular graph's characteristic polynomial ([RDKit definition](https://www.rdkit.org/docs/source/rdkit.Chem.GraphDescriptors.html#rdkit.Chem.GraphDescriptors.Ipc)). LogP was calculated with `Crippen.MolLogP` ([RDKit definition](https://www.rdkit.org/docs/source/rdkit.Chem.Crippen.html)). Aromatic rings were counted with `rdMolDescriptors.CalcNumAromaticRings` ([RDKit definition](https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html)).

Boxes extend from the first to the third quartile, and central lines indicate medians. Whiskers reach the most extreme observations within 1.5 interquartile ranges of the box. Every compound is shown as an individual point for every group, including WT and Δ6. The same marker size, opacity, edge width, and horizontal jitter rule are used throughout. Jitter depends only on strain and molecular structure, so the same compound has the same horizontal offset across the three properties. Values are not jittered vertically. Box-plot outlier markers are suppressed because those observations are already included in the complete point layer. Points can overlap, particularly for integer aromatic ring counts. No points are sampled, imputed, or clipped by the displayed axes.

The property panels are descriptive, and no independent-group significance tests are reported. All displayed values fall inside the axes. The existing strain colours are retained. Panel g occupies the same region as the earlier property comparison and now contains one row of three plots.

Compound counts were calculated from the supplied filtered 8K, 39K, ECBD, and CO-ADD files. Each canonical isomeric structure contributed once per dataset and strain. Existing binary calls were combined by logical OR within each structure. Missing binary calls were retained as missing. Nine records with unparseable structures were excluded across 39K and CO-ADD. None carried a positive call in the displayed assays. Counts represent the supplied dataset annotations, including the ECBD hit annotation, and do not assume identical experimental conditions across libraries.

## Caption additions

**d**, Updated numbers of unique hit structures in each supplied dataset. Each structure is counted once per strain within each dataset. Mutant counts include structures that are also WT hits. **f**, Violin distributions of mean growth inhibition. Widths use a common logarithmic probability density scale with a floor of $10^{-5}$ per percentage point and a Gaussian bandwidth of two percentage points. The inhibition axis is linear. Screen annotation n = 7,994. Density sample sizes are 7,671 WT structures and 7,985 structures per mutant. Bars and points indicate interquartile ranges and medians. The dashed line marks 80% inhibition, and annotations give hit counts. **g**, Physicochemical properties of hits in WT, Δ40, Δ3, and Δ6. The panels show $\log_{10}(\mathrm{IPC})$, logP, and number of aromatic rings. IPC denotes the information content descriptor. Each strain group contains all compounds with a positive binary activity label for that strain, with sample sizes shown below the groups. Shared hits appear in more than one group. Boxes indicate interquartile ranges, central lines indicate medians, and whiskers extend to the most extreme observations within 1.5 interquartile ranges. Every compound is shown as a point using the same plotting rule in all groups. Horizontal jitter aids visibility, and points may overlap. Green, orange, red, and purple denote WT, Δ40, Δ3, and Δ6, respectively. **h, i**, Original chemical-space maps, retained with their original coordinates and membership.

## Activity-label reconciliation

Three Δ6-only hits have missing continuous WT inhibition measurements but have explicit WT activity labels of 0 and Δ6 activity labels of 1. All three are included as Δ6-only hits.

The Δ6 violin and full-library count panel both include all 124 Δ6 hits.

The original replicate panels retain their annotation of n = 7,994 because their artwork was not regenerated. Panel f uses the same user-specified screen annotation, while its actual available-data counts are documented above. Chemical-space maps use the original coordinates and membership. Updated counts use the current supplied files.

## Unique-hit counts

| Dataset | WT | ΔOMP40 | ΔE-pump3 | ΔE-pump6 | ΔE-pump5 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Internal 8K | 8 | 5 | 39 | 124 | |
| Internal 39K | 8 | | | | |
| ECBD | 1 | | | | |
| CO-ADD | 7 | | | | 75 |

CO-ADD's 75 Δ5 hit structures correspond to 143 positive source records. Repeated records are not counted as different molecules.

## Reproduction

The current panel g is reproduced by `isolated_panel_g_2026_10_04/methods/reproduce_isolated_properties.py`, which reads the current 8K CSV directly using its absolute path. `validate_isolated_properties.py` independently checks source membership, descriptors, point counts, and SVG output. The earlier Figure 2d property list is not an input. The three original isolated SVGs remain available.

`methods/prepare_three_property_panel.py` arranges those same data and plots into the 257 by 164 point panel g region. It records an audit of all 528 compound points across the three properties. `methods/integrate_three_properties_in_illustrator.jsx` replaces the property panel in version 022 and saves `Figure1_v023_three_properties.ai`. `methods/export_three_property_figure.py` creates the vector PDF, PNG, and 300 dpi review PDF from the saved native Illustrator file and checks the artwork outside panel g against version 022. The original supplied Illustrator file remains unchanged.

The isolated-panel tables contain 125 unique hit structures and 176 strain-specific observations. These observations comprise 8 WT, 5 Δ40, 39 Δ3, and 124 Δ6 hits. Source identifiers, activity labels, calculated descriptors, overlapping membership, and plotting audits are included. The historical `tables/8K_Delta6_only_116_compounds.csv` remains available to document the 116 WT-inactive Δ6 hits, but the current Δ6 property group contains all 124 hits.

`methods/reproduce_revision.py` and its companion scripts reproduce the earlier count and inhibition panels and the historical six-property comparison. Running that script alone does not generate the current panel g. Earlier property summaries and independent-group comparisons refer to the historical WT, Δ6-only, and inactive groups and are not statistical analyses of version 023.

Scripts require no command-line arguments and write products to the configured absolute revision directory. The original replicate panels and chemical-space maps are retained, with the source-data distinctions documented above.

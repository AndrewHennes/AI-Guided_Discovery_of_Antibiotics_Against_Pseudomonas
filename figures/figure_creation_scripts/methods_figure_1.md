# Figure 1 revision, version 025

## Methods text

The new panels were generated directly from `Broad_8K_Antibiotic_Activity_Dataset_Filtered_.csv`. The earlier Figure 2d property table was not used to define compounds or activity groups. Molecular structures were parsed with RDKit 2025.09.5 and identified by canonical isomeric SMILES. Stereochemistry and disconnected fragments were retained. No additional protonation, tautomer, or salt standardization was performed. The supplied file contains 7,997 records representing 7,985 unique structures.

Property groups use the supplied binary activity labels, with positive calls combined by logical OR within each canonical structure and strain. Missing continuous measurements do not exclude compounds from these groups.

Each strain group contains every compound with a positive binary activity label in that strain. The groups contain 8 WT, 5 Δ40, 39 Δ3, and 124 Δ6 structures. Compounds active in multiple strains appear in every corresponding group. No exclusions are applied based on activity in another strain. All eight WT hits are also Δ3 and Δ6 hits, and all 39 Δ3 hits are Δ6 hits. Four WT hits are also Δ40 hits. The four groups contain 176 observations representing 125 unique structures. Figure labels use only WT, Δ40, Δ3, and Δ6 with their sample sizes. No inactive comparison group is shown.

Mean percentage growth inhibition was taken from the corresponding continuous measurement columns. Each strain's violin used all available measurements for that strain, with repeated measurements averaged within a canonical structure. The WT distribution contained 7,671 measured structures and each mutant distribution contained 7,985 measured structures. Missing continuous values were not imputed from binary activity labels. The displayed n = 7,994 is the screen annotation. Density sample sizes are reported separately above.

Gaussian kernel density estimates used a common bandwidth of two percentage points. Violin half-widths use a shared base-10 logarithmic density transformation to make low-density tails visible. For probability density $p(x)$ per percentage point, the half-width is $w(x) = 0.405 \max[\log_{10}(p(x)/p_0), 0] / \log_{10}(p_{\max}/p_0)$, where $p_0 = 10^{-5}$ and $p_{\max}$ is the largest density over all four strains. Density values below the common floor have zero displayed width. Each violin is mirrored around its strain centre, and widths are not rescaled separately by strain. The inhibition axis remains linear. Measurements below zero and above 100% were retained. Black bars indicate interquartile ranges and white points indicate medians. The dashed line marks 80% inhibition. Annotations report eight WT hits, five ΔOMP40 hits, 39 ΔE-pump3 hits, and 124 ΔE-pump6 hits. Hit percentages remain available in the analysis table using the actual measured structure counts, but are not displayed on panel f.

Topological polar surface area (TPSA) was calculated with `rdMolDescriptors.CalcTPSA(mol, includeSandP=False)`. LogP was calculated with `Crippen.MolLogP`. Amine count is the total number of primary, secondary, and tertiary amine nitrogen atoms, including protonated forms and amines attached to aromatic carbon atoms. Aromatic ring nitrogens, amides, thioamides, sulfonamides, amidines, guanidines, cyanamides, imines, nitriles, nitro groups, amine oxides, and quaternary ammonium groups are excluded. The exact query is provided in `shared/figure1_amine_descriptor.py`, along with 33 functional-group checks. This descriptor counts structural amine sites and does not predict their protonation state at assay pH. Separate primary, secondary, and tertiary amine counts are included in the data table.

Boxes show the first and third quartiles, with a line at the median. Whiskers reach the most extreme observations within 1.5 interquartile ranges of the box. Every compound is overlaid with identical marker styling and a deterministic horizontal jitter rule. Box-plot outlier markers are suppressed because those observations are already included in the complete point layer. No vertical jitter or point sampling is used. Values are not clipped. Observed ranges are 0 to 206.07 square ångströms for TPSA, 0 to 3 for amine count, and approximately 0.230 to 7.634 for logP. The axes span −8 to 220, −0.2 to 3.25, and 0 to 8, respectively. Small margins below zero for TPSA and amine count keep zero-valued points fully visible. The original green, orange, red, and purple strain colours are retained.

Group membership was checked independently using sets of positively labelled source records. Descriptor values for all 125 structures were checked against the verified version 025 property table. The plotting audit checks all 528 points across the three panels.

Compound counts were calculated from the supplied filtered 8K, 39K, ECBD, and CO-ADD files. Each canonical isomeric structure contributed once per dataset and strain. Existing binary calls were combined by logical OR within each structure. Missing binary calls were retained as missing. Nine records with unparseable structures were excluded across 39K and CO-ADD. None carried a positive call in the displayed assays. Counts represent the supplied dataset annotations, including the ECBD hit annotation, and do not assume identical experimental conditions across libraries.

## Caption additions

**d**, Updated numbers of unique hit structures in each supplied dataset. Each structure is counted once per strain within each dataset. Mutant counts include structures that are also WT hits. **f**, Violin distributions of mean growth inhibition. Widths use a common logarithmic probability density scale with a floor of $10^{-5}$ per percentage point and a Gaussian bandwidth of two percentage points. The inhibition axis is linear. Screen annotation n = 7,994. Density sample sizes are 7,671 WT structures and 7,985 structures per mutant. Bars and points indicate interquartile ranges and medians. The dashed line marks 80% inhibition, and annotations give hit counts. **g**, Physicochemical properties of active compounds in WT, Δ40, Δ3, and Δ6. The groups contain all compounds with a positive binary activity label in the corresponding strain, giving n = 8, 5, 39, and 124, respectively. Compounds active in multiple strains appear in each corresponding group. Activity is defined by the supplied binary labels at 80% growth inhibition. The panels show topological polar surface area (TPSA, $\mathrm{\AA}^2$), number of primary, secondary, and tertiary amines combined, and Wildman–Crippen logP. Amine counts include protonated forms and exclude amides and aromatic ring nitrogens. Boxes show interquartile ranges, central lines show medians, and whiskers extend to the most extreme observations within 1.5 interquartile ranges of the box. Every compound is shown as a point using the same horizontal jitter rule and marker style in all groups. Points may overlap, especially for integer amine counts. All values are included within the displayed axes. Green, orange, red, and purple denote WT, Δ40, Δ3, and Δ6, respectively. **h, i**, Original chemical-space maps, retained with their original coordinates and membership.

## Activity-label reconciliation

The historical set of 116 WT-inactive Δ6 hits includes three compounds with missing continuous WT inhibition measurements and explicit WT activity labels of 0. They remain eligible under the binary-label grouping rules. The current panel g includes every Δ6 hit, giving 124 compounds regardless of its WT or Δ3 activity.

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

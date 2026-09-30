# Solvent accessibility of the candidate epitope

Computed output of `analysis/03_solvent_accessibility.py`. Do not hand-edit.

Normalisation maxima are from Table 1 of Tien MZ, Meyer AG, Sydykova DK,
Spielman SJ, Wilke CO (2013), "Maximum Allowed Solvent Accessibilites of
Residues in Proteins", PLOS ONE 8(11): e80635,
doi:10.1371/journal.pone.0080635 — transcribed from the journal's
manuscript XML on 30 September 2026.

The classification cutoffs (RSA >= 0.25 exposed, <= 0.05 buried) are
**conventional, not physical constants**. Nothing changes in the molecule
at 0.25. Residues near a cutoff are flagged BORDERLINE.

```
========================================================================
SOLVENT ACCESSIBILITY OF THE CANDIDATE EPITOPE
========================================================================

Input:  data/structures/6aru_receptor_only.pdb  (Fab REMOVED)
Method: Shrake-Rupley probe rolling, Bio.PDB.SASA.ShrakeRupley
Norm:   Tien et al. 2013, PLOS ONE 8(11):e80635, Table 1
        theoretical column for headline values,
        empirical column as a sensitivity check
Numbering: UniProt = PDB + 24 (established empirically by step 02)

The Fab is removed deliberately. Cetuximab sits on the surface we are
measuring, so computing on the complex would report our epitope as
buried when it is only covered by an antibody absent from our assay.

Computed SASA for 609 residues in the receptor chain.

1. Every residue in 415-466

   RSA is the number to read. 'anchor' marks our eight candidates.
   'BORDERLINE' means within 0.05 of a classification cutoff, so the
   label should not be read as precise.

   | UniProt | aa | PDB# | SASA A^2 | RSA | class | flags |
   |---|---|---|---|---|---|---|
   | 415 | T | 391 |    76.3 | 0.444 | exposed |  |
   | 416 | D | 392 |    22.5 | 0.117 | partial | ANCHOR |
   | 417 | L | 393 |     1.2 | 0.006 | buried | BORDERLINE |
   | 418 | H | 394 |     7.2 | 0.032 | buried | ANCHOR, BORDERLINE |
   | 419 | A | 395 |     0.0 | 0.000 | buried |  |
   | 420 | F | 396 |     0.0 | 0.000 | buried |  |
   | 421 | E | 397 |    46.8 | 0.210 | partial | ANCHOR, BORDERLINE |
   | 422 | N | 398 |    70.2 | 0.360 | exposed |  |
   | 423 | L | 399 |     1.1 | 0.005 | buried | BORDERLINE |
   | 424 | E | 400 |    95.7 | 0.429 | exposed | ANCHOR |
   | 425 | I | 401 |    18.1 | 0.092 | partial | BORDERLINE |
   | 426 | I | 402 |     0.0 | 0.000 | buried |  |
   | 427 | R | 403 |    67.3 | 0.246 | partial | BORDERLINE, NORM-SENSITIVE(exposed if empirical) |
   | 428 | G | 404 |     0.0 | 0.000 | buried |  |
   | 429 | R | 405 |   165.5 | 0.604 | exposed |  |
   | 430 | T | 406 |    69.6 | 0.405 | exposed |  |
   | 431 | K | 407 |    51.2 | 0.217 | partial | BORDERLINE |
   | 432 | Q | 408 |    36.5 | 0.162 | partial |  |
   | 433 | H | 409 |   145.1 | 0.648 | exposed | ANCHOR |
   | 434 | G | 410 |    54.1 | 0.520 | exposed |  |
   | 435 | Q | 411 |    89.9 | 0.399 | exposed |  |
   | 436 | F | 412 |    40.9 | 0.171 | partial |  |
   | 437 | S | 413 |     1.2 | 0.008 | buried | BORDERLINE |
   | 438 | L | 414 |     1.2 | 0.006 | buried | BORDERLINE |
   | 439 | A | 415 |     7.1 | 0.055 | partial | BORDERLINE |
   | 440 | V | 416 |     1.2 | 0.007 | buried | BORDERLINE |
   | 441 | V | 417 |    45.9 | 0.264 | exposed | BORDERLINE |
   | 442 | S | 418 |    67.0 | 0.432 | exposed |  |
   | 443 | L | 419 |     7.2 | 0.036 | buried | BORDERLINE |
   | 444 | N | 420 |    46.6 | 0.239 | partial | BORDERLINE |
   | 445 | I | 421 |     1.1 | 0.006 | buried | BORDERLINE |
   | 446 | T | 422 |    33.8 | 0.196 | partial |  |
   | 447 | S | 423 |     8.6 | 0.055 | partial | BORDERLINE |
   | 448 | L | 424 |     4.6 | 0.023 | buried | BORDERLINE |
   | 449 | G | 425 |     2.4 | 0.023 | buried | BORDERLINE |
   | 450 | L | 426 |     1.2 | 0.006 | buried | BORDERLINE |
   | 451 | R | 427 |    37.3 | 0.136 | partial |  |
   | 452 | S | 428 |    34.8 | 0.225 | partial | BORDERLINE |
   | 453 | L | 429 |     7.5 | 0.037 | buried | BORDERLINE |
   | 454 | K | 430 |   126.0 | 0.534 | exposed |  |
   | 455 | E | 431 |    77.6 | 0.348 | exposed | ANCHOR |
   | 456 | I | 432 |     0.0 | 0.000 | buried |  |
   | 457 | S | 433 |    27.1 | 0.175 | partial |  |
   | 458 | D | 434 |    47.8 | 0.248 | partial | ANCHOR, BORDERLINE, NORM-SENSITIVE(exposed if empirical) |
   | 459 | G | 435 |     5.9 | 0.057 | partial | BORDERLINE |
   | 460 | D | 436 |    35.9 | 0.186 | partial | ANCHOR |
   | 461 | V | 437 |     1.2 | 0.007 | buried | BORDERLINE |
   | 462 | I | 438 |    29.0 | 0.147 | partial |  |
   | 463 | I | 439 |     0.0 | 0.000 | buried |  |
   | 464 | S | 440 |    29.3 | 0.189 | partial |  |
   | 465 | G | 441 |    28.7 | 0.276 | exposed | BORDERLINE |
   | 466 | N | 442 |     0.0 | 0.000 | buried |  |

   Epitope summary (52 resolved): 13 exposed, 19 partial, 20 buried
   Mean RSA across the block: 0.173

2. VERDICT ON EACH OF THE EIGHT ANCHORS

   This determines how much material we have left to work with.

   | anchor | role | SASA A^2 | RSA (theor.) | RSA (emp.) | VERDICT | usable? |
   |---|---|---|---|---|---|---|
   | D416 | acidic -> binder HIS |    22.5 | 0.117 | 0.120 | PARTIAL | YES |
   | H418 | target HIS -> binder D/E |     7.2 | 0.032 | 0.034 | BURIED *BORDERLINE | no |
   | E421 | acidic -> binder HIS |    46.8 | 0.210 | 0.218 | PARTIAL *BORDERLINE | YES |
   | E424 | acidic -> binder HIS |    95.7 | 0.429 | 0.447 | EXPOSED | YES |
   | H433 | target HIS -> binder D/E |   145.1 | 0.648 | 0.672 | EXPOSED | YES |
   | E455 | acidic -> binder HIS |    77.6 | 0.348 | 0.362 | EXPOSED | YES |
   | D458 | acidic -> binder HIS |    47.8 | 0.248 | 0.256 | PARTIAL *BORDERLINE,NORM-SENSITIVE | YES |
   | D460 | acidic -> binder HIS |    35.9 | 0.186 | 0.192 | PARTIAL | YES |

   Anchors surviving as exposed or partially exposed: 7 of 8
   Survivors: D416, E421, E424, H433, E455, D458, D460
   Lost to burial or absence: H418

3. How solid are these verdicts?

   Two ways a verdict could be an artefact of an arbitrary choice:

   a) Which reference table we normalise by.
      1 anchor(s) change class between the
      theoretical and empirical columns: D458
      Their verdicts are NOT robust and should be treated as
      ambiguous rather than settled.

   b) Where the cutoffs sit (0.25 / 0.05 are conventions, not physics).
      3 anchor(s) sit within 0.05 of a
      cutoff: H418, E421, D458
      Read these as 'somewhere near the boundary', not as the label.

4. CAVEAT NOT IN THE ORIGINAL PLAN: glycan occlusion

   EGFR is a glycoprotein -- sugar chains are attached to it at specific
   points. Those chains are large, and the receptor-only file used above
   contains protein atoms ONLY, because step 02 stripped everything that
   was not a standard amino acid. So the surface measured above is the
   BARE protein. If a sugar chain sits over our epitope, the real
   accessible surface is smaller than computed here.

   Checking the original complex for sugar atoms near the epitope:
   Sugar residues present in 6ARU: 13
   Epitope residues within 5 A of a sugar atom: 3
     T415: 4.44 A from NAG (chain E)
     N444: 1.44 A from NAG (chain E)
     T446: 3.25 A from NAG (chain E)

   No anchor is within 5 A of a sugar. Glycan occlusion is
   not a concern for the anchors specifically.

   LIMIT: a crystal structure resolves only the innermost, most
   ordered sugars. Real glycan chains extend considerably further
   and are mobile. Absence of a modelled sugar is weak evidence of
   absence. Flagged as a residual risk, not resolved.

5. WHAT THIS MEANS FOR THE DESIGN

   7 of 8 anchors are reachable. The pairing rule needs
   three or four pairs stacked, because the pH switch is partial
   rather than binary, so this is enough material to work with --
   PROVIDED they cluster within reach of a single binder, which is
   step 05's question and is not answered here.

   NOTE ON WHAT EXPOSURE DOES AND DOES NOT TELL US. An exposed anchor
   is reachable. It is not therefore a good contact point: it still has
   to point in a compatible direction, sit in a pocket a designed
   backbone can present a partner to, and cluster with the others.
   Exposure is a filter that removes impossible options, not evidence
   that the remaining ones work.

Wrote data/derived/03-epitope-rsa.csv
Wrote data/derived/03-anchor-verdicts.csv

========================================================================
RESULT: 7 of 8 anchors usable (D416, E421, E424, H433, E455, D458, D460).
Clustering is step 05's question. This step only removed the buried.
========================================================================
```

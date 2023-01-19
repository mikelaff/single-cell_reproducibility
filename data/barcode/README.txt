WTK1-WTK6 all use version 1 chemistry from parse biosciences.

Barcode files pulled from the parse biosciences pipeline scripts.

From their pipeline:
KIT_CHEM_DEFS = """
    kit       chem  nwells  bc1     bc2      bc3    ktype
    custom    NA    1       NA       NA       NA    special
    WT_mini   v1    12      n24_v4   v1       v1    normal
    WT_mini   v2    12      n24_v4   v1       v1    normal
    WT        v1    48      v2       v1       v1    normal
    WT        v2    48      n96_v4   v1       v1    normal
    WT_mega   v1    96      n192_v4  v1       v1    normal
    WT_mega   v2    96      n192_v4  v1       v1    normal
"""

We are using WT v1 48, therefore the order of barcodes comes from the files:
bc_data_v2.csv - bc_data_v1.csv - bc_data_v1.csv

Inside bc_data_v2.csv are 48 barcodes (A1-D12) for oligo dT and 48 barcodes
(A1-D12) for random heximer. To use the splitp algorithm, the file 
IVIV_oligo_hex_bc_mapping.txt was generated to map between oligo dT and
random heximer.

FILES:
bc_data_v1.csv
round 2 and round 3 barcode list for 96 wells A1-H12.

bc_data_v2.csv
round 1 barcode list for 48 wells A1-D12 oligo dT and round 1
barcode list for 48 wells A1-D12 random heximer.

IVIV_oligo_hex_bc_mapping.txt
mapping between oligo dT and random heximer barcodes.

fullList_barcode_combination.txt
all barcode combinations in the order: bc_data_v2.csv - bc_data_v1.csv - bc_data_v1.csv
using only oligo dT barcodes. this file is used by alevin-fry to create a permit list, 
to define the list of cell barcodes to quantify.

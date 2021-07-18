"""
Test dna synthesis functions.
"""

import pandas as pd
import tempfile
from mutagenesis_visualization.main. import (
    generate_primers, _primerdesign, _create_primers_list, create_variants
)

def test_generate_primers():
    dna = 'GGCAATGCGcccccaATGaaaaaaTAAaaACGGGGTTTTaaa'
    start = ('GGCAATGCGccccca')
    end = ('aaACGGGGTTTTaaa')

    #create a temporary directory using the context manager
    with tempfile.TemporaryDirectory() as tmpdirname:
        df_prim = generate_primers(
            dna,
            start,
            end,
            output_file=tmpdirname+'/primers.xlsx',
            codon='NNS',
            length_primer=18,
            tm=60,
            return_df=True
        )
    assert len(df_prim) == 9, 'error the primer length does not match'
    assert df_prim.iloc[
        1, 1
    ] == 'NNSGCGCCCCCAATGAAAAAATAAAAACGGG', 'error the first primer is not mutated at Met'
    assert df_prim.iloc[
        8, 1
    ] == 'CGCCCCCAATGAAAAAANNSAAACGGGGTTTTAA', 'error the last primer is not mutated at the stop codon'

def test_primerdesign():
    dna = 'GGCAATGCGcccccaATGaaaaaaTAAaaACGGGGTTTTaaa'  # 42 codons
    primers0 = _primerdesign(
        dna, codon='NNS', codonposition=0, length_primer=None, tm=60
    )
    primers1 = _primerdesign(
        dna, codon='NNS', codonposition=1, length_primer=None, tm=60
    )
    primers2 = _primerdesign(
        dna, codon='NNS', codonposition=2, length_primer=None, tm=60
    )
    primers3 = _primerdesign(
        dna, codon='NNS', codonposition=3, length_primer=None, tm=60
    )
    primers21 = _primerdesign(
        dna, codon='NNS', codonposition=21, length_primer=18, tm=None
    )
    assert type(primers3) == tuple, 'error the output is not a tuple'

    assert primers0 == (
        'NNSAATGCGcccccaATGaaaaaaTAAaaACGG', 'CCGTttTTAttttttCATtgggggCGCATTSNN'
    ), 'error the output is not the forward and reverse primers at codon position 0'
    assert primers1 == (
        'NNSATGCGcccccaATGaaaaaaTAAaaACGGG', 'CCCGTttTTAttttttCATtgggggCGCATSNN'
    ), 'error the output is not the forward and reverse primers at codon position 1'
    assert primers2 == (
        'NNSTGCGcccccaATGaaaaaaTAAaaACGGG', 'CCCGTttTTAttttttCATtgggggCGCASNN'
    ), 'error the output is not the forward and reverse primers at codon position 2'
    assert primers3 == (
        'NNSGCGcccccaATGaaaaaaTAAaaACGGG', 'CCCGTttTTAttttttCATtgggggCGCSNN'
    ), 'error the output is not the forward and reverse primers at codon position 3'

    assert primers21 == (
        'AATGCGcccccaATGaaaNNSTAAaaACGGGGTTTTaaa',
        'tttAAAACCCCGTttTTASNNtttCATtgggggCGCATT'
    ), 'error the codon is not replaced after position 21'


# In[4]:


def test_create_primers_list():
    dna = 'GGCAATGCGcccccaATGaaaaaaTAAaaACGGGGTTTTaaa'

    primerslist = _create_primers_list(
        dna, start_codon=3, end_codon=6, codon='NNS', length_primer=15, tm=60
    )
    assert type(primerslist[0]) == list, 'error the type is not a list'


def test_create_variants():
    dna = 'AAGAAGAAG'
    codon_list1 = ['GAA', 'TAA', 'AAA']
    variants1 = create_variants(
        dna, codon_list1, output_file=None, return_df=True
    )
    assert variants1.shape == (
        10, 1
    ), 'error there are an incorrect number of variants'
    assert variants1.iloc[
        0, 0] == 'AAGAAGAAG', 'error the first output is not wild type'
    assert variants1.iloc[
        1, 0
    ] == 'GAAAAGAAG', 'error the second output is not replaced with the first codon'
    assert variants1.iloc[
        2, 0
    ] == 'TAAAAGAAG', 'error the third output is not replaced with the second codon'

    #create a temporary directory using the context manager
    with tempfile.TemporaryDirectory() as tmpdirname:
        codon_list2 = ['CAA', 'CCC', 'TTT', 'AAA']
        variants2 = create_variants(
            dna, codon_list2, output_file=tmpdirname+'/test.fasta', return_df=True
        )
        variants2 = create_variants(
            dna, codon_list2, output_file=tmpdirname+'/test.xlsx', return_df=True
        )
    assert variants2.iloc[
        1, 0
    ] == 'CAAAAGAAG', 'error the second output is not replaced with the first codon'
    assert variants2.shape == (
        13, 1
    ), 'error there are an incorrect number of variants'

## Sequence Objects
```python
from Bio.Seq import Seq

```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index,letter))
    
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
print(len(my_seq))
```

    5



```python
# We can also print the lenght of each sequence 
print(len(my_seq))
```

    5



```python
print(my_seq[0])
```

    G



```python
Seq("AAAA"). count("AA")
```




    2




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
len(my_seq)
```




    32




```python
my_seq.count("G")
```




    9




```python
100 * (my_seq.count("G")+ my_seq.count("C")) / len(my_seq)
```




    46.875




```python
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
gc_fraction(my_seq)
```




    0.46875




```python
my_seq[4:12]
```




    Seq('GATGGGCC')




```python
my_seq[0::3]
```




    Seq('GCTGTAGTAAG')




```python
my_seq[1::3]
```




    Seq('AGGCATGCATC')




```python
my_seq[2:3]
```




    Seq('T')




```python
my_seq[::-1]
```




    Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG')




```python
str(my_seq)
```




    'GATCGATGGGCCTATATAGGATCGAAAATCGC'




```python
fasta_format_string = "Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    Name
    GATCGATGGGCCTATATAGGATCGAAAATCGC
    



```python
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
contigs = [Seq("ATG"), Seq("ATCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N" *10)
```


```python
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCGNNNNNNNNNNTTGCA')


```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```

    Seq('acgtACGT')




```python
"GTAC" in dna_seq
```


    False


```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq.complement()
```




    Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')




```python
my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')




```python
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
coding_dna = Seq("ATGGC")
```


```python
coding_dna
```




    Seq('ATGGC')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCGATAG')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna

```




    Seq('CTATCGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCGAUAG')




```python
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCGAUAG')




```python
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCGATAG')




```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCGAUAG')




```python
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAD')




```python
coding_dna.translate(table= "Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAD')




```python
coding_dna.translate(table = 2)
```




    Seq('MAIVMGRWKGAD')




```python
coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
coding_dna.translate(table =2, to_stop=True)
```




    Seq('MAIVMGRWKGAD')




```python
gene = Seq("GTGAAAAAGATCAATCTATCGTACTGCATTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCC")
```


```python
gene.translate(table = "Bacterial")
```


    Seq('VKKINLSYCISLVLVAPMAAQAAEITLV')


```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKINLSYCISLVLVAPMAAQAAEITLV')




```python
gene.translate(table = "Bacterial", cds = True)
```
    <ipython-input-83-9b7bfe3f3060> in <module>
    ----> 1 gene.translate(table = "Bacterial", cds = True)
    

    ~/anaconda3/lib/python3.7/site-packages/Bio/Seq.py in translate(self, table, stop_symbol, to_stop, cds, gap)
       1446 
       1447         return self.__class__(
    -> 1448             _translate_str(str(self), table, stop_symbol, to_stop, cds, gap=gap)
       1449         )
       1450 


    ~/anaconda3/lib/python3.7/site-packages/Bio/Seq.py in _translate_str(sequence, table, stop_symbol, to_stop, cds, pos_stop, gap)
       2791         if n % 3 != 0:
       2792             raise CodonTable.TranslationError(
    -> 2793                 f"Sequence length {n} is not a multiple of three"
       2794             )
       2795         if str(sequence[-3:]).upper() not in stop_codons:




```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
seq1 = "ACGT"
```


```python
unknown_seq = Seq(None,10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAAGGATAGT"}, length = 159345973)
```


```python
seq[1000:1020]
```




    Seq(None, length=20)




```python
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length =10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
from Bio.Seq import MutableSeq
```




```python
mutable_seq = MutableSeq(my_seq)
```
\

![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/0236dd0f-e02e-4d3c-9773-ea8cb4161d37)



## Sequence Annotation

```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
simple_seq_r.id = "AC12345"
```


```python
simple_seq_r.description = "Made up sequence for the VDB Computational Biology Class"
```


```python
print(simple_seq_r.description)
```

    Made up sequence for the VDB Computational Biology Class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
simple_seq_r.annotations["evidence"] = "None. This is just an example"
```


```python
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example



```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
simple_seq_r.letter_annotations
```




    {'phred_quality': [40, 40, 38, 30]}




```python
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```
record

```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.dbxrefs
```




    []




```python
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'NC_005816.1'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations
```




    {}




```python
len(record.annotations)
```




    13




```python
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
record.dbxrefs
```




    ['Project:58037']




```python
len(record.features)
```




    41




```python
from Bio import SeqFeature
```


```python
start_pos = SeqFeature.AfterPosition(5)
```


```python
end_pos = SeqFeature.BetweenPosition(9, left =8, right =9)
```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>5:(8^9)]



```python
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
int(my_location.end)
```




    9




```python
int(my_location.start)
```




    5




```python
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)




```python
shifted = record [2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
len(shifted)
```




    9609




```python
len(shifted.features)
```




    40




```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
shifted.dbxrefs
```




    []




```python
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
rc = record.reverse_complement(id="Testing")
```


```python
rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
print("%s %i %i %i %i" %(rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0



![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/69a77b68-198c-4693-bc5c-c2afc99dbd03)


## Sequence Input/Output

```python
from Bio import SeqIO
```


```python
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    gi|2765657|emb|Z78532.1|CCZ78532
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    gi|2765656|emb|Z78531.1|CFZ78531
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    gi|2765655|emb|Z78530.1|CMZ78530
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    gi|2765654|emb|Z78529.1|CLZ78529
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    gi|2765652|emb|Z78527.1|CYZ78527
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    gi|2765651|emb|Z78526.1|CGZ78526
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    gi|2765650|emb|Z78525.1|CAZ78525
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    gi|2765649|emb|Z78524.1|CFZ78524
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    gi|2765648|emb|Z78523.1|CHZ78523
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    gi|2765647|emb|Z78522.1|CMZ78522
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    gi|2765646|emb|Z78521.1|CCZ78521
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    gi|2765645|emb|Z78520.1|CSZ78520
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    gi|2765644|emb|Z78519.1|CPZ78519
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    gi|2765643|emb|Z78518.1|CRZ78518
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    gi|2765642|emb|Z78517.1|CFZ78517
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    gi|2765641|emb|Z78516.1|CPZ78516
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    gi|2765640|emb|Z78515.1|MXZ78515
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    gi|2765639|emb|Z78514.1|PSZ78514
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    gi|2765638|emb|Z78513.1|PBZ78513
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    gi|2765637|emb|Z78512.1|PWZ78512
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    gi|2765636|emb|Z78511.1|PEZ78511
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    gi|2765635|emb|Z78510.1|PCZ78510
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    gi|2765634|emb|Z78509.1|PPZ78509
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    gi|2765633|emb|Z78508.1|PLZ78508
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    gi|2765632|emb|Z78507.1|PLZ78507
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    gi|2765631|emb|Z78506.1|PLZ78506
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    gi|2765630|emb|Z78505.1|PSZ78505
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    gi|2765629|emb|Z78504.1|PKZ78504
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    gi|2765628|emb|Z78503.1|PCZ78503
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    gi|2765627|emb|Z78502.1|PBZ78502
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    gi|2765626|emb|Z78501.1|PCZ78501
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    gi|2765625|emb|Z78500.1|PWZ78500
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    gi|2765624|emb|Z78499.1|PMZ78499
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    gi|2765623|emb|Z78498.1|PMZ78498
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    gi|2765622|emb|Z78497.1|PDZ78497
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    gi|2765621|emb|Z78496.1|PAZ78496
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    gi|2765620|emb|Z78495.1|PEZ78495
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    gi|2765619|emb|Z78494.1|PNZ78494
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    gi|2765618|emb|Z78493.1|PGZ78493
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    gi|2765617|emb|Z78492.1|PBZ78492
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    gi|2765616|emb|Z78491.1|PCZ78491
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    gi|2765615|emb|Z78490.1|PFZ78490
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    gi|2765614|emb|Z78489.1|PDZ78489
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    gi|2765613|emb|Z78488.1|PTZ78488
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    gi|2765612|emb|Z78487.1|PHZ78487
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765611|emb|Z78486.1|PBZ78486
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    gi|2765610|emb|Z78485.1|PHZ78485
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    gi|2765609|emb|Z78484.1|PCZ78484
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    gi|2765608|emb|Z78483.1|PVZ78483
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    gi|2765607|emb|Z78482.1|PEZ78482
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    gi|2765606|emb|Z78481.1|PIZ78481
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    gi|2765605|emb|Z78480.1|PGZ78480
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    gi|2765604|emb|Z78479.1|PPZ78479
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    gi|2765603|emb|Z78478.1|PVZ78478
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    gi|2765602|emb|Z78477.1|PVZ78477
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    gi|2765601|emb|Z78476.1|PGZ78476
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    gi|2765600|emb|Z78475.1|PSZ78475
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    gi|2765599|emb|Z78474.1|PKZ78474
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    gi|2765598|emb|Z78473.1|PSZ78473
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    gi|2765597|emb|Z78472.1|PLZ78472
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    gi|2765596|emb|Z78471.1|PDZ78471
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765595|emb|Z78470.1|PPZ78470
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    gi|2765594|emb|Z78469.1|PHZ78469
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    gi|2765593|emb|Z78468.1|PAZ78468
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    gi|2765592|emb|Z78467.1|PSZ78467
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    gi|2765591|emb|Z78466.1|PPZ78466
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    gi|2765590|emb|Z78465.1|PRZ78465
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    gi|2765589|emb|Z78464.1|PGZ78464
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    gi|2765588|emb|Z78463.1|PGZ78463
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    gi|2765587|emb|Z78462.1|PSZ78462
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    gi|2765586|emb|Z78461.1|PWZ78461
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765585|emb|Z78460.1|PCZ78460
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    gi|2765584|emb|Z78459.1|PDZ78459
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    gi|2765583|emb|Z78458.1|PHZ78458
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    gi|2765582|emb|Z78457.1|PCZ78457
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    gi|2765581|emb|Z78456.1|PTZ78456
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765580|emb|Z78455.1|PJZ78455
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    gi|2765579|emb|Z78454.1|PFZ78454
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    gi|2765578|emb|Z78453.1|PSZ78453
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    gi|2765577|emb|Z78452.1|PBZ78452
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    gi|2765576|emb|Z78451.1|PHZ78451
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    gi|2765575|emb|Z78450.1|PPZ78450
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    gi|2765574|emb|Z78449.1|PMZ78449
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    gi|2765573|emb|Z78448.1|PAZ78448
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    gi|2765572|emb|Z78447.1|PVZ78447
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    gi|2765571|emb|Z78446.1|PAZ78446
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    gi|2765570|emb|Z78445.1|PUZ78445
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    gi|2765569|emb|Z78444.1|PAZ78444
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    gi|2765568|emb|Z78443.1|PLZ78443
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    gi|2765567|emb|Z78442.1|PBZ78442
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    gi|2765566|emb|Z78441.1|PSZ78441
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    gi|2765565|emb|Z78440.1|PPZ78440
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    Z78532.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    Z78531.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    Z78530.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    Z78529.1
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    Z78527.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    Z78525.1
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    Z78524.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    Z78523.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    Z78522.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    Z78521.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    Z78520.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    Z78519.1
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    Z78518.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    Z78517.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    Z78516.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    Z78515.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    Z78514.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    Z78513.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    Z78512.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    Z78511.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    Z78510.1
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    Z78509.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    Z78508.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    Z78507.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    Z78506.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    Z78505.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    Z78504.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    Z78503.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    Z78502.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    Z78501.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    Z78500.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    Z78499.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    Z78498.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    Z78497.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    Z78496.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    Z78495.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    Z78494.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    Z78493.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    Z78492.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    Z78491.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    Z78490.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    Z78489.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    Z78488.1
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    Z78487.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78486.1
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    Z78485.1
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    Z78484.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    Z78483.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    Z78482.1
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    Z78481.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    Z78480.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    Z78479.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    Z78478.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    Z78477.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    Z78476.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    Z78475.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    Z78474.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    Z78473.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    Z78472.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    Z78471.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78470.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    Z78469.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    Z78468.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    Z78467.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    Z78466.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    Z78465.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    Z78464.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    Z78463.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    Z78462.1
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    Z78461.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78460.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    Z78459.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    Z78458.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    Z78457.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    Z78456.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78455.1
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    Z78454.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    Z78453.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    Z78452.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    Z78451.1
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    Z78450.1
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    Z78449.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    Z78448.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    Z78447.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    Z78446.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    Z78445.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    Z78444.1
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    Z78443.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    Z78442.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    Z78440.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
identifiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")]
```


```python
identifiers
```




    ['Z78533.1',
     'Z78532.1',
     'Z78531.1',
     'Z78530.1',
     'Z78529.1',
     'Z78527.1',
     'Z78526.1',
     'Z78525.1',
     'Z78524.1',
     'Z78523.1',
     'Z78522.1',
     'Z78521.1',
     'Z78520.1',
     'Z78519.1',
     'Z78518.1',
     'Z78517.1',
     'Z78516.1',
     'Z78515.1',
     'Z78514.1',
     'Z78513.1',
     'Z78512.1',
     'Z78511.1',
     'Z78510.1',
     'Z78509.1',
     'Z78508.1',
     'Z78507.1',
     'Z78506.1',
     'Z78505.1',
     'Z78504.1',
     'Z78503.1',
     'Z78502.1',
     'Z78501.1',
     'Z78500.1',
     'Z78499.1',
     'Z78498.1',
     'Z78497.1',
     'Z78496.1',
     'Z78495.1',
     'Z78494.1',
     'Z78493.1',
     'Z78492.1',
     'Z78491.1',
     'Z78490.1',
     'Z78489.1',
     'Z78488.1',
     'Z78487.1',
     'Z78486.1',
     'Z78485.1',
     'Z78484.1',
     'Z78483.1',
     'Z78482.1',
     'Z78481.1',
     'Z78480.1',
     'Z78479.1',
     'Z78478.1',
     'Z78477.1',
     'Z78476.1',
     'Z78475.1',
     'Z78474.1',
     'Z78473.1',
     'Z78472.1',
     'Z78471.1',
     'Z78470.1',
     'Z78469.1',
     'Z78468.1',
     'Z78467.1',
     'Z78466.1',
     'Z78465.1',
     'Z78464.1',
     'Z78463.1',
     'Z78462.1',
     'Z78461.1',
     'Z78460.1',
     'Z78459.1',
     'Z78458.1',
     'Z78457.1',
     'Z78456.1',
     'Z78455.1',
     'Z78454.1',
     'Z78453.1',
     'Z78452.1',
     'Z78451.1',
     'Z78450.1',
     'Z78449.1',
     'Z78448.1',
     'Z78447.1',
     'Z78446.1',
     'Z78445.1',
     'Z78444.1',
     'Z78443.1',
     'Z78442.1',
     'Z78441.1',
     'Z78440.1',
     'Z78439.1']




```python
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record.id)
```

    gi|2765658|emb|Z78533.1|CIZ78533



```python
print(first_record.description)
```

    gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
second_record = next(record_iterator)
```


```python
print(second_record.id)
```

    gi|2765657|emb|Z78532.1|CCZ78532



```python
print(second_record.description)
```

    gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
records = list(SeqIO.parse("ls_orchid.gbk.txt", "genbank"))
```


```python
print("Found %i records" % len(records))
```

    Found 94 records



```python
print("The last record")
last_record = records [-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
print("The first record")
first_record = records [0]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    The first record
    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740



```python
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")
]
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    all_species.append(seq_record.description.split()[1])
```


```python
print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
first_record.id
```




    'gi|2765658|emb|Z78533.1|CIZ78533'




```python
first_record
```




    SeqRecord(seq=Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), id='gi|2765658|emb|Z78533.1|CIZ78533', name='gi|2765658|emb|Z78533.1|CIZ78533', description='gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', dbxrefs=[])




```python
first_record.description = first_record.id + " " + "Mutations induced randomly"
```


```python
print(first_record.format("fasta"[:200]))
```

    >gi|2765658|emb|Z78533.1|CIZ78533 Mutations induced randomly
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCC
    CGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCC
    CAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAA
    CGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTG
    AATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCA
    GGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCG
    GCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCG
    GCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTG
    GCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCC
    TTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGG
    GGCACCCGCTGAGTTTACGC
    



```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
rec1 = SeqRecord (
    Seq("MMYQ0GCFAGGTVLRLAKDLAENNRGARVLWCSEITAVTFRGPSETHLDSMVGQALFGD"
        "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
        "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
        "SSAC" ,
       ),
    
    
 id = "gil14150838|gb\AAK54658.1|AF376133 1",
description = "chalcone synthase [Cucumis sativus]",
)
```


```python
rec2 = SeqRecord(
    Seq("YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEMAPSLDARQ"
 "DMVWEIPKLGKEAAVKAIKEWGQ",
       ),
    id = "gi | 13919613 |gb AAK33142.1|" ,
    
description = "chalcone synthase [Fragaria vesca subsp. bracteata]"
)




```


```python
rec3 = SeqRecord(
Seq("MVTVEEFRRAOCAEGPATVMAIGTATPSNCVDOSTYPDYYFRITNSEHKVELKEKFKRMC"
    "EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARODIVVVEVPKLGKEAA0KAIKEWGOP"
    "KSKITHLVFCTTSGVDMPGCDYOLTKLLGLRPSVKRFMMYO0GCFAGGTVLRMAKDLAEN"
    "NKGARVLVVCSEITAVTFRGPNDTHLDSLVGOALFGDGAAAVIIGSDPIPEVERPLFELV"
    "SAAOTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAF0PLGISDWNSLFW"
    "IAHPGGPAILDQVELKLGLKQEKLKATKVLSNYGNMSSACVLFILDEMRKASAKEGLGT"
    "TGEGLEWGVLFGFGPGLTVETVVLHSVAT" ,
   ),

    id="gi|13925890|gb|AAK49457.1|", 
    description="chalcone synthase [Nicotiana tabacum]",
)
```


```python
my_records = [rec1, rec2, rec3]
```


```python
from Bio import Seq
SeqIO.write(my_records, "my example.faa", "fasta")
```




    3




```python
records= SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
count= SeqIO.write(records, "my_example.fasta", "fasta")
print("Converted %i records" % count)
```

    Converted 94 records



```python
for record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(record.id)
    print(record.seq.reverse_complement())
```

    Z78533.1
    GCGTAAACTCAGCGGGTGCCCCCGCCTGACCTGGGGTCACATCCGAATGGCGGTCAACCGCCCTCCATTGGGGTTCGGAAGGGTTCTCCTGCCTGTCCGACAAGCACGACAACATGGGGGATTCGCACGGCAGCTGCTGCCAGCATCCGTCCACCTCTTGCCGGGTTCCGGGCCATCAAAACACCAGCTCTTGGACCCGCCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACCACGCCGGCCTGGCTGTATGCCGGGCAAGCATTGGCAGGAGAGAGACGAAGCGCGACGCCCAAGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCAAGATTCCCGTTTGCGAGAGTCATCAAAATTCATTGGGCACGCGACAGCACGCCGCCGCTCCGGGTTTTGGGGAAGACAATGCCATTCGCCGGTGATGCTTTCATATGGCTTGGCGCCCAAACTGCGCCGGGCTAGAGGTTCAAACCCGCCATGGACGCTCCCGAGGCGGCCCAACAACAAATCAGGGTCACCACGGGAGCAATGCCCCCGGTGAGCTGAGTACACCGGTCCTCCGGATTCACTCGATCGTTTATTCCACGGTCTCATCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78532.1
    GCCTCAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTGAATGGAAATCAACTGCCCAATGGTTATTTTAGCTCCATTGGGGTTCAATTAGGTTCTTGTGTAGGTTCGAAAAAATACAACAACATGGGGGATTCAAATAGCAGCCTTATGACTGTTAGCATTCTCCACCTCGTGCCACATTCCTACCCATCAAAGCAACAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCATTCACATCCGTATAATGCCAGCTTAGCGATATGCCAAGCAAGCATTGGTAGGAGAGACGCAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTAGCTGCGTTCTTCATCGATGTTAGAGCCAAGATATCCGTTGCGAGAGTCATAAAAATTCACTGGCATGCTCAGTAGCATACTGCCCCTCTGATTTTTTCTGACAATAATGTCATTCATCAGTGATGCTTTATATGACTTGGCGCAAACTGCACCGTACTAAAGTTCAAACCTGCCATGAAAGCTCTTGAGGAGGCCCAACAACAAAGCAGGGTCACGACAAAAGCAGTGCCACGACGAGCTGAGTTACCACAGGTCCTCCAGATTCACTCGATCATATATTCTGTTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78531.1
    TTACTCACGGGTGGCCCGCCTGACTGGGGTCGCATCTGAATGGAAATCACCGCCCGAGGGCGGTTTTGCGTCCACTGGGGGTTCAAACAGGTTCTTCTGTAGGCCTGACGAGTACGACAACACGGGGGATTCGAACGGCAGCCTTGCGGCTGCCAGCATCCTCCACCTACTGCCAAGTTCCGGGCCATCAGAACACCGATGCTTAGACCCGCCGCACCTAGGCGCAAGGGGCCAATCTTTCACGTCCTCACAATGACAGACTAGCCGCATGCCAATCAAGCATTATCAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAGATATCCCGTTGCCGAGAGTCGTCATAAATTCATTGGCACGCGCAACAGACGCCGCCCCTCCGAGTTTTTCTTGACAAAAATGCCATCCATCGGTGACGCTCCATATGACTTGGCGCAAACTGCGCCGCGCTAGACGTTCAAACCCGCCATGAAAGTTCCCGAGGCGGCCCGGTCGCAAACCGGGTTCACCACGAGAGCAAAGCCACGGTGAGCCGTGTAACCACGGGTCCTCCGGATTCACTCGATCGTATGTTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78530.1
    ATGGGGTGGCCCCCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAACATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGCGTTCTTCATCGATGCAAGAGTCAAGATATCCGTCTGCCGAGAGTCATCAAAATTCATTGGAATGAACATTAACACACCACACCTCCGACGTTGTGTGGGGCAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGGACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGATATCCCTAGCGAGCCAAATTACCACAAGTCCTCCAGATTCACTCAATCGTTTATTATGTTGTTTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78529.1
    TTTTAGCTCAGCTGGTGGCCCGCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAATATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGGGATTCTTCATCGATGCAAGAGCCAAGATATCCCTCGGCGAGAGTCATCAACAATTCATTGGCATGACCAATAACACACCACACCTCCGACTTTTTTGACAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGCACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGAAATCCCTAGCGAGCCAAATAACCACAAGTCCTCCAGATTCACTCAATCGTATATTCTGCTGTCTCAACAATGTCCTTCGGCAGCTCGCCGT
    Z78527.1
    GGGTGGCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACAATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACCCAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAAGAGCCAGATATCCGTTGGCCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAACACACTGCCACTCCAACTTTTGGTGACAGTAGTGCCATTTGCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGCCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGCGAGCCGAGTAACCACAGGTCATCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78526.1
    ACATAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGGCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAGCACACTGCCACTCCAACTTTTGGTGACAATAATGCCATTGTTCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGGCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGGGAGCCGAGTAACCACAGGTCCTCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78525.1
    TGCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAGTCAACCATCCAAGGGTCATTTTGCCTCCACTGGGGTTCAAACAAGTTATTATATAGGTCCGATAAATACGATAACATGGGGGATTTAAATGACAACCTTGTGGCAGCCAATGTCCTCCACCTCCTGCCAAGTTTCGGGCCATCAAAACATTCCTTAGACCCAACGCGCCAAGGCACAAGGGGCCAATCTTTCGCATCCGCACAATGCAAGCCTAGATATATGGACAGGCATTGGCAAGAGAGACACAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCCGATGCAAGAGCCAAGATATCCGCTTGTCGAGAGTCATCAAAATTATATTCGGCATGCACAGGAGGACACCGCCCCTCCAACTTTTTTGACAATAATGTCATTTGTCGGTGATGCTTCATATGACTTGGCACAAACTGTGCCGTGCTAGAGTTTGAAACCTGCCATGAAAGCTCCCGAGGAGGCCCAACGACAAATTTGGGTCACCACAAAAGCAAAGCCCTCGGCAAGCCGAATAACCACAGGTCCTCCGGATTCACTCGATGTATATTCTGCTATCTCAACA
    Z78524.1
    GCTTAATCTCAGGGTGGCCACCTGACCTGGGGTTGCATCTGAATGAAAATCAACCGCCCAAGGGTCATTTTTCCTTCATTGGGGTTCAAACAAGTTCTTTTATAGGTTCAACAAATACGACAACATGGGGAATTCAAATGGCAGCCTTGCAGCTACCAACATTCTCCACCTCCTGCCGAGTTCCGAGCTATCAAAATACTAATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATATTTCATCCGCACAATGCTAACCTAGCTATATGCTGGGCAAGCATTGGNAGGAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGGATGTAAAGAGCCAAGTTCCGTTGCGAGAGTCATCAAAAAATTCATTAGCATGCACAACAGCATGCCGCCCCTCCGACTTTTTTAACAACAATGCCATTCATCAGGGATGCTTATATGACTTGGCACAAACTGCGCCGTGCTAAAGGTTTGAACCCGCCATGAAAGCTCCTAAGGAGGCCCAATTAAGGTCACCACAAAAGCAAAGCACTGACGAGCCAAGTAACCACATGTCCTCCATATTCACTCAATCATATATTCTACTATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78523.1
    CTTAAACTCAGCGGTGGCCCCGCCTGACCTGGGTCGCAATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGGTTCAAACGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCGCATCCACGCAATGCAAACCAATCTATATGATAGGAAAGCATTGATAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACATATGTCGAGAGTCAATAAAAAATTCATTATGCATGCATGCAAGAGCACGCTTCCACTCTAGCTTTTGGGACAATAATGTCATTCCTCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCACACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAGTTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTGCCGCAGGTTCACCTACGGAAACCTGGTTACG
    Z78522.1
    CTCAGCGGGTGGCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCCCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACGCAATGCAAACCTGTCTATATGACGTATAACCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTTCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCGATCACACCACGATATGTCGAGAGTCATAAAAAATTCATTTGCATGTATGCAATAGCACGCTTCCACTCCAACTTTTTGGACAATAATGTCATTCATCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGGGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78521.1
    GGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGAGTTGTTTGCCTCCATTGGGATTCAAACATGTTTTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGTACCTAGGCACAAGGGGCCATTCTTTCACATCCACGCAATGCAAACCTATCTATATGACATGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTTGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAATTCATTTGCATGCATGCAATAGCACGCTTCCACTCCAACTTTTTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGGCAAATTAGGGTCACCGCCAAAGCCAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78520.1
    AAACTCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCCAGGTCCGACAAACACGACAATATGGGGGATTCACATGACAACCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCGATATGACAGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCACATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTTCGCTGCGTTCTTCATCCGATGCAAGAGCCAAGATATCCCTTGTCGAGAGTCATAAAATATTCTATTGCATGCATGCAATAGCACGCTTCTACTCCAACTTTTTGGACAATAATGTCATTCATCGGTGATGATTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78519.1
    TAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCAACAAACGCGACAATATGGGGGATTCAAATGATAGCCTTATATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACATTAAGTTCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCTATATGATGGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAAATTCATTTGCATGCATGCAGTAGCACGCTTCCACTCCAACTTTTGGACAATAATGTCATTCATCGGCAATGCTTCATATGACTTGGTGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCAGATTCACTCGATCATAT
    Z78518.1
    GGAAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCAAGGGTTGTTTTGCCTCCATAGGGGTTCAAATAGGTTCTTCTGTGGGTCCGGCAAATACGACAACATGTGGGATTTGAACGACAGCCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTTGGGCCATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGAGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGTATGCGTAACACACACCGTCCATCCGAGTTTGGGGGACAATAATGTAATTCGTCGGCGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGCTCAAACCCGCCATAAAAGCTCTCGAGGATGCCCAACTACAAATCAGGGTCACCACAAAAGTTAAGCCTTCGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCGAATATTCTACTATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78517.1
    GCTTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAATGACAGTCTTGCGACTATCAACATGATCCACCTCCTGCCGGGGTTCGGGCCACCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGTCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGGAAGAGCAAGATATCCGTTGGCCGAGAGTCATCAAAAATTTATTGGCATGCACAACAGCACACCGCCCATCCGACTTTTGTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCAAACTGCACCGTGATAGAGGTTCAAACCCGCCATAAAATCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCTGGATTCACTCGATCGTATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78516.1
    TTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCATATCTGAATGGCAATCAATCACTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAACGATAGTCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTCGGGACATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTAAGATATCCGTTGACGAGAGTCATCAAAAAATTATTGGTATGCACAATAGCACACCGCCCATCCGACTTTGGTGACAATAATGTCATTCGTCGGTGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGTTCAAACCCGCCATAAAAGCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCATATATTATACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78515.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCTTCCAAATGGCAATCAACCGCCCATGGGTTGGGTTGATGTGGCACCAATGGGGTTCCGAAGGGTCCTACGAATTATGATGGCCCATCAAATCGGACAACATACGGCATTCGCACAACAGCTTTGGCTCAGCCATTGCCCATCCACCTCTTGTCGAATTTTGGACCATCAAAAGCCCGAGGTCTTAGACCCATTGAACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGACCGAACGGTATTCTTGGACAACCATTGGTGGGAGAGACACAGTACACAACGCCCAGGCAGGCGTGCCCTTAGCCTGACGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTAAAAATGTAACTTGGGGGGCACGAATCAAGTGCCACCCCCACCGGCTACAAGGGAACACCCCTATGCCGATTCTTGTCTCAAAAGACTTGGCGCGAAGCTGCGCCGGGCTAGAAGTTTGATCGCCTTCAAGAACACTCCCAAGGAGGCCCGATGGCAGATCAAAGGCCACCACAGTAGCGAAAGCCCCTGGGAGATCGAAGTACACCGGTCCTCCAGATTAACTCGATCGTTTATTCTACGGTCTCAGCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78514.1
    TAGCGGGTAGTCTCACCTGACCTGGGGTTGCATCCTAATGGCCGTCAACCGCCCATGGGTTGGGTCGACGTGCCTCCGACGGGGTTCGAAAGGGATCATTCCGGCCCATCAAATACGACAACATAGGGCATTCGCACAACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGACGAGAGTTGTCAGACGAATCTTAAGGGGGCACAGCAGCACGCCGGCCCTCCCCGGTCCATCCATGGAGGACGAAGACCCTTGTCGGTTCTATCTCATTTGACTTGGAGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATAGAGAACTCTCCCAAGGAGGTTCGATGGGAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTATCTCGATCGTATATTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78513.1
    CTCAGCGGTAGCTCACTGACTGGGTTGCATCCTAATGGCGTCACCGCCCATGGGTTGGGTCGACGTGCCTCCGAAGGGGTTCGAAAGGGATCGTTCTGGCACATCAAATACGACAACATAGGGCATTCGCACAAAAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACAGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCGAAAATTCTTTCGGGCACAGCAGCATGCCGGCCTCCCCGGCTCCATGGAGGACGATGCCCTTGCCGGTTCTATCTCATTTGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATGGAGAAATCTCCCAATGAGGCTCGATGGCAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTAACTCGATCGTGTATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78512.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGGATCCAAATGACCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAGCTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCGAGAGTTGTCAAAAAATAATTTCATTGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGGGGACGATGCCCTCCGCCGGTTCCATCTCATTTGATTTGGGGCGAAACTGCGCCGGGCCAAGGGTTTCATTTGCCATCAAGAATTCCTCCCAGAGGGTTCGACGGAAATTCACGGTCACCACAGTAGCCAAGGCCCCTGGGAGACCAAACTACACCGGCCCTCCGGATTAATTCGATCGTTTTTTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78511.1
    TCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGCTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTTCATCGGATGCAAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAAAAATTCATGGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAATTCTCCCAAGGAGGTTCGACGGAAATTCACGGCCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78510.1
    TAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAAAACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAATCTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCCATCCGATGCAAAGAGGCGAGATATCCGTTGCGAGAGTTGTCAAAAAATCCGAGGGGGGAACGGAAGAATGCCGGCCCTCCCCGGTTCCATGGGGGGGGACGCCCCCTGCCGGTCCTATCTCATTTGATTGGGGGGGAAATTGCGCCGGGCCAAGGGTCCATTGGCCACCAAGAATTCTCCCAGGGGGGTTCGATGGAAATTCACGGCCACCAAGGGGGGAAAGCCCCTGGGGAGACCAAATTAAACCGGTCCTCCGGTTTAATTCGATCGTTTTTTCGGGGGCTTAAAAAGGAATCCTCCCGAAGGTCACCTCGGAACCCTGGTTAG
    Z78509.1
    TCCTAAAATTAACGGGTGCCTTAACTTGCCGGGGTTAACCAAATGCCCGTCACCCCCCATTGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAAAACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAACCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCTGCTTATCACATTTAAACCTTGTTACGTCCGTTGCCGAGAGTTGTCAGAAAAATTCATGGGGGGGCACGGAGCATGCCGGCCCTCCCCGCTCCATGGAGGACGACGCCCTCCGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78508.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAGTGGCCGTCACCGCCCATGGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAATTCATTGGGGGCACGGCAGCATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78507.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAATGGCCGTCGACCGCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGTGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAAATTCATTGGGGGACGGAAGCACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGAAATTCACGGTCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78506.1
    TCAGCGGGTGGCCTCACCTGACTGGGTTGGCATCCAATGGCCGTCAAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTATGAGCTGAGGGCTCCTTGAACGAGAGTTGTCAGCCCGATTCAAAGGGGGTACGGCGGCATGCCGGTCCTCCCCGGTTCCATGGAGGACGAAGCCCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78505.1
    AAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGACCGAACAGTATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTGACGAGAGTTGTCACAAAGATTCATTGGGGGTATGGTGGTATGCCGGGCCTCCCCGGTTCCATGGAGGACGAAGACCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCACGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAATATCATGGTCACCGCAGTAGCCACCGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78504.1
    TTACTCAGCGGTGGCTCACTGACTGGGGTTGCATCCAATGGCCGTCACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGTTTTGTCGCAGCCACCGTCCGTCCACCTCTTGTCGGATTTGGTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGTCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAGATTCATTGGGGGCACGGCGGGCATGCCGGCCCTCCCCGGCTCCATGGAGGGACGATGCCCTCTGACGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78503.1
    TTATCTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAAAGAGCTGAGACTCCGTTGCCGAGAGTTGTCAAAATAATTCATTGGGGGTACGGTAGCATGCCGGCCCTCCCCGGTTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGTAAATCACGGACACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAAGGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78502.1
    GCGTAAACTCACGGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTACGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTTCGCTGCGTCCTTCCATCGATGCAAGAGCCGAGATTCCGGCCGAGAGTTGTCAATATAAAATTCATTGGGGGGACGGTAGTATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGACCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGTCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGGGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78501.1
    TCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCTCGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACAACACTTATCACAATTTCGCTGCGTCCTTCCATCCGATGCAAGGAGCCGAGATTTCCCGTTTGGCCGAGAGTTGCCCAAAAAAAAATTCATTGGGGGCACGGCAACACGCCGGCCCTCCCCGGCTCCATGGAGGACGATTCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTTCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78500.1
    CTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCGAGATTCCCGTTGGCCGAGAAGTTTTCCCAAAGAAAAATTCATTGGGGGCACGGCAGGACGCCGGCCCTCCCCGGCTCCATGGAGGGCGATGCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCAATTATTTTGCGGTCTCAACAATGAGCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78499.1
    GGTTGCTCACCTGACCTGGGGTCGCATCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACACTTATCGCATTTCGCTGCGTCCTCCATCCGATGCAAAGAGCCGAGATATCCGTTGGCTGAGAGTTGTCAAAAAAATAATTGGGAGAGGTCAAGGCCACATGCCGTCCCCTCCTCGAATTCAACGAAGGGGGTATGTCCTTCACCATTTATGTGTCATATGACTTGGTGCAAAACTGCGACGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGAAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGCGATCTCAACAATGATCCCTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78498.1
    GCTTAAACTCAGCGGGTTGCTCACTGACTGGGGTCGCAACAAATGGCCATCAACTGATCACGGGTTGATGTGCCTCCAGTGGGGTTCACAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAGTTTTGTTGTAATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGTACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCCAAAAATAATTTGGAGAGGGGGAGTGATTGTAGGGTCAAGGGCAAAATGCCGCCCCCTCCTTCGCCATTATTGTGTCATATGACTTGGGGCAAAACTGCCCCGGGCAAGGGTAAGATCGCCGGCAAGAAAACTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCACCCATGGGTGACCGGAGTAAACTGATCCTCTGGAGTAACTCGATCAATTATTATGTGATCTCAACAATGACCTTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78497.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTGTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCGCGTTGCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACACCCGCACCGGCTGAACAGTATGCCTGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAGAGCTGAGATCTCCGTTGCTGAGAGTTGTTTCACAAAAAATAATTTGGAGAGGGGGAGGGAGTGTAGGCTCAAGGCCACATGCCCTTCACCTCCTTGAATTCACAAAGGGCCGTGCCCTTCACCATTTATGTGTCATGTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCGCCCATGGGCGACCAAAGTAAACTGATCCTCTGGATTCACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78496.1
    GCTTAAACTCAGCTGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGGTGCGTCCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTCAAAAAATTAATTGGAGAGGTCAAGGCCACATGCCGGCCCCTCCTCGAATTCACGAAGGGATATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGCAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78495.1
    CACTCAGCGGGTTGCTCACTGACCTGGGGTCGCAACCGAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTGGATTATGCCGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAACAGTGTTGTTGTAGCCTGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCTGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTGAGGGAGAGAGATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTCATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGGATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTCAGAAAAATACTTTGGAGAGGGGGAGAGCGAGTGCAGTGTCAAGGCCACATGCCGCCCCCCCTCCTTGAATACACGAAGGGCTATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACCGCGCCGGACAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCTCCATGGCAACTCTAGGTCACCGCAATAGCAAGCGCCCATGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78494.1
    CTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTACGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGACAGGCGTTCCCTTGGCCTGTTGGCCTCGGGCCCAACTTGCGTTCAAAGACTCGATGTTAACGGGATTCTGCAATTCACACCTGAGATATCCGAAGTTGAGAGTGTTACCTCAATAATGGGGAGTTGGTCAAGGCAGTATGCCGTCCCCCTTCATTAATTATGGGTCATTTGATTTGGCGCATTACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGTTCCCAAGGAGGCCCGATGGTAAATCTAGGTCACTGCAACAGTAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGACCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78493.1
    GGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCACACCGGATGACCAGTACGTTTGGACAGCCTCATTAAGGGAGAGATATGAGTCGCAATGTCCAGGCAGGCGTCCCCTTGGGCTGATGGCCTCGCGCGCAACTTGCGTTCAAAGATTCGATGGTTCACGGGATTCTGCAAGGCACACCATTTATCGAATCTCGCTGCGTTCCTTCATCGATGCATGAGCCGAGATATCCGTTGCTGAGAGTTGTTACCTAAATACTTGGGGGCTGGTCAAGGAAGAATGCCGCCCCCCTTCACCACTTATGTAAAATTTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGCTCCCAGGGAGGCCCGATGGAAATTCTAGGTCACTGCAACAGAAAATGCCCATGGGTGACCAAAGTAAACCGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78492.1
    TATTAAACTCAGCGGTGTGCCTCACCTGACTGGGATCGCAACCAAATGGCCAATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCTGATTAAGCTGGCCCATGTAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAATTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCGTGCAGCTCTTAGACCCCCCGTCCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTCCACGGGATTCTGCAAATCACACCAGTATCGCATTTCGCTGGGTTCTACATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGAGTGGGTCAAGGCAGTACGCTCCCCCCCTTCACCAATTATGTGACATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78491.1
    GCTTATCTCAGCGGGTTGCTCACCTGACTGGGATCGCAACAAATGGCCATTCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTCCAAAGGGTCTTCTGATTATGCTGGTCCATGTAACACAACAACCCGGGGCATTCGCACAACATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAATATAATTTGGAGCGGGTCAAGGCAGCACGCCTCCCCCCAACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCTAAGGAGGCCCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78490.1
    TCAGCTGGTTGCTCACCTGACTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGACGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAAATTATGCTGGCCCATCTAATACGACAACGCGGGGCATTCGCACAACACTTGTTGCAGAACAGTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCGTTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCGAGGCAGCATGCCGCCCCCTTCACCAATTATGTGCCATATGACTTGGCGCAAAACTGCGCCGGACTAGGGTTCAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78489.1
    GCCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCATGGGGTTCAAAAGGGTCTTCAGATATGCTGGGCCCATCTAATACGACAACTCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAATCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACCCACATCCCCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATCTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTCGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCTCGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78488.1
    AGCCGGTTGCTCACCTGACTGGGGTCGCAACAAATGTCCATCAACTGATCATGGGTTGATGGGCCTCCGATGGGGTTCAAAAGGGTCTTCAGATTATGATGGCCCATCTAATACGACAACCCGGGGGATCCGCACAACAGTTTGGTTGTGGCGTTGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGGCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTCCATCCCGTTGATGAGAGTTGTTAGAAAATAATTAGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCACCAATTATGTGTCGTATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTAAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAACAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTGCGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACAG
    Z78487.1
    TTACTCAGCGGTTGCTCACCTGACTGGGGTCGCAACAAGTGGCCCCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTACGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGCATTTAGGACCATCGATAGCCCATGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCGAACTCACATCCGCACCGGCTGAACAATATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTAAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGAGTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78486.1
    TCAGCGGGTTGCCTCACCTGACCTGGGGTCGCTACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGGCACTGCAATAGCAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATAATCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGATTCACCTACGGAAACCTCGTGACG
    Z78485.1
    TACTAAATCAGGGGTTGCTCAGCTGACTGGGGTCGCAACACATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAGGGGTCTTTAGATTGTGCTGGCCCGTCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGGCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGCCCAAACTCACATCCGCACCGGCTGCACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCCCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGATATCTCAACAATGATCCATCCGCAGATTCACCTTCGGACACCAGGTTCAG
    Z78484.1
    AAAACTCAGGGAGTTGCTTCACCTGACTTGGGGTCGCAACCCAATGGACCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTGTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGGGAGGGTAGGAAACATGCCGCCCCTTCACAATTATGTGTCATATGACTTGGGCCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGGAGGCTCGATGGTAAACTCTCGGTCACTGCAATAGCCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCCCAGGTTCACCTACGGAAACCTTGTTACG
    Z78483.1
    TGCTAAACTCAGGGGTTGCTTCACTTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGATATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATTGTTCACGGGATTCTGCAATTCACACCATTTATGATTTCGGTTTGGTCTTCATCCATGCAAGAGGCCAGATTTCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78482.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACGCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGNTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGAAGAGGCGAGATATCCGTTGCNGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTACGGTTTAGATCGCCAGCAAGAAAGCTCCCAGGAGGCTCGATGGCAAATCTCGGTCACTGCAGTAGA
    Z78481.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTCCGTTGCTGAGAGTTGTTAAAAAAATGATTTGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78480.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCGCACAACATTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78479.1
    ACTCAGAGGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCGTCAACTGATCTTGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAGACTCGATGGTTCACGGGATCTCTGTATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGGCAAGGCAAGGCAGCATCCCTCCCCTTCCAGTAATGTCTCATATAAATTGGCGCAAATCTGCGCCGGGCAAGGGTTTAGTTCGGCTCGATGGCAAATCTAGGTCACTGAAACAGCAAATGCCCATGGGTGACCAAAGTAACCTGATCCTCCTGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78478.1
    GCCTAAACTCAGCGGGTGCCTCATCTGACCTGGGTCGCAACCAAATGTCCTCAACTGATAATGGGTTGATGGGCCTCCAATGGGGGTCAAAAGGGTCTTTAGATTATGCTGACCCATCTAAACGACAACCGGGGCATTCGCACAACAGTTCTGTTGTCACATAGCTCGTCCACCTCTTGCCGTATTTAGGACCATCCAACCCATACAGCTCTTAGACCCCCGTACCAAAGGACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGCAGCCTCATTAAGGGAGAGATATGTCTCACAATTTCCAGGCAGGCGTCCCCTTGACCTGATGGCCTCGGACGCAACATGCGTTCAAAACTCGATGTTCAGGGATCTGCTATCACGCCATTATCAATATTTGGGGGAGCCATGGCAGCATGCCGCCCCTTCCAGTTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAACTCTAGGCCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCCCCAGATTAACTCGATCAATTATTATGTGATCTCAACACTGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78477.1
    GCATAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCCGCAACTTGCGTTCCAAAGACTCGATGGTTCAGGGATTCTTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAGAGCCGAGATACTCGTTGCTGAGAGTTGGTTAACAAATATTTGGGGAGGGCAAGGCAGCATGCCGCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGCAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78476.1
    GGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGATTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCGCAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTGGATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCGGATGGCCTAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGCCAAGGCAGCATGCCGCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78475.1
    ACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCATTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGACGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78474.1
    AAGCTCAACGGGTTGCTTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAGAAGGGTCTGTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGTGTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTCAAAAATATAATTTGGGGAGGGTCAAGGCAGGATGCCGCCCTTCCAATTATGTGTCATATGATTGGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGAAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTACGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78473.1
    CCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGGCAAGAACAAGGGGCCAAACTCACATCCGCACCGACTGAACAGTATGTATGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAACAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCACCAAGAAAGCTCCCATGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78472.1
    GCTTAAACTCAGCGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATTGGCCTCTAAAGGGGTTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGAACTACAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAAGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATAGACAGCCTCATTAAAGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATTCCGTTGCTGAGAGTTATCAAAAAATAATTTGGGGAGGGTCTTAGACTGCATGCCGCCCCCTTCCAATTATGTGTGAAATGACTTGGCGCAAGACTGCGCCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78471.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTACTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGATGTTGAGAGTTGTCAAAAAATTAATTTGGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCATTATGTGTCGGGGTGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGACAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACCGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78470.1
    AACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTAGTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTAAAAAAATAATTTGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78469.1
    AACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGCACTATCATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCTAGACCGCATGCCGCCCCCTTCCAATTATGTGTGATATGACTTGGCGCAAGACTGCGCCGGGCAACGAATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78468.1
    AACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCACACAATAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTCCATCAAAAGCCCAGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGAGTTCCGGAGCTGAGAGTCGTTAAAAAAATGATATGGGGAGGGTCAAGGCAGCATGCTGCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTCAGATCGCCATCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATTCCCATGAGTGACCAAAGTAAACTGATCCTCCAGATTAACTCAATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78467.1
    TCAGCGGGTTGCCTCACCTGACCTGGGTCGCAAACATATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAGAGGGTCTTTAGATTATGCTGGCCCAGCTAATACGACAACCCGGGTCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCACCAAAGGCACATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTACCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78466.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACGCAAGAACAAGGGGCCAAACTCACATCCGCACAGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCAAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTACAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGAAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78465.1
    GCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACTCGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTTCGTTCAAAGACTCGATGGTTTACGGGATTCTTCAATTCACACCAATTATCGCATTTCGCTGCGTTCCTCATCGATTCAAGAACCGAGAAATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78464.1
    GCTTAAACTCAGCGGGTTGCTCACCTGACCTGGGGTCGCACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATTACGAAAGCCCGGGGAATCCGAACAGAAATTTGGTTGAAGAATAGTCCGCCCACCTCTTGCCGTTTTTAGGACCACCAAAAGCCCATGAAGTTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGAACCGGTTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATTTGATTCGAAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGAACTTGGCGTTCAAAGACTCGATGGTTCACGGGATTCTGAAATTCACACCATTTATCGGATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCTCCGGGCAAGGGTTTAGATCTCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGATCAATTATTATGTGATCTCAACAATGACCCTTCCGCTCACCTACGGAAACCTTGTTACG
    Z78463.1
    GCTTAAACTCAGCGGGTTGCTCACCTGATCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGTATTCGCACAGCAATTTTGTTGCAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAACATGCCGCCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTCCCCCGGGCGAGGGTTTAGATCCCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGAACAATGATTATGTGATCTCAACAATGAACCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78462.1
    ATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTGCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGGTCACATCCGGAGACCTCGTGACG
    Z78461.1
    TTAACTCAGCGGCTTGCTCACTGACTGGGTCGCAACAATGGCATCAACTGATCATGGGTTGATGGGCTCCAATGGGGTTCAAAGGGTCTTCAGATTACGGTGGCCCATCTTATACGACAACCCGGGGGCTTTGGACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGGATGTATGGACAGGCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGGAATTCACACCATTTATCGGATTTCGGTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTGAAAAAATTATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGTCATTTGACTGGGGGAAAAATTGCGCCGGGTAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTCCAGCAGCAACTGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78460.1
    TAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTGGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78459.1
    AAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAATTTATCGCAATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78458.1
    CAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCGTACGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78457.1
    CTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACGATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGGTCCTTCATCCGATGAAAGAGCCGAGATATCCGTTGTTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGCCATTTGACTGGGGGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTACAACAGAAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78456.1
    GCTAAACTCAGGGTTGCCTCACCTGACCTGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGACCGGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTCCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78455.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGAAGCTCTTAGACCCCCCGTGCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCAACCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGGATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGAAAGGGTTTAGTTCTCGCCAAAAAGAAAGTTCCCAAGGGGGTTCGATGGAAATTCTAGGTCACTCCAAAAGAAATTGTTCATGGGTGACCAAAGTAAAAAGTTCCTCCAGATTAACTCGATCAATTTTTATGTGATCTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTGGTTACG
    Z78454.1
    GTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGNGGTTGATGGGCCTCCAATGGGTTCAAAAGGGGCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAGAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGTCTGATGGCCTCGGCCGCAACTTGCGTTCACAGACTCGATGGTTCACGGGATTCTGCAAATCACACCATTTATCGCATGAGTTCCGTTGATGAGAGTTGTTAACAATAGTGGGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCGATTATGTGTCAAATGACTTGGAGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGTAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78453.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAGTTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78452.1
    TGCTAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACAGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCATCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATTCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAAAAAGAAAGCTCCCAAGGAGGCTTGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAATAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78451.1
    GCTCAGCTGGTTGCTCACTGACTGGTGTCGCATCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCTATAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGGCTTCGCACAACAATTTTATTGTAGGGATAGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAATCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGGTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATTTGACTCGCGATGCCCAGGGAGGGGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGATTTCGCTGGTTCTTCATCGATGAAAGAGCGAGATTTCCGTTGTTGAGAGTTGCAAAAAATACTTTGGGGAGGGTCAGGGCAGGATGCCCGCCCCCTACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCACGGGTTTAGATCTCGCCAGCAAGAAAGCTCCCAGGGGGGCTCCATGGAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTACACCTACGGAAACCTTGTTACG
    Z78450.1
    CTCAGGGGTTGCTCAGCTGATCTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGTCCTCCAATGGGGTTCAACAGGGTCTACAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGCACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCAGTTTTGAGGACCATCATATGCCCATGCAGTTCTTAGACCCCCCGGACCAAAGAACAAGGGGCCACACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAGGGGAGAGATATGACTCGGAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGGCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGTAAATGCTCATGGATGACCAAAGTAAACAGATCCTCCAGCTTAACTCGATCAATTATTATGTGATATCAGCAATGATCCTTCC
    Z78449.1
    GCATAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTAGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGGAGGGTCAGGGCAGGATGCCACCCTTCACCAATTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78448.1
    CCTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCANCCAAATGGCCATCAACTGATCATGGGCTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTNTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGCGTTCTTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAAGGGAAGGATGCCACCCCCCTTCACCAATTATGTGTCATACGACTTGGCGCAAAACTGCGCCGGGAAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78447.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATAGGATGCCACCCCCTTCACCACTTATGTGTCATTTGACTGGGGGCAAAATTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAAAGCTCCCAGGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGGAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78446.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGAGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGGCCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCTGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTATCGCAATTTCGCTGGTCCTCCATCGATGAAGAGCCGAGTATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGGAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78445.1
    ACAGCTGTTGCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTGCACGGGATTCTGCAATTCACACCATTTATCGCATTTGCCGAAGACTGAGTCTCCGTTGCTGAGAGTTGTTAAAAAAATAGTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCAAACGACTTGGAGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78444.1
    AATGGCCATCAACTGACATGGGTTGATGGGCTCCAATGGGGTCCAAAAGGGTCTTCAGATATCGGTGGCCCATCTTATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCATCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTAAGGGAAGAATGCCACCCCCTCCACCAATTTTGTGTAATTTGATTGGGGGAAAAATTGAGGGTTTAGATATCGCCAACAAGGATGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGTTCACCCTACGGAAACCTTGTTACG
    Z78443.1
    CCTCACCTGACTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTATGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGGTGTGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTGGCCCAAAACTCTGCCGGGCAAGGGTTTAGATATCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78442.1
    ACTCAGCGGGTTGCTCAGCTGACTGGGGTCGCAACAAATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTTGTTGTGGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCGGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTTTGGATAGCCTCGTTAGGGGAGAGATATGATTCCCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCAAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78441.1
    CTCAGCGCGTTGCTCAGCTGACTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCACAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAACTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATATCGGCAATGACCTTCC
    Z78440.1
    TGCTAAACTCAGGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGCTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCAACGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGATATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAATGCTCCCAAGGAGGGTCGATGGAAATTCTAGGTCACTACAACAGAAATTGCCCATGGGTGACCAAAATATGCAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGGAGGTCCACCTACGGAAACCTTGTTACG
    Z78439.1
    GGCCCAACTAAACGGCCAACCGGGGCATTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTTCCGTATTTAGGACCATCCAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATTACCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCCAACTCCGCCGGGCAGGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATG



```python
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
]
```


```python
len(records)
```




    94




```python
records = [
    rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
]
```


```python
len(records)
```




    18




```python
records = (
rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
)
```


```python
records = (
rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
)

SeqIO.write(records, "rev_comp.fasta", "fasta")
```




![image](https://github.com/jameciajones/Python_Portfolio/assets/134228924/2e427311-27da-40b3-b03f-3507a214a066)

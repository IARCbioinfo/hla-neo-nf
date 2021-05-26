
# alleles

## match xHLA to pvactools

A total of 10230 alleles are available in xHLA, of those 28 are not valid in pvactools.
The file xHLA2PVAC_alleles.txt convert xHLA to pvactools alleles. 

### command
```
perl match_xHLA_to_pvactools_alleles.pl -a pvacseq_valid_alleles.txt -b xHLA.all.alleles.txt >  xHLA2PVAC_alleles.txt
```

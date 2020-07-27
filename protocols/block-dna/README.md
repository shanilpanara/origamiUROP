# Generation of Block Copolymer DNA pieces

This folder contains code to generate block copolymer ssDNA/dsDNA.

## Usage
To create a system containing a strand of 30bp with double-stranded parts of length 10bp and single-stranded parts of length 5bp, written to the files `oxdna.out.conf` and `oxdna.out.top`, in a box of size 50.0, run the following command:

```
python fiona_gen.py -n 30 -ds 10 -ss 5 -f out -b 50
```

Using Ovito to view, we can visualise a snapshot of this starting configuration:

![Snapshot](fiona-dna.png)
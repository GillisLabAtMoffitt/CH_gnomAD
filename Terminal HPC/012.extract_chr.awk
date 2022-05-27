{ if (NR > 1) { 
    # all data lines 
    split($39, a, ":");  
    print a[1], $0; 
  } else {
    # header line
    print "#X.CHROM", $0;
  } 
}

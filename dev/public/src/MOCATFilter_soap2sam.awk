#!/usr/bin/awk -f

# This code is part of the MOCAT analysis pipeline
# Code is (c) Copyright EMBL and BGI, 2012-2013
# This code is released under GNU GPL v3.

# USAGE:  <some command> | ./soap2sam.awk -v MIN_LEN=30 -v MIN_AS=90 
BEGIN { 
    OFS = "\t";
    total = 0;
    lengthfilter = 0;
    asfilter = 0;
}
{
  total = total + 1;
  if ($6 < MIN_LEN) next;
  lengthfilter = lengthfilter + 1;

  
  strMM = $NF;
  mm = 0;
  n = length(strMM);
  
  for (i = 1; i <= n; i++) 
    {
      c = substr(strMM, i, 1);
      if (index("ATCG", c) != 0) mm++;
    }
  
  as = 100-(mm/$6)*100;
  if (as < MIN_AS) next;
  asfilter = asfilter + 1;
  
  mapq = int(as * 10); 
  
  if ($7 == "+") flag = "0"; else flag = "16";
  
  print $1, flag, $8, $9, 255, $(NF-1), "*\t0\t0", $2, $3, "AS:i:" mapq, "NM:i:" mm;
  
}
END {
  OFS = "";
  print "XXXXXXXX XXXXXX MOCATFilter - length & quality filter : [STATS] total_reads_in=", total, " | total_reads_after_length_filter_only=", lengthfilter, " | total_reads_out=", asfilter >"/dev/stderr"; 
}


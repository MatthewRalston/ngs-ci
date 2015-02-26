Library Complexity Index

* Summary
This is a Ruby gem that calculates an alternative to the "pileup" format. The calculation is an average of average overlaps among reads at a particular base in the genome. Offers both stranded and unstranded calculations.

* To-do list
1. Create gem environment
2. Create master class with options a la transrate
   * including FR, RF, ??, or F strand specific options
3. Create bam processing class that
   * increments through bases (x)
	 * calls stranded or unstranded method
	 * adds results to list
   * stranded method (strand chemistry)
	 * calls methods by strand chemistry
	 * returns calculation class results
   * strand chemistry methods
     * each has specific SAM flags used to acquire reads
     * e.g. for FR chemistry
	     * F reads are acquired according to strand
		 * R reads are acquired and assigned according to strand of mate
   * unstranded method
	 * calls samtools to acquire reads from base "x"
	 * converts to bed and sorts
	 * returns calculation class results
4. Create calculation class that
   * Increments through reads (i)
	 * Acquires reads overlapping read "i"
     * Calculates average overlap and adds to list
   * Averages all overlaps and returns
6. Prints either a Nx2 matrix or 1xN matrix to file.

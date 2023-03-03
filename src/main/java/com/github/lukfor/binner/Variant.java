package com.github.lukfor.binner;

public class Variant {

	public double pval;
	
	public String chrom;
	
	public long pos;
	
	public int num_significant_in_peak;
	
	
	@Override
	public String toString() {
		return chrom + ":" + pos + " --> " + pval;
	}
	
}

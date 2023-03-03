package com.github.lukfor.binner;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

public class Bin {

	public String chrom;

	public long startpos;

	public Set<Double> qval = new TreeSet<Double>();

	private double size = 0.1;

	public List<Double[]> getLines() {
		List<Double[]> lines = new Vector<>();
		// int start =

		Double first = qval.iterator().next();
		Double[] line = new Double[] { first, first };
		lines.add(line);
		for (double p : qval) {
			if (line[1] + size * 1.1 >= p) {
				line[1] = p;
			} else {
				line = new Double[] { p, p };
				lines.add(line);
			}
		}

		return lines;
	}

}

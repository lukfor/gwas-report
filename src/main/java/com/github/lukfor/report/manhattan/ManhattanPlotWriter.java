package com.github.lukfor.report.manhattan;

import java.io.File;

import com.github.lukfor.binner.Bin;
import com.github.lukfor.binner.ChromBin;
import com.github.lukfor.binner.Variant;

import genepi.io.table.writer.CsvTableWriter;

public class ManhattanPlotWriter {

	public static void saveAsFile(ManhattanPlot data, File file) {

		CsvTableWriter writer = new CsvTableWriter(file.getAbsolutePath());
		writer.setColumns(new String[] { "CHR", "BP", "P", "peak" });

		int unbinned = 0;
		int binned = 0;
		int peak = 0;

		int bins = 0;

		for (ChromBin chromBin : data.getBins().values()) {
			
			for (Variant variant : chromBin.getUnbinnedVariants()) {
				writer.setString("CHR", variant.chrom);
				writer.setInteger("BP", (int) variant.pos);
				writer.setDouble("P", variant.pval);
				writer.setString("peak", "no");
				writer.next();
				unbinned++;
			}

			for (Variant variant : chromBin.getPeakVariants()) {
				writer.setString("CHR", variant.chrom);
				writer.setInteger("BP", (int) variant.pos);
				writer.setDouble("P", variant.pval);
				writer.setString("peak", "" + variant.num_significant_in_peak);
				writer.next();
				peak++;
			}
			
			for (Bin bin : chromBin.getBins().values()) {
				bins++;
				for (double p : bin.qval) {
					writer.setString("CHR", chromBin.chrom);
					writer.setInteger("BP", (int) bin.startpos);
					writer.setDouble("P", p);
					writer.setString("peak", "no");
					writer.next();
					binned++;
				}
			}
		}

		writer.close();

		System.out.println("Wrote to file " + file.getAbsolutePath() + " Unbinned: " + unbinned + ". Peaks: " + peak
				+ ". Binned: " + binned + " in " + bins + " bins");

	}

}

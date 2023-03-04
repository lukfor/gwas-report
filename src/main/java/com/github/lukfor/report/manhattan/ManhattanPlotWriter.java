package com.github.lukfor.report.manhattan;

import java.io.File;

import com.github.lukfor.binner.Bin;
import com.github.lukfor.binner.ChromBin;
import com.github.lukfor.binner.Variant;

import genepi.io.table.writer.CsvTableWriter;

public class ManhattanPlotWriter {

	public static void saveAsFile(ManhattanPlot data, File file) {

		CsvTableWriter writer = new CsvTableWriter(file.getAbsolutePath());
		writer.setColumns(new String[] { "CHR", "BP", "P", "type", "y1", "y2" });

		int unbinned = 0;
		int binned = 0;
		int peak = 0;

		int bins = 0;

		for (ChromBin chromBin : data.getBins().values()) {

			for (Variant variant : chromBin.getUnbinnedVariants()) {
				writer.setString("CHR", variant.chrom);
				writer.setInteger("BP", (int) variant.pos);
				writer.setDouble("P", variant.pval);
				writer.setString("type", "unbinned");
				writer.next();
				unbinned++;
			}

			for (Variant variant : chromBin.getPeakVariants()) {
				writer.setString("CHR", variant.chrom);
				writer.setInteger("BP", (int) variant.pos);
				writer.setDouble("P", variant.pval);
				writer.setString("type", "peak");
				writer.next();
				peak++;
			}

			for (Bin bin : chromBin.getBins().values()) {
				bins++;
				for (Double[] line : bin.getLines()) {
					writer.setString("CHR", chromBin.chrom);
					writer.setInteger("BP", (int) bin.startpos);
					writer.setDouble("P", -1);
					writer.setString("type", "bin");
					writer.setDouble("y1", line[0]);
					writer.setDouble("y2", line[1]);
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

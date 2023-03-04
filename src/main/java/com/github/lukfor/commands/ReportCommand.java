package com.github.lukfor.commands;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import com.github.lukfor.App;
import com.github.lukfor.binner.Binner;
import com.github.lukfor.binner.Variant;
import com.github.lukfor.report.Report;
import com.github.lukfor.report.manhattan.ManhattanPlot;

import genepi.io.table.reader.CsvTableReader;
import picocli.CommandLine.Command;
import picocli.CommandLine.Help.Visibility;
import picocli.CommandLine.Option;

@Command(name = "report", version = App.VERSION)
public class ReportCommand implements Callable<Integer> {

	@Option(names = { "--input" }, description = "Input filename", required = true)
	private String input;

	@Option(names = { "--chr" }, description = "Chromosome column in input file", required = true)
	private String chr = "CHROM";

	@Option(names = { "--position", "--pos" }, description = "Position column in input file", required = true)
	private String position = "GENPOS";

	@Option(names = { "--pvalue", "--pval" }, description = "PValue column in input file", required = true)
	private String pval = "LOG10P";

	@Option(names = { "--title" }, description = "Custom title of report", required = false)
	private String title = null;

	// TODO: add flag to use p-values (not -log10)
	// --neg-log-pvalue
	private boolean negLogPvalue = true;

	// TODO: add flag --optimize
	private boolean optimize = true;

	// TODO: add flag --annotations
	private boolean annotations = true;

	@Option(names = {
			"--sep" }, description = "Separator of input file", required = false, showDefaultValue = Visibility.ALWAYS)
	private char separator = '\t';

	@Option(names = {
			"--ref" }, description = "Ref allele column in input file", required = false, showDefaultValue = Visibility.ALWAYS)
	private String ref = null;

	@Option(names = {
			"--alt" }, description = "Alt allele column in input file", required = false, showDefaultValue = Visibility.ALWAYS)
	private String alt = null;

	@Option(names = {
			"--rsid" }, description = "RsID column in input file", required = false, showDefaultValue = Visibility.ALWAYS)
	private String rsid = null;

	@Option(names = {
			"--gene" }, description = "Gene column in input file", required = false, showDefaultValue = Visibility.ALWAYS)
	private String gene = null;

	@Option(names = { "--output" }, description = "Output filename", required = true)
	private String output;

	public void setInput(String input) {
		this.input = input;
	}

	public void setChr(String chr) {
		this.chr = chr;
	}

	public void setPosition(String position) {
		this.position = position;
	}

	public void setPval(String pval) {
		this.pval = pval;
	}

	public void setOutput(String output) {
		this.output = output;
	}

	public void setNegLogPvalue(boolean negLogPvalue) {
		this.negLogPvalue = negLogPvalue;
	}

	public void setOptimize(boolean optimize) {
		this.optimize = optimize;
	}

	public void setAnnotations(boolean annotations) {
		this.annotations = annotations;
	}

	public void setTitle(String title) {
		this.title = title;
	}

	@Override
	public Integer call() throws Exception {

		System.out.println("Process file '" + input + "'...");

		CsvTableReader reader = new CsvTableReader(input, separator, true);

		Binner binner = new Binner();

		int count = 0;
		while (reader.next()) {
			count++;
			Variant variant = new Variant();
			variant.chrom = reader.getString(chr);
			variant.pos = reader.getInteger(position);
			if (negLogPvalue) {
				variant.pval = reader.getDouble(pval);
			} else {
				variant.pval = -Math.log10(reader.getDouble(pval));
			}
			// TODO: load rsid and gene when set

			binner.process_variant(variant);
		}

		// TODO: last peak/bin is not closed...
		// binner.stop();

		reader.close();

		System.out.println("Processed " + count + " variants.");

		ManhattanPlot data = new ManhattanPlot(optimize);
		data.setBins(binner.getBins());
		data.setShowAnnotations(annotations);
		data.setPeaks(new ArrayList<Variant>(binner.getPeaks()));
		data.setUnbinnedVariants(new ArrayList<Variant>(binner.getUnbinnedVariants()));

		Report report = new Report(data);
		if (title != null && !title.isEmpty()) {
			report.setTitle(title);
		}
		report.saveAsFile(new File(output));

		return 0;
	}

}

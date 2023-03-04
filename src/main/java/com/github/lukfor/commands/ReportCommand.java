package com.github.lukfor.commands;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import com.github.lukfor.App;
import com.github.lukfor.binner.Binner;
import com.github.lukfor.binner.Variant;
import com.github.lukfor.report.Report;
import com.github.lukfor.report.manhattan.ManhattanPlot;
import com.github.lukfor.report.manhattan.ManhattanPlotWriter;
import com.github.lukfor.util.AnnotationType;
import com.github.lukfor.util.OutputFormat;

import genepi.io.table.reader.CsvTableReader;
import genepi.io.table.reader.ITableReader;
import picocli.CommandLine.Command;
import picocli.CommandLine.Help.Visibility;
import picocli.CommandLine.Option;

@Command(name = "report", version = App.VERSION)
public class ReportCommand implements Callable<Integer> {

	@Option(names = { "--input" }, description = "Input filename", required = true)
	private String input;

	@Option(names = { "--chr" }, description = "Chromosome column in input file", required = false)
	private String chr = "CHROM";

	@Option(names = { "--position", "--pos" }, description = "Position column in input file", required = false)
	private String position = "GENPOS";

	@Option(names = { "--pvalue", "--pval" }, description = "PValue column in input file", required = false)
	private String pval = "LOG10P";

	@Option(names = { "--title" }, description = "Custom title of report", required = false)
	private String title = null;

	// TODO: add flag to use p-values (not -log10)
	// --neg-log-pvalue
	private boolean negLogPvalue = true;

	// TODO: add flag --optimize
	private boolean optimize = true;

	@Option(names = {
			"--annotation" }, description = "Show annotation labels in plot", required = false, showDefaultValue = Visibility.ALWAYS)
	private AnnotationType annotation = AnnotationType.NONE;

	@Option(names = {
			"--sep" }, description = "Separator of input file", required = false, showDefaultValue = Visibility.ALWAYS)
	private char separator = '\t';

	@Option(names = {
			"--ref" }, description = "Ref allele column in input file", required = false, showDefaultValue = Visibility.ALWAYS)
	private String ref = "ALLELE0";

	@Option(names = {
			"--alt" }, description = "Alt allele column in input file (effect allele)", required = false, showDefaultValue = Visibility.ALWAYS)
	private String alt = "ALLELE1";

	@Option(names = {
			"--beta" }, description = "Beta column in input file", required = false, showDefaultValue = Visibility.ALWAYS)
	private String beta = "BETA";

	@Option(names = {
			"--rsid" }, description = "RsID column in input file", required = false, showDefaultValue = Visibility.ALWAYS)
	private String rsid = "ID";

	@Option(names = {
			"--gene" }, description = "Gene column in input file", required = false, showDefaultValue = Visibility.ALWAYS)
	private String gene = null;

	@Option(names = { "--output" }, description = "Output filename", required = true)
	private String output;

	@Option(names = {
			"--format" }, description = "Output format", required = false, showDefaultValue = Visibility.ALWAYS)
	private OutputFormat format = OutputFormat.HTML;

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

	public void setSeparator(char separator) {
		this.separator = separator;
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

	public void setAnnotation(AnnotationType annotation) {
		this.annotation = annotation;
	}

	public void setRsid(String rsid) {
		this.rsid = rsid;
	}

	public void setGene(String gene) {
		this.gene = gene;
	}

	public void setRef(String ref) {
		this.ref = ref;
	}

	public void setAlt(String alt) {
		this.alt = alt;
	}

	public void setBeta(String beta) {
		this.beta = beta;
	}

	public void setTitle(String title) {
		this.title = title;
	}

	public void setFormat(OutputFormat format) {
		this.format = format;
	}

	@Override
	public Integer call() throws Exception {

		System.out.println("Process file '" + input + "'...");

		CsvTableReader reader = new CsvTableReader(input, separator, true);

		if (!reader.hasColumn(chr)) {
			System.out.println("Error: Column '" + chr + "' not found.");
			return 1;
		}

		if (!reader.hasColumn(chr)) {
			System.out.println("Error: Column '" + position + "' not found.");
			return 1;
		}

		if (!reader.hasColumn(chr)) {
			System.out.println("Error: Column '" + pval + "' not found.");
			return 1;
		}

		rsid = checkColumn(reader, rsid, "rsid");
		beta = checkColumn(reader, beta, "beta");
		gene = checkColumn(reader, gene, "gene");
		alt = checkColumn(reader, alt, "alt");
		ref = checkColumn(reader, ref, "ref");

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
			if (rsid != null && !rsid.isEmpty()) {
				variant.id = reader.getString(rsid);
			}
			if (beta != null && !beta.isEmpty()) {
				variant.beta = reader.getString(beta);
			}
			if (gene != null && !gene.isEmpty()) {
				variant.gene = reader.getString(gene);
			}
			if (alt != null && !alt.isEmpty()) {
				variant.alt = reader.getString(alt);
			}
			if (ref != null && !ref.isEmpty()) {
				variant.ref = reader.getString(ref);
			}

			binner.process_variant(variant);
		}

		binner.close();

		reader.close();

		System.out.println("Processed " + count + " variants.");

		ManhattanPlot data = new ManhattanPlot(optimize);
		data.setBins(binner.getBins());
		data.setAnnotation(annotation);
		data.setPeaks(new ArrayList<Variant>(binner.getPeaks()));
		data.setUnbinnedVariants(new ArrayList<Variant>(binner.getUnbinnedVariants()));

		if (format == OutputFormat.HTML) {
			Report report = new Report(data);
			if (title != null && !title.isEmpty()) {
				report.setTitle(title);
			}
			report.saveAsFile(new File(output));
		}

		if (format == OutputFormat.CSV) {
			ManhattanPlotWriter.saveAsFile(data, new File(output));
			return 0;
		}

		if (format == OutputFormat.JSON) {
			System.out.println("Not yet implemented.");
			return 1;
		}

		return 0;
	}

	private String checkColumn(ITableReader reader, String column, String param) {
		if (column != null && !column.isEmpty()) {
			if (reader.hasColumn(column)) {
				return column;
			} else {
				System.out.println(
						"Warning: Column '" + column + "' not found. Use '--" + param + " <COLUMN>' to change it.");
			}
		}
		return null;
	}

}

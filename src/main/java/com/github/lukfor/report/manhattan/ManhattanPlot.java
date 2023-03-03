package com.github.lukfor.report.manhattan;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import com.github.lukfor.binner.Bin;
import com.github.lukfor.binner.ChromBin;
import com.github.lukfor.binner.Chromosome;
import com.github.lukfor.binner.Variant;

public class ManhattanPlot {

	private Map<Integer, ChromBin> bins;

	private String[] colors = new String[] { "#779ECB", "#03254c" };

	public static final double bin_size = 3_000_000;

	public static final int gap = 10;

	public static final int point_size = 5;

	private boolean asShapes = true;

	private boolean annotations = true;

	private List<Variant> peaks;

	public ManhattanPlot(boolean asShapes) {
		this.asShapes = asShapes;
	}

	public void setBins(Map<Integer, ChromBin> bins) {
		this.bins = bins;
	}

	public Map<Integer, ChromBin> getBins() {
		return bins;
	}

	public void setPeaks(List<Variant> peaks) {
		for (Variant peak : peaks) {
			int index = Chromosome.getOrder(peak.chrom);
			ChromBin bin = bins.get(index);
			bin.addPeakVariant(peak);
		}
		this.peaks = peaks;
	}
	
	 public List<Variant> getPeaks() {
		return peaks;
	}

	public void setUnbinnedVariants(List<Variant> unbinnedVariants) {
		for (Variant variant : unbinnedVariants) {
			int index = Chromosome.getOrder(variant.chrom);
			ChromBin bin = bins.get(index);
			bin.addUnbinnedVariant(variant);
		}
	}

	private Map<String, Object> getTrace(int chrIndex, long offset, String color) {
		Map<String, Object> trace = new HashMap<String, Object>();
		List<Double> x = new Vector<Double>();
		List<Double> y = new Vector<Double>();
		ChromBin bin = bins.get(chrIndex);
		for (Bin singleBin : bin.getBins().values()) {
			for (double p : singleBin.qval) {
				x.add((singleBin.startpos / bin_size) + offset);
				y.add(p);
			}
		}
		trace.put("x", x);
		trace.put("y", y);
		trace.put("mode", "markers");
		trace.put("type", "scatter");
		trace.put("name", bin.chrom);

		Map<String, Object> marker = new HashMap<>();
		marker.put("color", color);
		marker.put("size", point_size);

		trace.put("marker", marker);
		trace.put("hoverinfo", "none");
		return trace;
	}

	private List<Object> getShapes(int chrIndex, long offset, String color) {
		List<Object> shapes = new Vector<>();
		ChromBin bin = bins.get(chrIndex);
		for (Bin singleBin : bin.getBins().values()) {
			for (Double[] line : singleBin.getLines()) {
				Map<String, Object> shape = new HashMap<String, Object>();
				shape.put("type", "line");
				shape.put("xref", "x");
				shape.put("yref", "y");
				// + binsize / 2
				shape.put("x0", (singleBin.startpos / bin_size) + offset);
				shape.put("x1", (singleBin.startpos / bin_size) + offset);
				shape.put("y0", line[0]);
				shape.put("y1", line[1]);

				Map<String, Object> lineStyle = new HashMap<>();
				lineStyle.put("color", color);
				lineStyle.put("width", point_size);
				lineStyle.put("linecap", "round");

				shape.put("line", lineStyle);

				shapes.add(shape);
			}
		}
		return shapes;
	}

	public List<Object> getData() {
		List<Object> traces = new Vector<Object>();
		long offset = 0;
		int index = 0;
		for (int chr : bins.keySet()) {
			String color = colors[index % colors.length];
			if (!asShapes) {
				traces.add(getTrace(chr, offset, color));
			}
			traces.add(getTraceUnbinned(chr, offset, color));
			offset += (bins.get(chr).getLength() / bin_size) + gap;
			index++;
		}
		return traces;
	}

	private List<Object> getAnnotations() {
		List<Object> traces = new Vector<Object>();
		long offset = 0;
		for (int chr : bins.keySet()) {
			traces.addAll(getAnnotationsByChromosome(chr, offset));
			offset += (bins.get(chr).getLength() / bin_size) + gap;
		}
		return traces;
	}

	private Map<String, Object> getXAxis() {
		Map<String, Object> axis = new HashMap<>();
		axis.put("fixedrange", true);
		axis.put("tickmode", "array");
		List<Object> tickvals = new Vector<Object>();
		List<Object> ticktext = new Vector<Object>();
		long offset = 0;
		for (int chr : bins.keySet()) {
			ChromBin bin = bins.get(chr);
			ticktext.add(bin.chrom);
			tickvals.add(offset + bins.get(chr).getLength() / bin_size / 2);
			offset += (bins.get(chr).getLength() / bin_size) + gap;
		}
		axis.put("tickvals", tickvals);
		axis.put("ticktext", ticktext);
		axis.put("showgrid", false);
		return axis;
	}

	private Map<String, Object> getYAxis() {
		Map<String, Object> axis = new HashMap<>();
		axis.put("fixedrange", true);
		return axis;
	}

	private List<Object> getShapes() {
		List<Object> traces = new Vector<Object>();
		long offset = 0;
		int index = 0;
		for (int chr : bins.keySet()) {
			String color = colors[index % colors.length];
			traces.addAll(getShapes(chr, offset, color));
			offset += (bins.get(chr).getLength() / bin_size) + gap;
			index++;
		}
		System.out.println("Created " + traces.size() + " shapes");

		return traces;
	}

	private Map<String, Object> getTraceUnbinned(int chrIndex, long offset, String color) {
		Map<String, Object> trace = new HashMap<String, Object>();
		List<Double> x = new Vector<Double>();
		List<Double> y = new Vector<Double>();
		List<String> text = new Vector<String>();
		ChromBin bin = bins.get(chrIndex);
		for (Variant variant : bin.getUnbinnedVariants()) {
			x.add((variant.pos / bin_size) + offset);
			y.add(variant.pval);
			// TODO: add variant id
			text.add(variant.chrom + ":" + variant.pos + "<br>SNP: <br>Gene: ");

		}
		for (Variant variant : bin.getPeakVariants()) {
			x.add((variant.pos / bin_size) + offset);
			y.add(variant.pval);
			// TODO: add variant id
			text.add(variant.chrom + ":" + variant.pos + "<br>SNP: <br>Gene: ");

		}
		trace.put("x", x);
		trace.put("y", y);
		trace.put("text", text);
		trace.put("mode", "markers");
		trace.put("type", "scatter");
		trace.put("name", bin.chrom);
		trace.put("hovertemplate", "%{text}");

		Map<String, Object> marker = new HashMap<>();
		marker.put("size", point_size);

		marker.put("color", color);

		trace.put("marker", marker);
		return trace;
	}

	public List<Object> getAnnotationsByChromosome(int chrIndex, long offset) {
		List<Object> annotations = new Vector<Object>();
		ChromBin bin = bins.get(chrIndex);
		for (Variant variant : bin.getPeakVariants()) {
			Map<String, Object> annotation = new HashMap<>();
			annotation.put("x", (variant.pos / bin_size) + offset);
			annotation.put("y", variant.pval);
			annotation.put("xref", "x");
			annotation.put("yref", "y");
			annotation.put("text", variant.chrom + ":" + variant.pos);
			annotation.put("ax", 0);
			annotation.put("showarrow", true);
			annotation.put("arrowhead", 0);
			annotation.put("bordercolor", "#aaaaaa");
			annotation.put("ay", -40);
			annotations.add(annotation);
		}
		return annotations;
	}

	protected Map<String, Object> getLayout() {
		Map<String, Object> layout = new HashMap<String, Object>();
		layout.put("xaxis", getXAxis());
		layout.put("yaxis", getYAxis());
		layout.put("showlegend", false);
		layout.put("hovermode", "closest");
		if (annotations) {
			layout.put("annotations", getAnnotations());
		}
		layout.put("autosize", false);
		layout.put("width", 1300);
		layout.put("height", 800);
		if (asShapes) {
			layout.put("shapes", getShapes());
		}
		return layout;
	}

	public void setShowAnnotations(boolean annotations) {
		this.annotations = annotations;
	}

}

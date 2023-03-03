package com.github.lukfor.binner;

import java.util.HashMap;
import java.util.Map;

import com.github.lukfor.util.IPoppedCallback;
import com.github.lukfor.util.MaxPriorityQueue;

public class Binner {

	private Variant peak_best_variant;

	private ChrPosition peak_last_chrpos;

	private MaxPriorityQueue peak_pq;

	private MaxPriorityQueue unbinned_variant_pq;

	private double _qval_bin_size;

	private int num_significant_in_current_peak;

	private Map<Integer, ChromBin> bins;

	private double manhattan_peak_pval_threshold = 6;

	private double manhattan_peak_variant_counting_pval_threshold = -Math.log10(5e-8); // sign

	private int manhattan_peak_sprawl_dist = 200_000;

	private int manhattan_peak_max_count = 500;

	private int manhattan_num_unbinned = 500;

	private int bin_length = 3_000_000;

	public Binner() {

		peak_best_variant = null;
		peak_last_chrpos = null;
		peak_pq = new MaxPriorityQueue(500);
		unbinned_variant_pq = new MaxPriorityQueue(500);
		bins = new HashMap<Integer, ChromBin>();
		_qval_bin_size = 0.05;
		num_significant_in_current_peak = 0;
	}

	public void process_variant(Variant variant) {

		double qval = variant.pval;
		if (qval > 40) {
			_qval_bin_size = 0.2;
		} else if (qval > 20) {
			_qval_bin_size = 0.1;
		}

		if (variant.pval > manhattan_peak_pval_threshold) {
			if (peak_best_variant == null) {
				peak_best_variant = variant;
				peak_last_chrpos = new ChrPosition(variant.chrom, variant.pos);
				num_significant_in_current_peak = variant.pval > manhattan_peak_variant_counting_pval_threshold ? 1 : 0;
			} else if ((peak_last_chrpos.chr.equals(variant.chrom))
					&& (peak_last_chrpos.position + manhattan_peak_sprawl_dist > variant.pos)) {
				if (variant.pval > manhattan_peak_variant_counting_pval_threshold) {
					num_significant_in_current_peak += 1;
				}
				peak_last_chrpos = new ChrPosition(variant.chrom, variant.pos);
				if (variant.pval <= peak_best_variant.pval) {
					maybe_bin_variant(variant);
				} else {
					maybe_bin_variant(peak_best_variant);
					peak_best_variant = variant;
				}
			} else {
				peak_best_variant.num_significant_in_peak = num_significant_in_current_peak;
				num_significant_in_current_peak = variant.pval > manhattan_peak_variant_counting_pval_threshold ? 1 : 0;
				_maybe_peak_variant(peak_best_variant);
				peak_best_variant = variant;
				peak_last_chrpos = new ChrPosition(variant.chrom, variant.pos);
			}
		} else {
			maybe_bin_variant(variant);
		}
	}

	private void _maybe_peak_variant(Variant variant) {
		peak_pq.add_and_keep_size(variant, manhattan_peak_max_count, new IPoppedCallback() {
			@Override
			public void call(Variant variant) {
				maybe_bin_variant(variant);
			}
		});

	}

	private void maybe_bin_variant(Variant variant) {
		unbinned_variant_pq.add_and_keep_size(variant, manhattan_num_unbinned, new IPoppedCallback() {
			@Override
			public void call(Variant variant) {
				_bin_variant(variant);
			}
		});

	}

	private void _bin_variant(Variant variant) {
		int chrom_idx = Chromosome.getOrder(variant.chrom);
		ChromBin chrom_bin = bins.get(chrom_idx);
		if (chrom_bin == null) {
			chrom_bin = new ChromBin(chrom_idx + "");
			bins.put(chrom_idx, chrom_bin);
		}
		long pos_bin_id = Math.floorDiv(variant.pos, bin_length);

		Bin bin = chrom_bin.getBins().get(pos_bin_id);
		if (bin == null) {
			bin = new Bin();
			bin.chrom = variant.chrom;
			bin.startpos = pos_bin_id * bin_length;
			chrom_bin.addBin(pos_bin_id, bin);
		}
		bin.qval.add(_rounded(variant.pval));
	}

	private Double _rounded(double xval) {
		double x = ((int) Math.floor(xval / _qval_bin_size)) * _qval_bin_size + _qval_bin_size / 2;
		return x;
	}

	public MaxPriorityQueue getUnbinnedVariants() {
		return unbinned_variant_pq;
	}

	public MaxPriorityQueue getPeaks() {
		return peak_pq;
	}

	public Map<Integer, ChromBin> getBins() {
		return bins;
	}

}

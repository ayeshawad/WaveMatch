// src/cpp_binned_coverage.cpp

#include <Rcpp.h>
extern "C" {
#include "htslib/sam.h"
#include "htslib/hts.h"
}

// Helper to add one alignment interval into bin coverage
inline void add_alignment_to_bins(int32_t pos,
                                  int32_t aln_end,
                                  int step_bp,
                                  Rcpp::NumericVector &cov_vec) {
  if (pos < 0 || aln_end <= pos) {
    return;
  }
  int start_bin = pos / step_bp;
  int end_bin   = (aln_end - 1) / step_bp;
  int n_bins    = cov_vec.size();

  if (start_bin >= n_bins || end_bin < 0) {
    return;
  }
  if (start_bin < 0) start_bin = 0;
  if (end_bin >= n_bins) end_bin = n_bins - 1;

  for (int b = start_bin; b <= end_bin; ++b) {
    cov_vec[b] += 1.0;
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_binned_coverage_bam(const std::string &bam_path,
                                            const std::string &chrom,
                                            int chrom_length,
                                            int step_bp) {
  if (chrom_length <= 0) {
    Rcpp::stop("chrom_length must be positive");
  }
  if (step_bp <= 0) {
    Rcpp::stop("step_bp must be positive");
  }

  int n_bins = (chrom_length + step_bp - 1) / step_bp;
  Rcpp::NumericVector cov_vec(n_bins);
  std::fill(cov_vec.begin(), cov_vec.end(), 0.0);

  samFile *fp = sam_open(bam_path.c_str(), "r");
  if (fp == NULL) {
    Rcpp::stop("htslib sam_open() could not open file: " + bam_path);
  }

  bam_hdr_t *hdr = sam_hdr_read(fp);
  if (hdr == NULL) {
    sam_close(fp);
    Rcpp::stop("htslib sam_hdr_read() could not parse header: " + bam_path);
  }

  int tid = bam_name2id(hdr, chrom.c_str());
  if (tid < 0) {
    // chromosome not in this BAM: return zeros
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return cov_vec;
  }

  // Try to use the index if available to speed up
  hts_idx_t *idx = sam_index_load(fp, bam_path.c_str());

  bam1_t *b = bam_init1();
  if (b == NULL) {
    if (idx) hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    Rcpp::stop("htslib bam_init1() failed");
  }

  if (idx != NULL) {
    // Use iterator over region [0, chrom_length)
    hts_itr_t *iter = sam_itr_queryi(idx, tid, 0, chrom_length);
    if (iter == NULL) {
      bam_destroy1(b);
      hts_idx_destroy(idx);
      bam_hdr_destroy(hdr);
      sam_close(fp);
      Rcpp::stop("htslib sam_itr_queryi() failed for region");
    }

    while (sam_itr_next(fp, iter, b) >= 0) {
      // skip unmapped, secondary, supplementary
      uint16_t flag = b->core.flag;
      if (flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
        continue;
      }
      int32_t pos = b->core.pos;        // 0-based start
      int32_t end = bam_endpos(b);      // 1-based end
      add_alignment_to_bins(pos, end, step_bp, cov_vec);
    }

    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
  } else {
    // Fallback: scan whole BAM and filter by tid
    while (sam_read1(fp, hdr, b) >= 0) {
      uint16_t flag = b->core.flag;
      if (flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
        continue;
      }
      if (b->core.tid != tid) {
        continue;
      }
      int32_t pos = b->core.pos;
      int32_t end = bam_endpos(b);
      add_alignment_to_bins(pos, end, step_bp, cov_vec);
    }
  }

  bam_destroy1(b);
  bam_hdr_destroy(hdr);
  sam_close(fp);

  return cov_vec;
}

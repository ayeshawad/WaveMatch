#include <Rcpp.h>
extern "C" {
#include "htslib/sam.h"
}

// [[Rcpp::export]]
bool cpp_inferPairedEnd(const std::string& bam_path,
                          int maxReads = 10000) {
  uint16_t const PROPER_PAIR_FLAG = 0x2;
  uint16_t flag = 0;

  samFile* fp = sam_open(bam_path.c_str(), "r");
  if (fp == NULL) {
    Rcpp::stop("htslib's `sam_open()` could not open file: " + bam_path);
  }
  
  bam_hdr_t* hdr = sam_hdr_read(fp);
  if (hdr == NULL) {
    sam_close(fp);
    Rcpp::stop("htslib's `sam_hdr_read()` could not parse header: " + bam_path);
  }
  
  bam1_t* b = bam_init1();
  if (b == NULL) {
    /* 
     * Note: hdr was verified non-NULL above, ditto fp,
     * if we fail here, we can safely free and raise
    */
    bam_hdr_destroy(hdr);
    sam_close(fp);
    Rcpp::stop("htslib's `bam_init1()` failed");
  }
  
  bool isPairedEnd_ = false;
  int scannedRecordsCount = 0;
  
  /* now we check each alignment's sam flag for 0x2 (proper pair) */
  while (sam_read1(fp, hdr, b) >= 0) {
    scannedRecordsCount = scannedRecordsCount + 1;
    flag = b->core.flag;
    /* here */
    if (flag & PROPER_PAIR_FLAG) {
      isPairedEnd_ = true;
      break;
    }
    
    if (maxReads > 0 && scannedRecordsCount >= maxReads) {
      break;
    }
  }
  
  /* no leaks! */
  bam_destroy1(b);
  bam_hdr_destroy(hdr);
  sam_close(fp);
  
  return isPairedEnd_;
}

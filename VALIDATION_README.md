# ENCODE Peak Validation Script

## Overview

The `validate_encode_peaks.R` script tests wavematch accuracy against ENCODE gold-standard peaks for H3K27ac in K562 cells.

## Files Used

- **ENCODE Peaks**: `ENCFF000BWW.narrowPeak` (58,937 total peaks, 1,819 on chr22)
- **ChIP BAM**: `ENCFF000BWX.bam` (main ChIP sample matching the peaks)
- **Test Chromosome**: `chr22` (small, fast for testing)

## Usage

### Quick Start

```bash
cd /Users/away/Documents/GitHub/StructuredPeaks
Rscript validate_encode_peaks.R
```

The script will:
1. ✅ Auto-detect files in `/Users/away/Downloads/encode_h3k27ac_k562`
2. ✅ Load ENCODE peaks and filter to chr22
3. ✅ Load track data from BAM file
4. ✅ Run wavematch with optimized algorithm (Consenrich-aligned)
5. ✅ Compute overlap metrics (precision, recall, F1, Jaccard)
6. ✅ Save results to CSV

### Output

The script prints:
- **Peak counts**: ENCODE vs detected peaks
- **Metrics**: Precision, Recall, F1-score, Jaccard index
- **Performance**: Load time, peak calling time
- **Peak widths**: Statistics comparing detected vs ENCODE peaks
- **Success criteria**: Pass/fail for each metric threshold

Results are saved to:
```
/Users/away/Downloads/encode_h3k27ac_k562/wavematch_validation_chr22.csv
```

## Algorithm Details

The script uses the **optimized wavematch algorithm** that:
- ✅ Uses split-halves null sampling (Consenrich-aligned)
- ✅ Uses natural cascade template lengths
- ✅ Uses optimized C++ convolution
- ✅ Recenters peaks at point sources
- ✅ Applies FDR correction per template/cascade group

## Success Criteria

- **Precision ≥ 0.90**: 90% of detected peaks match ENCODE
- **Recall ≥ 0.90**: 90% of ENCODE peaks are detected
- **F1 ≥ 0.90**: Overall accuracy metric
- **Jaccard ≥ 0.70**: Overlap similarity metric

## Customization

To test a different chromosome, edit the script:
```r
test_chromosome <- "chr1"  # Change from "chr22"
```

To use a different BAM file, edit:
```r
chip_bam <- bam_files[grepl("ENCFF121RHF", bam_files)]  # Use different replicate
```

## Notes

- The script automatically handles BED 0-based to GRanges 1-based conversion
- Overlap matching requires ≥10% overlap between peaks
- Results are saved in CSV format for further analysis

# wavematch

Wavelet-based peak calling for ChIP-seq and CUT&RUN data.

## Overview

wavematch detects structured peaks in genomic signal tracks using wavelet-based template matching. The algorithm uses discrete wavelet templates (Haar, Daubechies) to identify peaks and applies statistical filtering with split-halves null sampling for robust significance testing.

## Quick Start

### From BAM file

```r
library(wavematch)

# Load track data from BAM
track_data <- wavematchFromBam(
  bamFile = "sample.bam",
  genome = "hg38",
  chrom = "chr22",
  binSize = 25L
)

# Call peaks
peaks <- wavematchFromTrack(
  chromosome = track_data$chromosome,
  intervals = track_data$intervals,
  values = track_data$values,
  template_names = "haar",
  cascade_levels = 3L,
  iters = 20000L,
  alpha = 0.05
)
```

### From existing track

```r
peaks <- wavematchFromTrack(
  chromosome = "chr22",
  intervals = seq(0, by = 25, length.out = 10000),
  values = signal_values,
  template_names = "haar",
  cascade_levels = 3L,
  iters = 20000L,
  alpha = 0.05
)
```

## Key Features

- **Wavelet templates**: Haar and Daubechies wavelets for flexible peak detection
- **Split-halves null sampling**: Robust statistical testing without control samples
- **Optimized C++ implementation**: Fast convolution and peak detection
- **Point source recentering**: Accurate peak summit identification
- **FDR correction**: Multiple testing correction per template/cascade group

## Parameters

- `template_names`: Wavelet templates to use ("haar", "db2", "db3")
- `cascade_levels`: Template resolution (higher = broader peaks)
- `iters`: Number of null samples for statistical testing
- `alpha`: FDR threshold for significance
- `mergePeaks`: Whether to merge overlapping peaks
- `mergeGapBP`: Maximum gap between peaks to merge

## Citation

If you use wavematch in your research, please cite:

[Citation information]

## License

[License information]

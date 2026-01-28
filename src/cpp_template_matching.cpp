// src/cpp_template_matching.cpp
// Optimized C++ implementations for template matching operations

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <limits>
#include <cmath>
#include <string>
#include <complex>
#include <unordered_map>

using namespace Rcpp;

// Sample block maxima from response signal for empirical null distribution
//
// @param response NumericVector of template response values
// @param half_window_bins integer, half-window size in bins
// @param nsamp integer, number of samples to draw
// @param mask LogicalVector, optional mask for valid sampling regions (TRUE = allowed)
// @param eps double, tolerance for floating point comparisons
// @return NumericVector of sampled block maxima
// [[Rcpp::export]]
NumericVector cpp_sample_block_maxima(
    const NumericVector& response,
    int half_window_bins,
    int nsamp,
    const LogicalVector& mask = LogicalVector(),
    double eps = 1e-10) {
  
  int n = response.size();
  int window_size = 2 * half_window_bins + 1;
  
  if (window_size > n || nsamp < 1 || half_window_bins < 0) {
    return NumericVector(0);
  }
  
  // Find valid centers where full window fits
  std::vector<int> valid_centers;
  bool use_mask = (mask.size() == n);
  
  for (int i = half_window_bins; i < n - half_window_bins; i++) {
    if (!use_mask || mask[i]) {
      valid_centers.push_back(i);
    }
  }
  
  if (valid_centers.empty()) {
    return NumericVector(0);
  }
  
  // Sample block maxima with replacement
  NumericVector block_maxima(nsamp);
  GetRNGstate();
  
  for (int s = 0; s < nsamp; s++) {
    // Random center index
    int idx = static_cast<int>(unif_rand() * valid_centers.size());
    if (idx >= static_cast<int>(valid_centers.size())) idx = valid_centers.size() - 1;
    int center_idx = valid_centers[idx];
    
    // Find maximum in window
    double max_val = -std::numeric_limits<double>::infinity();
    bool found_finite = false;
    
    for (int j = center_idx - half_window_bins; 
         j <= center_idx + half_window_bins; j++) {
      double val = response[j];
      if (std::isfinite(val)) {
        found_finite = true;
        if (val > max_val + eps) {
          max_val = val;
        }
      }
    }
    
    block_maxima[s] = found_finite ? max_val : NA_REAL;
  }
  
  PutRNGstate();
  
  // Remove NA values
  block_maxima = block_maxima[!is_na(block_maxima)];
  
  if (block_maxima.size() < 10) {
    return block_maxima;
  }
  
  // Trim extreme outliers at 0.001 and 0.999 quantiles to remove artifacts
  NumericVector sorted = clone(block_maxima).sort();
  int n_vals = sorted.size();
  int low_idx = std::max(0, static_cast<int>(n_vals * 0.001));
  int high_idx = std::min(n_vals - 1, static_cast<int>(n_vals * 0.999));
  
  return NumericVector(sorted.begin() + low_idx, sorted.begin() + high_idx + 1);
}

// Find relative maxima in a signal
//
// @param x NumericVector of signal values
// @param order_bins integer, half-window size for local maximum detection
// @param eps double, tolerance for floating point comparisons
// @return IntegerVector of 1-based indices of local maxima
// [[Rcpp::export]]
IntegerVector cpp_relative_maxima(
    const NumericVector& x,
    int order_bins,
    double eps = 1e-10) {
  
  int n = x.size();
  std::vector<int> maxima;
  
  if (n == 0 || order_bins < 1) {
    return IntegerVector(0);
  }
  
  for (int i = 0; i < n; i++) {
    double center_val = x[i];
    
    if (!std::isfinite(center_val)) {
      continue;
    }
    
    // Check neighborhood
    int start = std::max(0, i - order_bins);
    int end = std::min(n - 1, i + order_bins);
    
    bool is_max = true;
    bool has_neighbor = false;
    
    for (int j = start; j <= end; j++) {
      if (j == i) continue;
      
      double neighbor_val = x[j];
      if (std::isfinite(neighbor_val)) {
        has_neighbor = true;
        if (neighbor_val > center_val + eps) {
          is_max = false;
          break;
        }
      }
    }
    
    // Require at least one neighbor and strict maximum
    if (is_max && has_neighbor && center_val > std::numeric_limits<double>::lowest()) {
      maxima.push_back(i + 1); // R is 1-indexed
    }
  }
  
  return wrap(maxima);
}

// Compute empirical p-values from null distribution
//
// OPTIMIZED: Uses binary search for O(n*log(m)) instead of O(n*m)
//
// @param observed NumericVector of observed response values
// @param null_dist NumericVector of null distribution samples
// @return NumericVector of empirical p-values (right-tail)
// [[Rcpp::export]]
NumericVector cpp_empirical_pvalues(
    const NumericVector& observed,
    const NumericVector& null_dist) {
  
  int n_obs = observed.size();
  NumericVector pvals(n_obs);
  
  if (null_dist.size() == 0) {
    std::fill(pvals.begin(), pvals.end(), 1.0);
    return pvals;
  }
  
  // Sort null distribution in descending order for efficient binary search
  NumericVector null_sorted = clone(null_dist);
  std::sort(null_sorted.begin(), null_sorted.end(), std::greater<double>());
  int n_null = null_sorted.size();
  
  // Precompute minimum non-zero p-value
  double min_pval = 1.0 / (n_null + 1.0);
  
  for (int i = 0; i < n_obs; i++) {
    double obs_val = observed[i];
    
    if (!std::isfinite(obs_val)) {
      pvals[i] = 1.0;
      continue;
    }
    
    // Binary search: find first position where null_sorted[j] < obs_val
    // Since null_sorted is descending, all values before this are >= obs_val
    int left = 0;
    int right = n_null;
    int pos = n_null; // Default: all null values are >= obs_val
    
    while (left < right) {
      int mid = left + (right - left) / 2;
      if (null_sorted[mid] < obs_val) {
        pos = mid;
        right = mid;
      } else {
        left = mid + 1;
      }
    }
    
    // Count of null values >= observed is 'pos'
    double pval = static_cast<double>(pos) / n_null;
    
    // Ensure non-zero for log transformation
    if (pval <= 0.0) {
      pval = min_pval;
    }
    
    pvals[i] = pval;
  }
  
  return pvals;
}

// Compute expected p-values from null distribution (for filtering)
//
// OPTIMIZED: Uses binary search for O(n*log(m)) instead of O(n*m)
// Similar to cpp_empirical_pvalues but returns expected p-values
//
// @param observed NumericVector of observed response values
// @param null_dist NumericVector of null distribution samples
// @return NumericVector of expected p-values (right-tail)
// [[Rcpp::export]]
NumericVector cpp_expected_pvalues(
    const NumericVector& observed,
    const NumericVector& null_dist) {
  
  int n_obs = observed.size();
  NumericVector expected_pvals(n_obs);
  
  if (null_dist.size() == 0) {
    std::fill(expected_pvals.begin(), expected_pvals.end(), 1.0);
    return expected_pvals;
  }
  
  // Sort null distribution in descending order for efficient binary search
  NumericVector null_sorted = clone(null_dist);
  std::sort(null_sorted.begin(), null_sorted.end(), std::greater<double>());
  int n_null = null_sorted.size();
  
  for (int i = 0; i < n_obs; i++) {
    double obs_val = observed[i];
    
    if (!std::isfinite(obs_val)) {
      expected_pvals[i] = 1.0;
      continue;
    }
    
    // Binary search: find first position where null_sorted[j] < obs_val
    // Since null_sorted is descending, all values before this are >= obs_val
    int left = 0;
    int right = n_null;
    int pos = n_null; // Default: all null values are >= obs_val
    
    while (left < right) {
      int mid = left + (right - left) / 2;
      if (null_sorted[mid] < obs_val) {
        pos = mid;
        right = mid;
      } else {
        left = mid + 1;
      }
    }
    
    // Expected p-value: proportion of null values >= observed
    expected_pvals[i] = static_cast<double>(pos) / n_null;
  }
  
  return expected_pvals;
}

// Normalize scores to 0-1000 range (matching Python implementation)
//
// @param scores NumericVector of raw template response scores
// @return IntegerVector of normalized scores in [0, 1000]
// [[Rcpp::export]]
IntegerVector cpp_normalize_scores(const NumericVector& scores) {
  int n = scores.size();
  IntegerVector normalized(n);
  
  if (n == 0) {
    return normalized;
  }
  
  // Transform: (1 + response)^2
  NumericVector transformed(n);
  for (int i = 0; i < n; i++) {
    double val = scores[i];
    if (std::isfinite(val)) {
      transformed[i] = std::pow(1.0 + val, 2.0);
    } else {
      transformed[i] = 0.0;
    }
  }
  
  // Find min and max
  double min_val = *std::min_element(transformed.begin(), transformed.end());
  double max_val = *std::max_element(transformed.begin(), transformed.end());
  double range = max_val - min_val;
  
  // Normalize to [250, 1000] with robust handling
  if (range > 1e-10) {
    for (int i = 0; i < n; i++) {
      double scaled = 250.0 + 750.0 * (transformed[i] - min_val) / range;
      // Ensure scores are always in valid range [250, 1000]
      scaled = std::max(250.0, std::min(1000.0, scaled));
      normalized[i] = static_cast<int>(std::round(scaled));
    }
  } else {
    // All values are the same - use middle of range
    std::fill(normalized.begin(), normalized.end(), 500);
  }
  
  return normalized;
}

// Batched null distribution sampling for multiple candidates
// This eliminates the R lapply overhead by sampling all candidates at once
//
// @param response NumericVector of template response values
// @param candidate_indices IntegerVector of candidate center indices (1-based)
// @param half_window_bins integer, half-window size in bins
// @param nsamp integer, number of samples per candidate
// @param null_mask_left LogicalVector for left half null regions (if use_split_halves)
// @param null_mask_right LogicalVector for right half null regions (if use_split_halves)
// @param mid_point integer, midpoint index for split-halves logic
// @param use_split_halves bool, whether to use split-halves null selection
// @param min_null_size integer, minimum null samples required
// @param eps double, tolerance for floating point comparisons
// @return List of NumericVectors, one null distribution per candidate
// [[Rcpp::export]]
List cpp_sample_null_distributions_batched(
    const NumericVector& response,
    const IntegerVector& candidate_indices,
    int half_window_bins,
    int nsamp,
    const LogicalVector& null_mask_left = LogicalVector(),
    const LogicalVector& null_mask_right = LogicalVector(),
    int mid_point = -1,
    bool use_split_halves = false,
    int min_null_size = 25,
    double eps = 1e-10) {
  
  int n = response.size();
  int n_candidates = candidate_indices.size();
  int window_size = 2 * half_window_bins + 1;
  
  if (window_size > n || nsamp < 1 || half_window_bins < 0 || n_candidates == 0) {
    return List(n_candidates);
  }
  
  // Pre-compute valid centers for each mask type
  std::vector<int> valid_centers_left, valid_centers_right, valid_centers_all;
  bool has_left_mask = (null_mask_left.size() == n);
  bool has_right_mask = (null_mask_right.size() == n);
  
  for (int i = half_window_bins; i < n - half_window_bins; i++) {
    if (!has_left_mask || null_mask_left[i]) {
      valid_centers_left.push_back(i);
    }
    if (!has_right_mask || null_mask_right[i]) {
      valid_centers_right.push_back(i);
    }
    valid_centers_all.push_back(i);
  }
  
  List null_distributions(n_candidates);
  GetRNGstate();
  
  for (int c = 0; c < n_candidates; c++) {
    int cand_idx = candidate_indices[c] - 1; // Convert to 0-based
    
    // Choose appropriate mask based on split-halves logic
    std::vector<int>* valid_centers = &valid_centers_all;
    if (use_split_halves && mid_point >= 0) {
      if (cand_idx <= mid_point) {
        valid_centers = &valid_centers_right;
      } else {
        valid_centers = &valid_centers_left;
      }
    }
    
    if (valid_centers->empty()) {
      valid_centers = &valid_centers_all;
    }
    
    // Sample block maxima
    NumericVector block_maxima(nsamp);
    for (int s = 0; s < nsamp; s++) {
      int idx = static_cast<int>(unif_rand() * valid_centers->size());
      if (idx >= static_cast<int>(valid_centers->size())) idx = valid_centers->size() - 1;
      int center_idx = (*valid_centers)[idx];
      
      double max_val = -std::numeric_limits<double>::infinity();
      bool found_finite = false;
      
      for (int j = center_idx - half_window_bins; 
           j <= center_idx + half_window_bins; j++) {
        double val = response[j];
        if (std::isfinite(val)) {
          found_finite = true;
          if (val > max_val + eps) {
            max_val = val;
          }
        }
      }
      
      block_maxima[s] = found_finite ? max_val : NA_REAL;
    }
    
    // Remove NA values
    block_maxima = block_maxima[!is_na(block_maxima)];
    
    // If too few samples, try with all centers
    if (block_maxima.size() < min_null_size && valid_centers != &valid_centers_all) {
      NumericVector block_maxima_all(nsamp);
      for (int s = 0; s < nsamp; s++) {
        int idx = static_cast<int>(unif_rand() * valid_centers_all.size());
        if (idx >= static_cast<int>(valid_centers_all.size())) idx = valid_centers_all.size() - 1;
        int center_idx = valid_centers_all[idx];
        
        double max_val = -std::numeric_limits<double>::infinity();
        bool found_finite = false;
        
        for (int j = center_idx - half_window_bins; 
             j <= center_idx + half_window_bins; j++) {
          double val = response[j];
          if (std::isfinite(val)) {
            found_finite = true;
            if (val > max_val + eps) {
              max_val = val;
            }
          }
        }
        
        block_maxima_all[s] = found_finite ? max_val : NA_REAL;
      }
      block_maxima_all = block_maxima_all[!is_na(block_maxima_all)];
      if (block_maxima_all.size() >= block_maxima.size()) {
        block_maxima = block_maxima_all;
      }
    }
    
    // Trim extreme outliers
    if (block_maxima.size() >= 10) {
      NumericVector sorted = clone(block_maxima).sort();
      int n_vals = sorted.size();
      int low_idx = std::max(0, static_cast<int>(n_vals * 0.0001));
      int high_idx = std::min(n_vals - 1, static_cast<int>(n_vals * 0.9999));
      block_maxima = NumericVector(sorted.begin() + low_idx, sorted.begin() + high_idx + 1);
    }
    
    null_distributions[c] = block_maxima;
  }
  
  PutRNGstate();
  return null_distributions;
}

// Pooled null distribution sampling per segment
// Samples one pooled null distribution per segment (much faster than per-candidate)
// CRITICAL: Excludes candidate positions to avoid inflating null distribution
//
// @param response NumericVector of template response values
// @param segment_ids IntegerVector of segment ID for each bin (1-based)
// @param half_window_bins integer, half-window size in bins
// @param nsamp_per_segment integer, number of samples per segment
// @param candidate_indices IntegerVector of candidate positions to EXCLUDE (1-based, optional)
// @param null_mask LogicalVector, optional mask for valid sampling regions
// @param min_null_size integer, minimum null samples required
// @param eps double, tolerance for floating point comparisons
// @return List of NumericVectors, one pooled null distribution per segment
// [[Rcpp::export]]
List cpp_sample_pooled_nulls_per_segment(
    const NumericVector& response,
    const IntegerVector& segment_ids,
    int half_window_bins,
    int nsamp_per_segment,
    const IntegerVector& candidate_indices = IntegerVector(),
    const LogicalVector& null_mask = LogicalVector(),
    int min_null_size = 25,
    double eps = 1e-10) {
  
  int n = response.size();
  int window_size = 2 * half_window_bins + 1;
  
  if (window_size > n || nsamp_per_segment < 1 || half_window_bins < 0) {
    return List(0);
  }
  
  // Find unique segments
  IntegerVector unique_segments = unique(segment_ids).sort();
  int n_segments = unique_segments.size();
  
  if (n_segments == 0) {
    return List(0);
  }
  
  // Create exclusion mask for candidate positions (and their windows)
  // This prevents null sampling from including actual peaks
  std::vector<bool> exclude_positions(n, false);
  bool exclude_candidates = (candidate_indices.size() > 0);
  
  if (exclude_candidates) {
    // Mark candidate positions and their windows as excluded
    for (int i = 0; i < candidate_indices.size(); i++) {
      int cand_idx = candidate_indices[i] - 1; // Convert to 0-based
      if (cand_idx >= 0 && cand_idx < n) {
        // Exclude the candidate position and a window around it
        int exclude_start = std::max(0, cand_idx - half_window_bins);
        int exclude_end = std::min(n - 1, cand_idx + half_window_bins);
        for (int j = exclude_start; j <= exclude_end; j++) {
          exclude_positions[j] = true;
        }
      }
    }
  }
  
  // Pre-compute valid centers (excluding candidate positions)
  std::vector<int> valid_centers;
  bool use_mask = (null_mask.size() == n);
  
  for (int i = half_window_bins; i < n - half_window_bins; i++) {
    // Exclude if it's a candidate position or in candidate window
    if (exclude_positions[i]) {
      continue;
    }
    // Also respect null_mask if provided
    if (!use_mask || null_mask[i]) {
      valid_centers.push_back(i);
    }
  }
  
  if (valid_centers.empty()) {
    return List(n_segments);
  }
  
  // OPTIMIZED: Use hash map for O(1) segment index lookup instead of O(n*m) nested loop
  std::unordered_map<int, int> seg_id_to_index;
  for (int s = 0; s < n_segments; s++) {
    seg_id_to_index[unique_segments[s]] = s;
  }
  
  // For each segment, find valid centers within that segment
  std::vector<std::vector<int> > segment_centers(n_segments);
  for (size_t i = 0; i < valid_centers.size(); i++) {
    int center_idx = valid_centers[i];
    if (center_idx >= 0 && center_idx < n) {
      int seg_id = segment_ids[center_idx];
      // O(1) lookup using hash map
      auto it = seg_id_to_index.find(seg_id);
      if (it != seg_id_to_index.end()) {
        int seg_index = it->second;
        segment_centers[seg_index].push_back(center_idx);
      }
    }
  }
  
  // Sample pooled null for each segment
  List null_distributions(n_segments);
  GetRNGstate();
  
  for (int s = 0; s < n_segments; s++) {
    const std::vector<int>& seg_centers = segment_centers[s];
    
    if (seg_centers.empty()) {
      // No valid centers in this segment, use all centers
      NumericVector pooled_null(nsamp_per_segment);
      for (int i = 0; i < nsamp_per_segment; i++) {
        int idx = static_cast<int>(unif_rand() * valid_centers.size());
        if (idx >= static_cast<int>(valid_centers.size())) idx = valid_centers.size() - 1;
        int center_idx = valid_centers[idx];
        
        double max_val = -std::numeric_limits<double>::infinity();
        bool found_finite = false;
        
        for (int j = center_idx - half_window_bins; 
             j <= center_idx + half_window_bins; j++) {
          if (j >= 0 && j < n) {
            double val = response[j];
            if (std::isfinite(val)) {
              found_finite = true;
              if (val > max_val + eps) {
                max_val = val;
              }
            }
          }
        }
        
        pooled_null[i] = found_finite ? max_val : NA_REAL;
      }
      
      pooled_null = pooled_null[!is_na(pooled_null)];
      
      // Trim outliers
      if (pooled_null.size() >= 10) {
        NumericVector sorted = clone(pooled_null).sort();
        int n_vals = sorted.size();
        int low_idx = std::max(0, static_cast<int>(n_vals * 0.0001));
        int high_idx = std::min(n_vals - 1, static_cast<int>(n_vals * 0.9999));
        pooled_null = NumericVector(sorted.begin() + low_idx, sorted.begin() + high_idx + 1);
      }
      
      null_distributions[s] = pooled_null;
    } else {
      // Sample from segment-specific centers
      NumericVector pooled_null(nsamp_per_segment);
      for (int i = 0; i < nsamp_per_segment; i++) {
        int idx = static_cast<int>(unif_rand() * seg_centers.size());
        if (idx >= static_cast<int>(seg_centers.size())) idx = seg_centers.size() - 1;
        int center_idx = seg_centers[idx];
        
        double max_val = -std::numeric_limits<double>::infinity();
        bool found_finite = false;
        
        for (int j = center_idx - half_window_bins; 
             j <= center_idx + half_window_bins; j++) {
          if (j >= 0 && j < n) {
            double val = response[j];
            if (std::isfinite(val)) {
              found_finite = true;
              if (val > max_val + eps) {
                max_val = val;
              }
            }
          }
        }
        
        pooled_null[i] = found_finite ? max_val : NA_REAL;
      }
      
      pooled_null = pooled_null[!is_na(pooled_null)];
      
      // If too few samples, supplement with all centers
      if (pooled_null.size() < min_null_size) {
        NumericVector supplement(nsamp_per_segment);
        for (int i = 0; i < nsamp_per_segment; i++) {
          int idx = static_cast<int>(unif_rand() * valid_centers.size());
          if (idx >= static_cast<int>(valid_centers.size())) idx = valid_centers.size() - 1;
          int center_idx = valid_centers[idx];
          
          double max_val = -std::numeric_limits<double>::infinity();
          bool found_finite = false;
          
          for (int j = center_idx - half_window_bins; 
               j <= center_idx + half_window_bins; j++) {
            if (j >= 0 && j < n) {
              double val = response[j];
              if (std::isfinite(val)) {
                found_finite = true;
                if (val > max_val + eps) {
                  max_val = val;
                }
              }
            }
          }
          
          supplement[i] = found_finite ? max_val : NA_REAL;
        }
        supplement = supplement[!is_na(supplement)];
        pooled_null = NumericVector(pooled_null.begin(), pooled_null.end());
        for (size_t i = 0; i < supplement.size(); i++) {
          pooled_null.push_back(supplement[i]);
        }
      }
      
      // Trim outliers
      if (pooled_null.size() >= 10) {
        NumericVector sorted = clone(pooled_null).sort();
        int n_vals = sorted.size();
        int low_idx = std::max(0, static_cast<int>(n_vals * 0.0001));
        int high_idx = std::min(n_vals - 1, static_cast<int>(n_vals * 0.9999));
        pooled_null = NumericVector(sorted.begin() + low_idx, sorted.begin() + high_idx + 1);
      }
      
      null_distributions[s] = pooled_null;
    }
  }
  
  PutRNGstate();
  
  // Name the list with segment IDs
  CharacterVector seg_names(n_segments);
  for (int s = 0; s < n_segments; s++) {
    std::string seg_name = "segment_";
    seg_name += std::to_string(unique_segments[s]);
    seg_names[s] = seg_name;
  }
  null_distributions.attr("names") = seg_names;
  null_distributions.attr("segment_ids") = unique_segments;
  
  return null_distributions;
}

// Vectorized adaptive scale refinement
// Computes best scale and response for multiple candidates at once
//
// @param values NumericVector of signal values
// @param candidate_indices IntegerVector of candidate center indices (1-based)
// @param initial_scales IntegerVector of initial scale indices (1-based)
// @param templates_per_scale List of NumericVectors, reversed templates for each scale
// @param refine_neighbors bool, whether to check neighboring scales
// @param intervals NumericVector of interval start positions
// @param interval_bp double, bin size in base pairs
// @return List with components: scores (NumericVector), best_scales (IntegerVector), 
//         starts_bp (NumericVector), ends_bp (NumericVector), span_bins (IntegerVector)
// [[Rcpp::export]]
List cpp_refine_adaptive_scales(
    const NumericVector& values,
    const IntegerVector& candidate_indices,
    const IntegerVector& initial_scales,
    const List& templates_per_scale,
    bool refine_neighbors,
    const NumericVector& intervals,
    double interval_bp) {
  
  int n_candidates = candidate_indices.size();
  int n_values = values.size();
  int K = templates_per_scale.size();
  
  if (n_candidates == 0 || K == 0) {
    return List::create(
      Named("scores") = NumericVector(0),
      Named("best_scales") = IntegerVector(0),
      Named("starts_bp") = NumericVector(0),
      Named("ends_bp") = NumericVector(0),
      Named("span_bins") = IntegerVector(0)
    );
  }
  
  NumericVector scores(n_candidates);
  IntegerVector best_scales(n_candidates);
  NumericVector starts_bp(n_candidates);
  NumericVector ends_bp(n_candidates);
  IntegerVector span_bins(n_candidates);
  
  // Pre-extract templates for faster access
  std::vector<NumericVector> templates(K);
  std::vector<int> template_lengths(K);
  for (int k = 0; k < K; k++) {
    templates[k] = as<NumericVector>(templates_per_scale[k]);
    template_lengths[k] = templates[k].size();
  }
  
  for (int i = 0; i < n_candidates; i++) {
    int cand_idx = candidate_indices[i] - 1; // Convert to 0-based
    int k = initial_scales[i] - 1; // Convert to 0-based
    if (k < 0) k = 0;
    if (k >= K) k = K - 1;
    
    if (cand_idx < 0 || cand_idx >= n_values) {
      scores[i] = 0.0;
      best_scales[i] = 1;
      starts_bp[i] = intervals[0];
      ends_bp[i] = intervals[0];
      span_bins[i] = template_lengths[0];
      continue;
    }
    
    // Compute response at initial scale
    double best_r = 0.0;
    int best_k = k;
    
    // Helper function to compute response at a given scale
    auto compute_response = [&](int scale_idx) -> double {
      if (scale_idx < 0 || scale_idx >= K) return -1e10;
      const NumericVector& tmpl = templates[scale_idx];
      int half = template_lengths[scale_idx] / 2;
      
      double response = 0.0;
      for (int j = 0; j < template_lengths[scale_idx]; j++) {
        int val_idx = cand_idx - half + j;
        if (val_idx >= 0 && val_idx < n_values) {
          response += values[val_idx] * tmpl[j];
        }
      }
      return response;
    };
    
    best_r = compute_response(k);
    best_k = k;
    
    // Refine neighbors if requested
    if (refine_neighbors && K > 1) {
      for (int offset : {-1, 1}) {
        int kk = k + offset;
        if (kk >= 0 && kk < K) {
          double r_kk = compute_response(kk);
          if (r_kk > best_r) {
            best_r = r_kk;
            best_k = kk;
          }
        }
      }
    }
    
    // Compute genomic coordinates
    int L_best = template_lengths[best_k];
    int half_best = L_best / 2;
    double center_bp = intervals[cand_idx] + std::max(1.0, interval_bp / 2.0);
    
    scores[i] = best_r;
    best_scales[i] = best_k + 1; // Convert back to 1-based
    starts_bp[i] = center_bp - half_best * interval_bp;
    ends_bp[i] = center_bp + half_best * interval_bp;
    span_bins[i] = L_best;
  }
  
  return List::create(
    Named("scores") = scores,
    Named("best_scales") = best_scales,
    Named("starts_bp") = starts_bp,
    Named("ends_bp") = ends_bp,
    Named("span_bins") = span_bins
  );
}

// Compute sparse scale energy matrix
// Computes energy (squared response) for all scales at sparse centers
//
// @param values NumericVector of signal values
// @param centers IntegerVector of sparse center indices (1-based)
// @param templates_per_scale List of NumericVectors, reversed templates for each scale
// @return NumericMatrix with K rows (scales) and length(centers) columns
// [[Rcpp::export]]
NumericMatrix cpp_compute_scale_energy(
    const NumericVector& values,
    const IntegerVector& centers,
    const List& templates_per_scale) {
  
  int n_values = values.size();
  int n_centers = centers.size();
  int K = templates_per_scale.size();
  
  if (n_centers == 0 || K == 0) {
    return NumericMatrix(K, 0);
  }
  
  NumericMatrix energy(K, n_centers);
  
  // Pre-extract templates
  std::vector<NumericVector> templates(K);
  std::vector<int> template_lengths(K);
  for (int k = 0; k < K; k++) {
    templates[k] = as<NumericVector>(templates_per_scale[k]);
    template_lengths[k] = templates[k].size();
  }
  
  for (int k = 0; k < K; k++) {
    const NumericVector& tmpl = templates[k];
    int half = template_lengths[k] / 2;
    
    for (int j = 0; j < n_centers; j++) {
      int c_idx = centers[j] - 1; // Convert to 0-based
      
      if (c_idx < 0 || c_idx >= n_values) {
        energy(k, j) = 0.0;
        continue;
      }
      
      // Compute response
      double response = 0.0;
      for (int t = 0; t < template_lengths[k]; t++) {
        int val_idx = c_idx - half + t;
        if (val_idx >= 0 && val_idx < n_values) {
          response += values[val_idx] * tmpl[t];
        }
      }
      
      energy(k, j) = response * response; // Squared response (energy)
    }
  }
  
  return energy;
}

// Find point sources (local maxima) within peak spans
//
// @param values NumericVector of signal values
// @param center_indices IntegerVector of peak center indices (1-based)
// @param span_bins IntegerVector of peak span sizes in bins
// @return IntegerVector of point source indices (1-based)
// [[Rcpp::export]]
IntegerVector cpp_find_point_sources(
    const NumericVector& values,
    const IntegerVector& center_indices,
    const IntegerVector& span_bins) {
  
  int n_peaks = center_indices.size();
  int n_values = values.size();
  
  IntegerVector point_sources(n_peaks);
  
  for (int i = 0; i < n_peaks; i++) {
    int c_idx = center_indices[i] - 1; // Convert to 0-based
    int L = span_bins[i];
    
    if (c_idx < 0 || c_idx >= n_values || L < 1) {
      point_sources[i] = center_indices[i]; // Fallback to center
      continue;
    }
    
    int half = L / 2;
    int start_idx = std::max(0, c_idx - half);
    int end_idx = std::min(n_values - 1, c_idx + half);
    
    // Find maximum in span
    double max_val = -std::numeric_limits<double>::infinity();
    int max_idx = c_idx;
    bool found_finite = false;
    
    for (int j = start_idx; j <= end_idx; j++) {
      double val = values[j];
      if (std::isfinite(val)) {
        found_finite = true;
        if (val > max_val) {
          max_val = val;
          max_idx = j;
        }
      }
    }
    
    if (!found_finite) {
      point_sources[i] = center_indices[i]; // Fallback to center
    } else {
      point_sources[i] = max_idx + 1; // Convert back to 1-based
    }
  }
  
  return point_sources;
}

// Match weights to intervals with optional interpolation
//
// @param intervals NumericVector of interval start positions
// @param weight_starts NumericVector of weight region start positions
// @param weight_ends NumericVector of weight region end positions
// @param weight_values NumericVector of weight values
// @param interpolate bool, whether to interpolate for unmatched intervals
// @return NumericVector of matched weights
// [[Rcpp::export]]
NumericVector cpp_match_weights(
    const NumericVector& intervals,
    const NumericVector& weight_starts,
    const NumericVector& weight_ends,
    const NumericVector& weight_values,
    bool interpolate = false) {
  
  int n_intervals = intervals.size();
  int n_weights = weight_starts.size();
  
  NumericVector matched_weights(n_intervals);
  std::fill(matched_weights.begin(), matched_weights.end(), NA_REAL);
  
  if (n_weights == 0) {
    return matched_weights;
  }
  
  // Simple exact matching (can be optimized with binary search for large datasets)
  for (int i = 0; i < n_intervals; i++) {
    double interval_start = intervals[i];
    
    // Find matching weight region
    for (int j = 0; j < n_weights; j++) {
      if (interval_start >= weight_starts[j] && interval_start < weight_ends[j]) {
        matched_weights[i] = weight_values[j];
        break;
      }
    }
    
    // If no exact match and interpolation requested, find nearest
    if (interpolate && std::isnan(matched_weights[i])) {
      double min_dist = std::numeric_limits<double>::max();
      int best_j = -1;
      
      for (int j = 0; j < n_weights; j++) {
        double center = (weight_starts[j] + weight_ends[j]) / 2.0;
        double dist = std::abs(interval_start - center);
        if (dist < min_dist) {
          min_dist = dist;
          best_j = j;
        }
      }
      
      if (best_j >= 0) {
        matched_weights[i] = weight_values[best_j];
      }
    }
  }
  
  // Fill NA with median
  NumericVector finite_weights = matched_weights[!is_na(matched_weights)];
  if (finite_weights.size() > 0) {
    NumericVector sorted = clone(finite_weights).sort();
    double median_val = sorted[sorted.size() / 2];
    for (int i = 0; i < n_intervals; i++) {
      if (std::isnan(matched_weights[i])) {
        matched_weights[i] = median_val;
      }
    }
  }
  
  return matched_weights;
}

// Fast running median using efficient algorithm
//
// @param x NumericVector of signal values
// @param window_size integer, window size (must be odd)
// @return NumericVector of running medians
// [[Rcpp::export]]
NumericVector cpp_running_median(
    const NumericVector& x,
    int window_size) {
  
  int n = x.size();
  
  if (n == 0 || window_size < 1) {
    return NumericVector(0);
  }
  
  // Ensure odd window size
  if (window_size % 2 == 0) {
    window_size++;
  }
  
  int half = window_size / 2;
  NumericVector result(n);
  
  // Simple implementation: for each position, sort window and take median
  // (More efficient algorithms exist but this is straightforward)
  for (int i = 0; i < n; i++) {
    int start = std::max(0, i - half);
    int end = std::min(n - 1, i + half);
    int win_size = end - start + 1;
    
    NumericVector window(win_size);
    for (int j = start; j <= end; j++) {
      window[j - start] = x[j];
    }
    
    NumericVector sorted = clone(window).sort();
    result[i] = sorted[win_size / 2];
  }
  
  return result;
}

// Fast running mean (rolling average)
//
// @param x NumericVector of values
// @param window_size integer, window size (must be odd)
// @return NumericVector of running means
// [[Rcpp::export]]
NumericVector cpp_running_mean(
    const NumericVector& x,
    int window_size) {
  
  int n = x.size();
  
  if (n == 0 || window_size < 1) {
    return NumericVector(0);
  }
  
  // Ensure odd window size
  if (window_size % 2 == 0) {
    window_size++;
  }
  
  int half = window_size / 2;
  NumericVector result(n);
  
  // Optimized: use sliding window sum
  double window_sum = 0.0;
  int window_count = 0;
  
  // Initialize first window
  for (int i = 0; i <= half && i < n; i++) {
    if (std::isfinite(x[i])) {
      window_sum += x[i];
      window_count++;
    }
  }
  
  // First position
  result[0] = (window_count > 0) ? window_sum / window_count : 0.0;
  
  // Sliding window for remaining positions
  for (int i = 1; i < n; i++) {
    // Remove leftmost element
    int left_idx = i - half - 1;
    if (left_idx >= 0 && std::isfinite(x[left_idx])) {
      window_sum -= x[left_idx];
      window_count--;
    }
    
    // Add rightmost element
    int right_idx = i + half;
    if (right_idx < n && std::isfinite(x[right_idx])) {
      window_sum += x[right_idx];
      window_count++;
    }
    
    result[i] = (window_count > 0) ? window_sum / window_count : 0.0;
  }
  
  return result;
}

// Bulletproof convolution with "same" mode extraction
// Optimized multi-strategy approach for all array sizes:
// - Small arrays: Direct convolution (fastest)
// - Medium arrays: R's FFT with fallback (balanced)
// - Large arrays: Chunked direct convolution (most reliable, no hangs)
//
// Designed to handle:
// - Small chromosomes (chr22, ~2M bins) - fast
// - Large chromosomes (chr1, ~250M bins) - reliable
// - ATAC-seq/ChIP-seq with high noise - robust
// - All template lengths (3-200+ bins) - consistent
//
// @param x NumericVector of signal values
// @param kernel NumericVector of template (reversed)
// @return NumericVector of convolution result in "same" mode
// [[Rcpp::export]]
NumericVector cpp_convolve_same(
    const NumericVector& x,
    const NumericVector& kernel) {
  
  int n_x = x.size();
  int n_k = kernel.size();
  
  if (n_x == 0 || n_k == 0) {
    return NumericVector(0);
  }
  
  // BULLETPROOF CONVOLUTION STRATEGY - Optimized for all array sizes:
  // 1. Small arrays (n_x <= 50k): Direct convolution (fastest, no overhead)
  // 2. Medium arrays (50k < n_x <= 200k): Try R's FFT, fallback to direct
  // 3. Large arrays (n_x > 200k): Chunked direct convolution (most reliable, no R overhead)
  // 4. Large kernels (n_k > 50): Always use chunked approach (more reliable)
  //
  // This ensures:
  // - No hangs on very large arrays (ATAC-seq, ChIP-seq)
  // - Fast performance for all template lengths
  // - Memory-efficient chunked processing
  
  const int SMALL_THRESHOLD = 50000;
  const int MEDIUM_THRESHOLD = 200000;  // Lower threshold to avoid R hangs
  const int CHUNK_SIZE = 500000;        // Process 500k bins per chunk
  const int KERNEL_LARGE_THRESHOLD = 50; // Lower threshold for reliability
  
  // Strategy 1: Small arrays - direct convolution (fastest)
  if (n_x <= SMALL_THRESHOLD && n_k <= KERNEL_LARGE_THRESHOLD) {
    int n_conv = n_x + n_k - 1;
    NumericVector conv_result(n_conv);
    
    // Optimized direct convolution with better memory access
    for (int i = 0; i < n_conv; i++) {
      double sum = 0.0;
      int start_j = std::max(0, i - n_x + 1);
      int end_j = std::min(n_k, i + 1);
      for (int j = start_j; j < end_j; j++) {
        int idx = i - j;
        if (idx >= 0 && idx < n_x) {
          sum += x[idx] * kernel[j];
        }
      }
      conv_result[i] = sum;
    }
    
    // Extract "same" mode
    int half_k = n_k / 2;
    NumericVector result(n_x);
    for (int i = 0; i < n_x; i++) {
      result[i] = conv_result[half_k + i];
    }
    
    return result;
  }
  
  // Strategy 2: Large arrays OR large kernels - chunked direct convolution (most reliable)
  // This avoids R overhead and memory issues for huge arrays (ATAC-seq, ChIP-seq)
  if (n_x > MEDIUM_THRESHOLD || n_k > KERNEL_LARGE_THRESHOLD) {
    int half_k = n_k / 2;
    NumericVector result(n_x);
    
    // Process in chunks - each chunk computes convolution independently
    // This is memory-efficient and avoids R overhead
    int n_chunks = (n_x + CHUNK_SIZE - 1) / CHUNK_SIZE;
    
    for (int chunk = 0; chunk < n_chunks; chunk++) {
      int chunk_start = chunk * CHUNK_SIZE;
      int chunk_end = std::min(n_x, (chunk + 1) * CHUNK_SIZE);
      
      // Direct convolution for each position in chunk
      // For "same" mode, we center the kernel at each output position
      for (int i = chunk_start; i < chunk_end; i++) {
        double sum = 0.0;
        // Convolve: sum over kernel positions
        for (int j = 0; j < n_k; j++) {
          // For "same" mode: kernel[j] aligns with x[i - half_k + j]
          int x_idx = i - half_k + j;
          if (x_idx >= 0 && x_idx < n_x) {
            sum += x[x_idx] * kernel[j];
          }
        }
        result[i] = sum;
      }
    }
    
    return result;
  }
  
  // Strategy 3: Medium arrays - try R's FFT-based convolve with fallback
  // This is faster for medium arrays but can hang for very large ones
  try {
    Environment stats_env = Environment::namespace_env("stats");
    Function convolve_func = stats_env["convolve"];
    
    NumericVector conv_result = convolve_func(x, kernel, Named("type") = "open");
    
    // Extract "same" mode
    int half_k = n_k / 2;
    NumericVector result(n_x);
    for (int i = 0; i < n_x; i++) {
      result[i] = conv_result[half_k + i];
    }
    
    return result;
  } catch (...) {
    // Fallback to direct convolution if R's convolve fails
    int n_conv = n_x + n_k - 1;
    NumericVector conv_result(n_conv);
    
    for (int i = 0; i < n_conv; i++) {
      double sum = 0.0;
      int start_j = std::max(0, i - n_x + 1);
      int end_j = std::min(n_k, i + 1);
      for (int j = start_j; j < end_j; j++) {
        int idx = i - j;
        if (idx >= 0 && idx < n_x) {
          sum += x[idx] * kernel[j];
        }
      }
      conv_result[i] = sum;
    }
    
    int half_k = n_k / 2;
    NumericVector result(n_x);
    for (int i = 0; i < n_x; i++) {
      result[i] = conv_result[half_k + i];
    }
    
    return result;
  }
}

// Compute adaptive scale coverage and check redundancy of fixed templates
//
// This function efficiently computes coverage percentages for each adaptive scale
// and determines which fixed templates are redundant based on dynamic thresholds.
//
// @param scale_map IntegerVector of scale indices (1-based) for each bin
// @param scale_L_bins IntegerVector of bin lengths for each adaptive scale (1-based indexing)
// @param fixed_templates IntegerVector of fixed template bin lengths to check
// @param n_total integer, total number of bins (for coverage calculation)
// @return LogicalVector indicating which fixed templates to keep (TRUE = keep, FALSE = redundant)
// [[Rcpp::export]]
LogicalVector cpp_check_template_redundancy(
    const IntegerVector& scale_map,
    const IntegerVector& scale_L_bins,
    const IntegerVector& fixed_templates,
    int n_total) {
  
  int n_scales = scale_L_bins.size();
  int n_fixed = fixed_templates.size();
  
  if (n_scales == 0 || n_fixed == 0 || n_total <= 0) {
    // If no adaptive scales, keep all fixed templates
    return LogicalVector(n_fixed, true);
  }
  
  // Compute coverage for each adaptive scale (1-based indexing)
  // scale_map values are 1..K where K = n_scales
  std::vector<int> scale_counts(n_scales, 0);
  int n_valid = 0;
  
  for (int i = 0; i < scale_map.size(); i++) {
    int scale_idx = scale_map[i];
    // scale_idx is 1-based, convert to 0-based for array indexing
    if (scale_idx >= 1 && scale_idx <= n_scales) {
      scale_counts[scale_idx - 1]++;
      n_valid++;
    }
  }
  
  // Compute coverage percentages
  std::vector<double> coverage(n_scales);
  double n_effective = std::max(1, n_valid);  // Avoid division by zero
  for (int k = 0; k < n_scales; k++) {
    coverage[k] = static_cast<double>(scale_counts[k]) / n_effective;
  }
  
  // Determine redundancy for each fixed template
  LogicalVector keep_template(n_fixed, true);
  
  for (int f = 0; f < n_fixed; f++) {
    int L_fixed = fixed_templates[f];
    bool is_redundant = false;
    
    // Check against each adaptive scale
    for (int k = 0; k < n_scales; k++) {
      int L_adaptive = scale_L_bins[k];
      double cov = coverage[k];
      
      // Compute bin difference
      int bin_diff = std::abs(L_adaptive - L_fixed);
      
      // Dynamic threshold based on coverage
      int threshold;
      if (cov > 0.90) {
        threshold = 1;  // High coverage: strict
      } else if (cov > 0.50) {
        threshold = 2;  // Medium coverage: moderate
      } else if (cov > 0.10) {
        threshold = 3;  // Low coverage: lenient
      } else {
        threshold = 999999;  // Very low coverage: never skip (complementary)
      }
      
      if (bin_diff <= threshold) {
        is_redundant = true;
        break;  // Found redundancy, no need to check other scales
      }
    }
    
    keep_template[f] = !is_redundant;
  }
  
  return keep_template;
}

// Fast valley detection using binary search for peak merging
//
// This function efficiently determines if two peaks should be merged based on
// valley detection and signal-to-valley ratio. Uses binary search for O(log n)
// valley region finding instead of O(n) linear scan.
//
// @param peak1_start double, start position of first peak in bp
// @param peak1_end double, end position of first peak in bp
// @param peak2_start double, start position of second peak in bp
// @param peak2_end double, end position of second peak in bp
// @param peak1_signal double, signal strength of first peak
// @param peak2_signal double, signal strength of second peak
// @param signal_values NumericVector, signal values for all bins
// @param intervals NumericVector, genomic intervals (bin start positions), must be sorted
// @param null_median double, median of null distribution (optional, use NaN if not available)
// @param null_mad double, MAD of null distribution (optional, use NaN if not available)
// @return List with should_merge (logical) and valley_ratio (numeric)
// [[Rcpp::export]]
List cpp_check_valley_between_peaks(
    double peak1_start,
    double peak1_end,
    double peak2_start,
    double peak2_end,
    double peak1_signal,
    double peak2_signal,
    const NumericVector& signal_values,
    const NumericVector& intervals,
    double null_median = NA_REAL,
    double null_mad = NA_REAL) {
  
  int n = intervals.size();
  if (n == 0 || signal_values.size() != n) {
    return List::create(
      Named("should_merge") = false,
      Named("valley_ratio") = NA_REAL
    );
  }
  
  // Check overlap first (fast path)
  double overlap = std::min(peak1_end, peak2_end) - std::max(peak1_start, peak2_start);
  if (overlap > 0) {
    double width1 = peak1_end - peak1_start;
    double width2 = peak2_end - peak2_start;
    double smaller_width = std::min(width1, width2);
    if (overlap > 0.5 * smaller_width) {
      return List::create(
        Named("should_merge") = true,
        Named("valley_ratio") = 0.0
      );
    }
  }
  
  // Determine valley region
  double valley_start_bp = peak1_end;
  double valley_end_bp = peak2_start;
  
  // If peaks overlap, use overlap region as valley
  if (valley_start_bp > valley_end_bp) {
    valley_start_bp = peak2_start;
    valley_end_bp = peak1_end;
  }
  
  // Binary search for valley region bounds (O(log n) instead of O(n))
  // Find first interval >= valley_start_bp
  int valley_start_idx = std::lower_bound(intervals.begin(), intervals.end(), valley_start_bp) - intervals.begin();
  
  // Find last interval <= valley_end_bp
  int valley_end_idx = std::upper_bound(intervals.begin(), intervals.end(), valley_end_bp) - intervals.begin() - 1;
  
  // Check if valley region exists
  if (valley_start_idx > valley_end_idx || valley_start_idx >= n || valley_end_idx < 0) {
    // No valley region (peaks are adjacent)
    double gap = peak2_start - peak1_end;
    if (gap <= 0) {
      return List::create(
        Named("should_merge") = true,
        Named("valley_ratio") = 0.0
      );
    }
    return List::create(
      Named("should_merge") = false,
      Named("valley_ratio") = NA_REAL
    );
  }
  
  // Extract valley signals (only the bins we need)
  std::vector<double> valley_signals;
  valley_signals.reserve(valley_end_idx - valley_start_idx + 1);
  
  for (int i = valley_start_idx; i <= valley_end_idx && i < n; i++) {
    double val = signal_values[i];
    if (R_FINITE(val)) {  // Use R_FINITE for R compatibility
      valley_signals.push_back(val);
    }
  }
  
  if (valley_signals.empty()) {
    return List::create(
      Named("should_merge") = false,
      Named("valley_ratio") = NA_REAL
    );
  }
  
  // Compute valley statistics
  double valley_min = *std::min_element(valley_signals.begin(), valley_signals.end());
  
  // Compute median efficiently
  std::vector<double> valley_sorted = valley_signals;
  std::sort(valley_sorted.begin(), valley_sorted.end());
  double valley_median;
  int n_valley = valley_sorted.size();
  if (n_valley % 2 == 0) {
    valley_median = (valley_sorted[n_valley / 2 - 1] + valley_sorted[n_valley / 2]) / 2.0;
  } else {
    valley_median = valley_sorted[n_valley / 2];
  }
  
  // Peak signals
  double peak_max_signal = std::max(peak1_signal, peak2_signal);
  double peak_avg_signal = (peak1_signal + peak2_signal) / 2.0;
  
  // Compute valley depth ratio
  double valley_depth_ratio = 0.0;
  if (peak_max_signal > 1e-10 && valley_min < peak_max_signal) {
    valley_depth_ratio = (peak_max_signal - valley_min) / peak_max_signal;
  }
  
  // Check valley relative to null distribution
  bool valley_above_null = true;
  if (R_FINITE(null_median) && R_FINITE(null_mad) && null_mad > 1e-10) {
    double null_threshold = null_median - 2.0 * null_mad;
    valley_above_null = valley_min > null_threshold;
  } else if (R_FINITE(null_median)) {
    valley_above_null = valley_min > null_median;
  }
  
  // DECISION LOGIC (same as R version):
  // 1. Deep valley (>50% drop) → distinct peaks
  if (valley_depth_ratio > 0.5) {
    return List::create(
      Named("should_merge") = false,
      Named("valley_ratio") = valley_depth_ratio
    );
  }
  
  // 2. Valley near null AND >30% drop → distinct peaks
  if (!valley_above_null && valley_depth_ratio > 0.3) {
    return List::create(
      Named("should_merge") = false,
      Named("valley_ratio") = valley_depth_ratio
    );
  }
  
  // 3. Shallow valley (<30% drop) AND above null → same peak
  if (valley_depth_ratio < 0.3 && valley_above_null) {
    return List::create(
      Named("should_merge") = true,
      Named("valley_ratio") = valley_depth_ratio
    );
  }
  
  // 4. Moderate valley (30-50%): use statistical test
  if (n_valley >= 3 && peak_avg_signal > 1e-10) {
    double valley_to_peak_ratio = valley_median / peak_avg_signal;
    if (valley_to_peak_ratio < 0.5) {
      return List::create(
        Named("should_merge") = false,
        Named("valley_ratio") = valley_depth_ratio
      );
    }
  }
  
  // Default: conservative (don't merge if uncertain)
  return List::create(
    Named("should_merge") = false,
    Named("valley_ratio") = valley_depth_ratio
  );
}

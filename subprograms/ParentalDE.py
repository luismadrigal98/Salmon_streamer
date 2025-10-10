#!/usr/bin/env python3

"""
Differential Expression Analysis Between Parental Lines

This program tests whether genes show differential expression between parental lines
(e.g., IM767 vs SF or IM767 vs SWB) using statistical tests.

@Author: Luis Javier Madrigal-Roca

@Date: 2025-10-10

"""

import sys
import os
import argparse
import numpy as np
from scipy import stats
import pandas as pd

# Add the parent directory to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def perform_ttest(parent1_values, parent2_values, gene_id):
    """
    Perform t-test between two parental lines for a gene.
    
    Args:
        parent1_values: List/array of expression values for parent 1
        parent2_values: List/array of expression values for parent 2
        gene_id: Gene identifier
        
    Returns:
        Dictionary with test results
    """
    # Remove any NaN or infinite values
    parent1_clean = np.array([x for x in parent1_values if np.isfinite(x)])
    parent2_clean = np.array([x for x in parent2_values if np.isfinite(x)])
    
    # Check if we have enough data
    if len(parent1_clean) < 2 or len(parent2_clean) < 2:
        return {
            'gene_id': gene_id,
            'parent1_mean': np.nan,
            'parent2_mean': np.nan,
            'parent1_sd': np.nan,
            'parent2_sd': np.nan,
            'parent1_n': len(parent1_clean),
            'parent2_n': len(parent2_clean),
            'fold_change': np.nan,
            'log2_fold_change': np.nan,
            't_statistic': np.nan,
            'p_value': np.nan,
            'status': 'insufficient_data'
        }
    
    # Calculate statistics
    parent1_mean = np.mean(parent1_clean)
    parent2_mean = np.mean(parent2_clean)
    parent1_sd = np.std(parent1_clean, ddof=1)
    parent2_sd = np.std(parent2_clean, ddof=1)
    
    # Calculate fold change (avoiding division by zero)
    if parent2_mean > 0:
        fold_change = parent1_mean / parent2_mean
        log2_fold_change = np.log2(fold_change) if fold_change > 0 else np.nan
    else:
        fold_change = np.nan
        log2_fold_change = np.nan
    
    # Perform t-test
    try:
        t_statistic, p_value = stats.ttest_ind(parent1_clean, parent2_clean, equal_var=False)
    except Exception as e:
        print(f"Warning: t-test failed for {gene_id}: {e}")
        return {
            'gene_id': gene_id,
            'parent1_mean': parent1_mean,
            'parent2_mean': parent2_mean,
            'parent1_sd': parent1_sd,
            'parent2_sd': parent2_sd,
            'parent1_n': len(parent1_clean),
            'parent2_n': len(parent2_clean),
            'fold_change': fold_change,
            'log2_fold_change': log2_fold_change,
            't_statistic': np.nan,
            'p_value': np.nan,
            'status': 'test_failed'
        }
    
    return {
        'gene_id': gene_id,
        'parent1_mean': parent1_mean,
        'parent2_mean': parent2_mean,
        'parent1_sd': parent1_sd,
        'parent2_sd': parent2_sd,
        'parent1_n': len(parent1_clean),
        'parent2_n': len(parent2_clean),
        'fold_change': fold_change,
        'log2_fold_change': log2_fold_change,
        't_statistic': t_statistic,
        'p_value': p_value,
        'status': 'tested'
    }


def perform_mannwhitney(parent1_values, parent2_values, gene_id):
    """
    Perform Mann-Whitney U test (non-parametric alternative to t-test).
    
    Args:
        parent1_values: List/array of expression values for parent 1
        parent2_values: List/array of expression values for parent 2
        gene_id: Gene identifier
        
    Returns:
        Dictionary with test results
    """
    # Remove any NaN or infinite values
    parent1_clean = np.array([x for x in parent1_values if np.isfinite(x)])
    parent2_clean = np.array([x for x in parent2_values if np.isfinite(x)])
    
    # Check if we have enough data
    if len(parent1_clean) < 2 or len(parent2_clean) < 2:
        return {
            'gene_id': gene_id,
            'parent1_median': np.nan,
            'parent2_median': np.nan,
            'parent1_n': len(parent1_clean),
            'parent2_n': len(parent2_clean),
            'u_statistic': np.nan,
            'p_value': np.nan,
            'status': 'insufficient_data'
        }
    
    # Calculate medians
    parent1_median = np.median(parent1_clean)
    parent2_median = np.median(parent2_clean)
    
    # Perform Mann-Whitney U test
    try:
        u_statistic, p_value = stats.mannwhitneyu(parent1_clean, parent2_clean, alternative='two-sided')
    except Exception as e:
        print(f"Warning: Mann-Whitney test failed for {gene_id}: {e}")
        return {
            'gene_id': gene_id,
            'parent1_median': parent1_median,
            'parent2_median': parent2_median,
            'parent1_n': len(parent1_clean),
            'parent2_n': len(parent2_clean),
            'u_statistic': np.nan,
            'p_value': np.nan,
            'status': 'test_failed'
        }
    
    return {
        'gene_id': gene_id,
        'parent1_median': parent1_median,
        'parent2_median': parent2_median,
        'parent1_n': len(parent1_clean),
        'parent2_n': len(parent2_clean),
        'u_statistic': u_statistic,
        'p_value': p_value,
        'status': 'tested'
    }


def benjamini_hochberg_correction(p_values):
    """
    Apply Benjamini-Hochberg FDR correction to p-values.
    
    Args:
        p_values: Array of p-values
        
    Returns:
        Array of adjusted p-values (q-values)
    """
    p_values = np.array(p_values)
    n = len(p_values)
    
    # Sort p-values and keep track of original indices
    sorted_indices = np.argsort(p_values)
    sorted_p_values = p_values[sorted_indices]
    
    # Calculate adjusted p-values
    adjusted_p_values = np.zeros(n)
    for i in range(n-1, -1, -1):
        if i == n-1:
            adjusted_p_values[i] = sorted_p_values[i]
        else:
            adjusted_p_values[i] = min(adjusted_p_values[i+1], 
                                      sorted_p_values[i] * n / (i + 1))
    
    # Restore original order
    q_values = np.zeros(n)
    q_values[sorted_indices] = adjusted_p_values
    
    return q_values


def parse_expression_file(expression_file, parent1_samples, parent2_samples):
    """
    Parse expression file and extract data for parental samples.
    
    Args:
        expression_file: Path to expression file (tab-delimited)
        parent1_samples: List of sample IDs for parent 1
        parent2_samples: List of sample IDs for parent 2
        
    Returns:
        Dictionary with gene expression data
    """
    print(f"Reading expression file: {expression_file}")
    
    # Read the file
    df = pd.read_csv(expression_file, sep='\t')
    
    # Get the gene ID column (usually first column)
    gene_col = df.columns[0]
    
    # Find which columns correspond to our samples
    parent1_cols = [col for col in df.columns if any(p1 in col for p1 in parent1_samples)]
    parent2_cols = [col for col in df.columns if any(p2 in col for p2 in parent2_samples)]
    
    if not parent1_cols:
        raise ValueError(f"No columns found matching parent 1 samples: {parent1_samples}")
    if not parent2_cols:
        raise ValueError(f"No columns found matching parent 2 samples: {parent2_samples}")
    
    print(f"Found {len(parent1_cols)} columns for parent 1: {parent1_cols}")
    print(f"Found {len(parent2_cols)} columns for parent 2: {parent2_cols}")
    
    # Extract data
    expression_data = {}
    for idx, row in df.iterrows():
        gene_id = row[gene_col]
        parent1_values = row[parent1_cols].values.astype(float)
        parent2_values = row[parent2_cols].values.astype(float)
        
        expression_data[gene_id] = {
            'parent1': parent1_values,
            'parent2': parent2_values
        }
    
    print(f"Loaded expression data for {len(expression_data)} genes")
    return expression_data


def main(args):
    """
    Main function to perform differential expression analysis between parental lines.
    """
    print("="*80)
    print("Parental Line Differential Expression Analysis")
    print("="*80)
    
    # Parse sample lists
    parent1_samples = args.parent1_samples
    parent2_samples = args.parent2_samples
    
    print(f"\nParent 1 ({args.parent1_name}): {parent1_samples}")
    print(f"Parent 2 ({args.parent2_name}): {parent2_samples}")
    
    # Parse expression data
    expression_data = parse_expression_file(args.expression_file, parent1_samples, parent2_samples)
    
    # Perform statistical tests
    print(f"\nPerforming {args.test_method} test for {len(expression_data)} genes...")
    
    results = []
    for gene_id, data in expression_data.items():
        if args.test_method == 'ttest':
            result = perform_ttest(data['parent1'], data['parent2'], gene_id)
        elif args.test_method == 'mannwhitney':
            result = perform_mannwhitney(data['parent1'], data['parent2'], gene_id)
        else:
            raise ValueError(f"Unknown test method: {args.test_method}")
        
        results.append(result)
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Apply multiple testing correction
    if args.test_method == 'ttest':
        valid_p_values = results_df['p_value'].notna()
        if valid_p_values.sum() > 0:
            print("\nApplying Benjamini-Hochberg FDR correction...")
            q_values = np.full(len(results_df), np.nan)
            q_values[valid_p_values] = benjamini_hochberg_correction(
                results_df.loc[valid_p_values, 'p_value'].values
            )
            results_df['q_value'] = q_values
            
            # Classify genes based on significance and direction
            results_df['significant'] = (results_df['q_value'] < args.fdr_threshold).astype(str)
            results_df.loc[results_df['q_value'].isna(), 'significant'] = 'NA'
            
            # Determine direction of change
            results_df['direction'] = 'none'
            significant_mask = results_df['q_value'] < args.fdr_threshold
            results_df.loc[significant_mask & (results_df['log2_fold_change'] > 0), 'direction'] = f'{args.parent1_name}_higher'
            results_df.loc[significant_mask & (results_df['log2_fold_change'] < 0), 'direction'] = f'{args.parent2_name}_higher'
    
    # Sort by p-value
    if args.test_method == 'ttest':
        results_df = results_df.sort_values('p_value')
    else:
        results_df = results_df.sort_values('p_value')
    
    # Save results
    output_file = os.path.join(args.output_dir, f"{args.parent1_name}_vs_{args.parent2_name}_DE_results.txt")
    os.makedirs(args.output_dir, exist_ok=True)
    results_df.to_csv(output_file, sep='\t', index=False, float_format='%.6g')
    print(f"\nResults saved to: {output_file}")
    
    # Print summary statistics
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    tested = (results_df['status'] == 'tested').sum()
    print(f"Total genes analyzed: {len(results_df)}")
    print(f"Genes successfully tested: {tested}")
    
    if args.test_method == 'ttest':
        significant = (results_df['q_value'] < args.fdr_threshold).sum()
        parent1_higher = ((results_df['q_value'] < args.fdr_threshold) & 
                         (results_df['log2_fold_change'] > 0)).sum()
        parent2_higher = ((results_df['q_value'] < args.fdr_threshold) & 
                         (results_df['log2_fold_change'] < 0)).sum()
        
        print(f"\nSignificant genes (FDR < {args.fdr_threshold}): {significant}")
        print(f"  - Higher in {args.parent1_name}: {parent1_higher}")
        print(f"  - Higher in {args.parent2_name}: {parent2_higher}")
        
        # Show top differentially expressed genes
        top_genes = results_df[results_df['q_value'] < args.fdr_threshold].head(10)
        if len(top_genes) > 0:
            print(f"\nTop 10 differentially expressed genes:")
            print(top_genes[['gene_id', 'parent1_mean', 'parent2_mean', 
                           'log2_fold_change', 'p_value', 'q_value', 'direction']].to_string(index=False))
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Perform differential expression analysis between parental lines',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # T-test between IM767 and SF parental lines
  python ParentalDE.py \\
    --expression-file combined_counts.txt \\
    --parent1-samples IM767-001 IM767-002 IM767-003 \\
    --parent2-samples SF-001 SF-002 SF-003 \\
    --parent1-name IM767 \\
    --parent2-name SF \\
    --output-dir DE_results

  # Mann-Whitney test (non-parametric)
  python ParentalDE.py \\
    --expression-file combined_counts.txt \\
    --parent1-samples IM767 \\
    --parent2-samples SF \\
    --parent1-name IM767 \\
    --parent2-name SF \\
    --test-method mannwhitney \\
    --output-dir DE_results
        """
    )
    
    # Required arguments
    parser.add_argument('--expression-file', required=True,
                       help='Path to expression data file (tab-delimited, genes in rows)')
    parser.add_argument('--parent1-samples', nargs='+', required=True,
                       help='Sample IDs or patterns for parent 1 (space-separated)')
    parser.add_argument('--parent2-samples', nargs='+', required=True,
                       help='Sample IDs or patterns for parent 2 (space-separated)')
    parser.add_argument('--parent1-name', required=True,
                       help='Name for parent 1 (e.g., IM767)')
    parser.add_argument('--parent2-name', required=True,
                       help='Name for parent 2 (e.g., SF, SWB)')
    
    # Optional arguments
    parser.add_argument('--test-method', choices=['ttest', 'mannwhitney'], default='ttest',
                       help='Statistical test method (default: ttest)')
    parser.add_argument('--fdr-threshold', type=float, default=0.05,
                       help='FDR threshold for significance (default: 0.05)')
    parser.add_argument('--output-dir', default='./DE_results',
                       help='Output directory for results (default: ./DE_results)')
    
    args = parser.parse_args()
    main(args)

#!/usr/bin/env python3
"""
Complete SV Analysis Visualization Suite
Generates all plots for structural variant analysis pipeline with proper data handling
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import numpy as np
import json
import sys
import os
import glob
from collections import defaultdict, Counter
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
import matplotlib.gridspec as gridspec
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Set style for clarity
plt.style.use('default')
sns.set_palette("Set2")
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['font.size'] = 10

def ensure_dir(directory):
    """Ensure directory exists"""
    if not os.path.exists(directory):
        os.makedirs(directory)

def find_vcf_files(output_dir):
    """Find all VCF files and map them to caller names"""
    vcf_files = {}
    
    # Look for specific patterns
    patterns = {
        'truth': ['truth.vcf'],
        'enhanced': ['enhanced_calls.vcf', 'complex-sv_calls.vcf', 'caller_calls.vcf'],
        'sniffles2': ['sniffles2_calls.vcf'],
        'cutesv': ['cutesv_calls.vcf'],
        'svim': ['svim_calls.vcf']
    }
    
    for caller, file_patterns in patterns.items():
        for pattern in file_patterns:
            filepath = os.path.join(output_dir, pattern)
            if os.path.exists(filepath):
                vcf_files[caller] = filepath
                print(f"Found {caller}: {filepath}")
                break
    
    return vcf_files

def parse_vcf(vcf_file):
    """Parse VCF file with enhanced DUP subtype handling"""
    svs = []
    
    if not os.path.exists(vcf_file):
        print(f"Warning: VCF file not found: {vcf_file}")
        return svs
    
    try:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                chrom = parts[0]
                pos = int(parts[1]) - 1
                sv_id = parts[2]
                alt = parts[4]
                info = parts[7]
                filter_field = parts[6]
                qual = float(parts[5]) if parts[5] != '.' and parts[5] != '0' else 0
                
                # Parse genotype if available
                genotype = './.'
                gt_category = 'Unknown'
                if len(parts) >= 10:
                    format_field = parts[8]
                    sample = parts[9]
                    
                    if 'GT' in format_field and sample:
                        format_parts = format_field.split(':')
                        sample_parts = sample.split(':')
                        
                        if len(format_parts) == len(sample_parts):
                            try:
                                gt_index = format_parts.index('GT')
                                genotype = sample_parts[gt_index]
                                
                                if genotype in ['0/1', '1/0']:
                                    gt_category = 'Heterozygous (0/1)'
                                elif genotype == '1/1':
                                    gt_category = 'Homozygous (1/1)'
                                elif genotype == '0/0':
                                    gt_category = 'Reference (0/0)'
                                else:
                                    gt_category = 'Unknown'
                            except (ValueError, IndexError):
                                pass
                
                # Parse INFO field
                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                    else:
                        info_dict[item] = True
                
                sv_type = info_dict.get('SVTYPE', 'Unknown')
                
                # Enhanced DUP subtype handling
                dup_subtype = None
                original_type = sv_type
                
                if sv_type == 'DUP' or 'DUP' in alt:
                    sv_type = 'DUP'
                    
                    if alt == '<DUP:TANDEM>' or 'DUP:TANDEM' in alt:
                        dup_subtype = 'TANDEM'
                    elif alt == '<DUP:INT>' or alt == '<DUP:INTERSPERSED>' or 'DUP:INT' in alt:
                        dup_subtype = 'INTERSPERSED'
                    elif alt == '<DUP>' and sv_type == 'DUP':
                        if 'DUPTYPE' in info_dict:
                            duptype = info_dict['DUPTYPE'].upper()
                            if 'INTERSPERSED' in duptype or 'INT' in duptype:
                                dup_subtype = 'INTERSPERSED'
                            else:
                                dup_subtype = 'TANDEM'
                        elif 'subtype' in info_dict:
                            subtype_val = info_dict['subtype'].upper()
                            if 'INTERSPERSED' in subtype_val:
                                dup_subtype = 'INTERSPERSED'
                            else:
                                dup_subtype = 'TANDEM'
                        elif 'INTERSPERSED' in info.upper():
                            dup_subtype = 'INTERSPERSED'
                        elif 'TANDEM' in info.upper():
                            dup_subtype = 'TANDEM'
                        else:
                            dup_subtype = 'TANDEM'
                    else:
                        dup_subtype = 'TANDEM'
                
                end = int(info_dict.get('END', pos + 1)) - 1
                size = abs(int(info_dict.get('SVLEN', end - pos)) if info_dict.get('SVLEN', '0') != '.' else end - pos)
                support = int(info_dict.get('SUPPORT', 0))
                af = float(info_dict.get('AF', 0))
                
                # Detect repeat information
                in_repeat = False
                repeat_type = None
                
                if 'REPEAT' in info_dict:
                    repeat_value = info_dict['REPEAT'].lower()
                    if repeat_value in ['yes', 'true', '1']:
                        in_repeat = True
                        repeat_type = info_dict.get('REPEAT_TYPE', 'Unknown')
                    elif repeat_value in ['tandem', 'homopolymer', 'low_complexity', 'low complexity']:
                        in_repeat = True
                        repeat_type = repeat_value
                    elif repeat_value != 'no' and repeat_value != 'false' and repeat_value != '0':
                        in_repeat = True
                        repeat_type = repeat_value
                
                # Check for complex SV
                is_complex = sv_type == 'CPX' or 'CPX_TYPE' in info_dict
                complex_components = info_dict.get('CPX_TYPE', '').split(',') if is_complex else []
                
                svs.append({
                    'chrom': chrom,
                    'pos': pos,
                    'end': end,
                    'id': sv_id,
                    'type': sv_type,
                    'original_type': original_type,
                    'alt_field': alt,
                    'size': size,
                    'support': support,
                    'af': af,
                    'filter': filter_field,
                    'qual': qual,
                    'in_repeat': in_repeat,
                    'repeat_type': repeat_type,
                    'is_complex': is_complex,
                    'complex_components': complex_components,
                    'dup_subtype': dup_subtype,
                    'dup_type': dup_subtype,
                    'genotype': genotype,
                    'gt_category': gt_category,
                    'key': f"{chrom}:{pos//1000}:{sv_type}"
                })
    except Exception as e:
        print(f"Error parsing {vcf_file}: {e}")
    
    return svs

def load_truvari_results(eval_dir):
    """Load Truvari evaluation results"""
    results = {}
    
    for caller in ['enhanced', 'complex-sv', 'sniffles2', 'cutesv', 'svim']:
        summary_file = f"{eval_dir}/truvari_{caller}/summary.json"
        
        if os.path.exists(summary_file):
            try:
                with open(summary_file, 'r') as f:
                    data = json.load(f)
                    results[caller] = {
                        'precision': float(data.get('precision', 0) or 0),
                        'recall': float(data.get('recall', 0) or 0),
                        'f1': float(data.get('f1', 0) or 0),
                        'tp': int(data.get('TP-base', data.get('TP', 0))),
                        'fp': int(data.get('FP', 0)),
                        'fn': int(data.get('FN', 0))
                    }
            except Exception as e:
                print(f"Could not load results for {caller}: {e}")
    
    return results

def load_runtime_data(output_dir):
    """Load runtime data"""
    runtime_file = f"{output_dir}/runtime/performance_summary.csv"
    runtime_data = {}
    memory_data = {}
    
    if os.path.exists(runtime_file):
        try:
            df = pd.read_csv(runtime_file)
            for _, row in df.iterrows():
                stage = row['Stage'].lower()
                if 'enhanced' in stage or 'caller' in stage:
                    runtime_data['enhanced'] = row['Runtime_seconds']
                    if 'Memory_GB' in row and pd.notna(row['Memory_GB']) and row['Memory_GB'] > 0:
                        memory_data['enhanced'] = row['Memory_GB']
                elif 'sniffles' in stage:
                    runtime_data['sniffles2'] = row['Runtime_seconds']
                    if 'Memory_GB' in row and pd.notna(row['Memory_GB']) and row['Memory_GB'] > 0:
                        memory_data['sniffles2'] = row['Memory_GB']
                elif 'cutesv' in stage:
                    runtime_data['cutesv'] = row['Runtime_seconds']
                    if 'Memory_GB' in row and pd.notna(row['Memory_GB']) and row['Memory_GB'] > 0:
                        memory_data['cutesv'] = row['Memory_GB']
                elif 'svim' in stage:
                    runtime_data['svim'] = row['Runtime_seconds']
                    if 'Memory_GB' in row and pd.notna(row['Memory_GB']) and row['Memory_GB'] > 0:
                        memory_data['svim'] = row['Memory_GB']
        except Exception as e:
            print(f"Could not load runtime data: {e}")
    
    return runtime_data, memory_data

def plot_1_performance_metrics(results, output_dir):
    """Plot 1: Performance metrics comparison"""
    try:
        fig, ax = plt.subplots(figsize=(12, 7))
        
        if not results:
            ax.text(0.5, 0.5, 'No Truvari evaluation results available\n\nRun Truvari bench to generate performance metrics', 
                   ha='center', va='center', fontsize=14)
            ax.set_title('Performance Metrics Comparison\n(Precision, Recall, F1 from Truvari evaluation)', 
                        fontsize=16, fontweight='bold')
        else:
            callers = list(results.keys())
            metrics = ['Precision', 'Recall', 'F1']
            
            x = np.arange(len(callers))
            width = 0.25
            colors = ['#3498db', '#2ecc71', '#e74c3c']
            
            for i, metric in enumerate(metrics):
                values = [results[c][metric.lower()] for c in callers]
                bars = ax.bar(x + i*width, values, width, label=metric, color=colors[i], alpha=0.8)
                
                for bar, val in zip(bars, values):
                    height = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                           f'{val:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
            
            ax.set_xlabel('SV Callers', fontsize=14, fontweight='bold')
            ax.set_ylabel('Score', fontsize=14, fontweight='bold')
            ax.set_title('Performance Metrics Comparison\n(Precision, Recall, F1 from Truvari evaluation)', 
                        fontsize=16, fontweight='bold')
            ax.set_xticks(x + width)
            ax.set_xticklabels([c.upper() for c in callers], fontsize=12)
            ax.legend(loc='upper right', fontsize=12)
            ax.set_ylim([0, 1.1])
            ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/1_performance_metrics.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 1_performance_metrics.png")
        return True
    except Exception as e:
        print(f"Error generating plot 1: {e}")
        return False

def plot_2_runtime_comparison(runtime_data, memory_data, output_dir):
    """Plot 2: Runtime comparison"""
    try:
        has_memory = len(memory_data) > 0
        
        if has_memory:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
        else:
            fig, ax1 = plt.subplots(figsize=(10, 7))
        
        if not runtime_data:
            ax1.text(0.5, 0.5, 'No runtime data available\n\nRuntime data is collected from\nperformance_summary.csv', 
                   ha='center', va='center', fontsize=14)
            ax1.set_title('Runtime Comparison', fontsize=14, fontweight='bold')
        else:
            callers = list(runtime_data.keys())
            runtimes = list(runtime_data.values())
            
            colors = ['#2ecc71', '#3498db', '#e74c3c', '#f39c12'][:len(callers)]
            
            bars = ax1.bar(range(len(callers)), runtimes, color=colors, alpha=0.7)
            ax1.set_xticks(range(len(callers)))
            ax1.set_xticklabels([c.replace('_', ' ').title() for c in callers], rotation=45, ha='right')
            ax1.set_ylabel('Runtime (seconds)', fontsize=12, fontweight='bold')
            ax1.set_title('Runtime Comparison\n(Time to complete SV calling)', fontsize=14, fontweight='bold')
            ax1.grid(True, alpha=0.3, axis='y')
            
            for bar, runtime in zip(bars, runtimes):
                ax1.text(bar.get_x() + bar.get_width()/2., runtime + max(runtimes)*0.01,
                        f'{runtime:.0f}s', ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        if has_memory:
            memory_callers = list(memory_data.keys())
            memory_values = list(memory_data.values())
            
            bars = ax2.bar(range(len(memory_callers)), memory_values, color=colors[:len(memory_callers)], alpha=0.7)
            ax2.set_xticks(range(len(memory_callers)))
            ax2.set_xticklabels([c.replace('_', ' ').title() for c in memory_callers], rotation=45, ha='right')
            ax2.set_ylabel('Memory Usage (GB)', fontsize=12, fontweight='bold')
            ax2.set_title('Memory Usage Comparison\n(Peak memory consumption)', fontsize=14, fontweight='bold')
            ax2.grid(True, alpha=0.3, axis='y')
            
            for bar, memory in zip(bars, memory_values):
                ax2.text(bar.get_x() + bar.get_width()/2., memory + max(memory_values)*0.01,
                        f'{memory:.1f}GB', ha='center', va='bottom', fontsize=10, fontweight='bold')
            
            plt.suptitle('Performance Metrics: Runtime and Memory Usage', fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/2_runtime_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 2_runtime_comparison.png")
        return True
    except Exception as e:
        print(f"Error generating plot 2: {e}")
        return False

def plot_3_tp_fp_fn_stacked(results, output_dir):
    """Plot 3: TP/FP/FN stacked bars"""
    try:
        fig, ax = plt.subplots(figsize=(12, 7))
        
        if not results:
            ax.text(0.5, 0.5, 'No Truvari evaluation results available\n\nRun Truvari bench to generate TP/FP/FN counts', 
                   ha='center', va='center', fontsize=14)
        else:
            callers = list(results.keys())
            tp = [results[c]['tp'] for c in callers]
            fp = [results[c]['fp'] for c in callers]
            fn = [results[c]['fn'] for c in callers]
            
            x = np.arange(len(callers))
            width = 0.6
            
            p1 = ax.bar(x, tp, width, label='True Positives', color='#2ecc71', alpha=0.9)
            p2 = ax.bar(x, fp, width, bottom=tp, label='False Positives', color='#e74c3c', alpha=0.9)
            p3 = ax.bar(x, fn, width, bottom=np.array(tp)+np.array(fp), label='False Negatives', color='#f39c12', alpha=0.9)
            
            for i, (t, f, n) in enumerate(zip(tp, fp, fn)):
                if t > 0:
                    ax.text(i, t/2, str(t), ha='center', va='center', fontsize=11, fontweight='bold', color='white')
                if f > 0:
                    ax.text(i, t + f/2, str(f), ha='center', va='center', fontsize=11, fontweight='bold', color='white')
                if n > 0:
                    ax.text(i, t + f + n/2, str(n), ha='center', va='center', fontsize=11, fontweight='bold', color='white')
            
            ax.set_xlabel('SV Callers', fontsize=14, fontweight='bold')
            ax.set_ylabel('Number of Calls', fontsize=14, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels([c.upper() for c in callers], fontsize=12)
            ax.legend(loc='upper right', fontsize=12)
            ax.grid(True, alpha=0.3, axis='y')
        
        ax.set_title('True Positives, False Positives, and False Negatives\n(From Truvari evaluation against truth set)', 
                    fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/3_tp_fp_fn_stacked.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 3_tp_fp_fn_stacked.png")
        return True
    except Exception as e:
        print(f"Error generating plot 3: {e}")
        return False

def plot_4_sv_types_comparison(truth_svs, enhanced_svs, output_dir):
    """Plot 4: SV types comparison"""
    try:
        fig, ax = plt.subplots(figsize=(14, 7))
        
        sv_types = ['DEL', 'INS', 'DUP', 'INV', 'BND', 'CPX']
        enhanced_pass = [sv for sv in enhanced_svs if sv['filter'] == 'PASS']
        
        truth_counts = {t: len([sv for sv in truth_svs if sv['type'] == t]) for t in sv_types}
        enhanced_counts = {t: len([sv for sv in enhanced_pass if sv['type'] == t]) for t in sv_types}
        
        x = np.arange(len(sv_types))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, [truth_counts[t] for t in sv_types], width, 
                       label='Simulated (Truth Set)', color='#3498db', alpha=0.9)
        bars2 = ax.bar(x + width/2, [enhanced_counts[t] for t in sv_types], width,
                       label='Detected (PASS calls)', color='#2ecc71', alpha=0.9)
        
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                           f'{int(height)}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        for i, sv_type in enumerate(sv_types):
            if truth_counts[sv_type] > 0:
                rate = enhanced_counts[sv_type] / truth_counts[sv_type] * 100
                ax.text(i, max(truth_counts[sv_type], enhanced_counts[sv_type]) + 5,
                       f'{rate:.1f}%', ha='center', va='bottom', fontsize=11,
                       fontweight='bold', color='red')
        
        ax.set_xlabel('SV Type', fontsize=14, fontweight='bold')
        ax.set_ylabel('Number of SVs', fontsize=14, fontweight='bold')
        ax.set_title('SV Type Detection Performance\n(Simulated vs Detected PASS calls by type)', 
                    fontsize=16, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(sv_types, fontsize=12)
        ax.legend(loc='upper right', fontsize=12)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/4_sv_types_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 4_sv_types_comparison.png")
        return True
    except Exception as e:
        print(f"Error generating plot 4: {e}")
        return False

def plot_5_duplication_pie_charts(truth_svs, enhanced_svs, output_dir):
    """Plot 5: Duplication type pie charts"""
    try:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
        
        # Truth duplications
        truth_dups = [sv for sv in truth_svs if sv['type'] == 'DUP']
        truth_tandem = 0
        truth_interspersed = 0
        
        for dup in truth_dups:
            subtype = dup.get('dup_subtype', dup.get('dup_type', ''))
            if 'INTERSPERSED' in str(subtype).upper():
                truth_interspersed += 1
            else:
                truth_tandem += 1
        
        # Plot 1: Truth Set
        if truth_tandem + truth_interspersed > 0:
            sizes1 = [truth_tandem, truth_interspersed]
            labels1 = [f'Tandem\n({truth_tandem})', f'Interspersed\n({truth_interspersed})']
            colors1 = ['#3498db', '#e74c3c']
            wedges1, texts1, autotexts1 = ax1.pie(sizes1, labels=labels1, colors=colors1, 
                                                  autopct='%1.1f%%', startangle=90, textprops={'fontsize': 11})
            for autotext in autotexts1:
                autotext.set_fontweight('bold')
        else:
            ax1.text(0.5, 0.5, 'No duplications\nin truth set', ha='center', va='center', fontsize=14)
        ax1.set_title('Truth Set Duplications\n(Simulated)', fontsize=14, fontweight='bold')
        
        # Enhanced caller duplications
        enhanced_dups = [sv for sv in enhanced_svs if sv['type'] == 'DUP' and sv['filter'] == 'PASS']
        enhanced_tandem = 0
        enhanced_interspersed = 0
        
        for dup in enhanced_dups:
            subtype = dup.get('dup_subtype', dup.get('dup_type', ''))
            if 'INTERSPERSED' in str(subtype).upper():
                enhanced_interspersed += 1
            else:
                enhanced_tandem += 1
        
        # Plot 2: Enhanced Caller
        if enhanced_tandem + enhanced_interspersed > 0:
            sizes2 = [enhanced_tandem, enhanced_interspersed]
            labels2 = [f'Tandem\n({enhanced_tandem})', f'Interspersed\n({enhanced_interspersed})']
            colors2 = ['#2ecc71', '#f39c12']
            wedges2, texts2, autotexts2 = ax2.pie(sizes2, labels=labels2, colors=colors2, 
                                                  autopct='%1.1f%%', startangle=90, textprops={'fontsize': 11})
            for autotext in autotexts2:
                autotext.set_fontweight('bold')
        else:
            ax2.text(0.5, 0.5, 'No duplications\ndetected', ha='center', va='center', fontsize=14)
        ax2.set_title('Enhanced Caller Duplications\n(PASS calls)', fontsize=14, fontweight='bold')
        
        # SVIM duplications
        ax3.text(0.5, 0.5, 'SVIM data\nnot available', ha='center', va='center', fontsize=14)
        ax3.set_title('SVIM Duplications\n(Not available)', fontsize=14, fontweight='bold')
        
        plt.suptitle('Duplication Type Distribution Comparison\n(Tandem vs Interspersed across callers)', 
                     fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/5_duplication_pie_charts.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 5_duplication_pie_charts.png")
        return True
    except Exception as e:
        print(f"Error generating plot 5: {e}")
        return False

def plot_6_repeat_svs_overview(truth_svs, enhanced_svs, output_dir):
    """Plot 6: Repeat SVs overview"""
    try:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
        
        # Truth statistics
        truth_total = len(truth_svs)
        truth_in_repeats = len([sv for sv in truth_svs if sv['in_repeat']])
        truth_not_in_repeats = truth_total - truth_in_repeats
        
        # Enhanced statistics (PASS only)
        enhanced_pass = [sv for sv in enhanced_svs if sv['filter'] == 'PASS']
        enhanced_total = len(enhanced_pass)
        enhanced_in_repeats = len([sv for sv in enhanced_pass if sv['in_repeat']])
        enhanced_not_in_repeats = enhanced_total - enhanced_in_repeats
        
        # Plot 1: Truth Set
        categories = ['Total', 'In Repeats', 'Not in Repeats']
        truth_counts = [truth_total, truth_in_repeats, truth_not_in_repeats]
        colors1 = ['#3498db', '#e67e22', '#95a5a6']
        
        bars1 = ax1.bar(categories, truth_counts, color=colors1, alpha=0.8)
        
        for bar, count in zip(bars1, truth_counts):
            ax1.text(bar.get_x() + bar.get_width()/2., count + max(truth_counts)*0.02,
                    str(count), ha='center', va='bottom', 
                    fontsize=12, fontweight='bold')
        
        if truth_total > 0:
            pct = (truth_in_repeats / truth_total * 100)
            ax1.text(0.5, 0.5, f'{pct:.1f}%\nin repeats', 
                    transform=ax1.transAxes, ha='center', va='center',
                    fontsize=14, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9))
        
        ax1.set_ylabel('Number of SVs', fontsize=12, fontweight='bold')
        ax1.set_title(f'Truth Set\n{truth_total} total SVs, {truth_in_repeats} in repeats', 
                     fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='y')
        
        # Plot 2: Enhanced Caller (PASS only)
        enhanced_counts = [enhanced_total, enhanced_in_repeats, enhanced_not_in_repeats]
        colors2 = ['#2ecc71', '#f39c12', '#7f8c8d']
        
        bars2 = ax2.bar(categories, enhanced_counts, color=colors2, alpha=0.8)
        
        for bar, count in zip(bars2, enhanced_counts):
            ax2.text(bar.get_x() + bar.get_width()/2., count + max(enhanced_counts)*0.02,
                    str(count), ha='center', va='bottom', 
                    fontsize=12, fontweight='bold')
        
        if enhanced_total > 0:
            pct = (enhanced_in_repeats / enhanced_total * 100)
            ax2.text(0.5, 0.5, f'{pct:.1f}%\nin repeats', 
                    transform=ax2.transAxes, ha='center', va='center',
                    fontsize=14, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9))
        
        ax2.set_ylabel('Number of SVs', fontsize=12, fontweight='bold')
        ax2.set_title(f'Enhanced Caller (PASS)\n{enhanced_total} PASS SVs, {enhanced_in_repeats} in repeats', 
                     fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='y')
        
        plt.suptitle('SV Distribution in Repeat Regions\n(Comparison of Truth Set vs Enhanced Caller)', 
                    fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/6_repeat_svs_overview.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 6_repeat_svs_overview.png")
        return True
        
    except Exception as e:
        print(f"Error generating plot 6: {e}")
        return False

def plot_7_detected_repeat_types(truth_svs, enhanced_svs, output_dir):
    """Plot 7: Simulated vs Detected SVs in Repeat Regions by Type"""
    try:
        fig, ax = plt.subplots(figsize=(14, 8))
        
        sv_types = ['DEL', 'INS', 'DUP', 'INV', 'BND', 'CPX']
        enhanced_pass = [sv for sv in enhanced_svs if sv['filter'] == 'PASS']
        
        # Count simulated and detected SVs in repeats
        simulated_counts = []
        detected_counts = []
        
        for sv_type in sv_types:
            simulated_type_repeats = len([sv for sv in truth_svs 
                                         if sv['type'] == sv_type and sv['in_repeat']])
            simulated_counts.append(simulated_type_repeats)
            
            detected_type_repeats = len([sv for sv in enhanced_pass 
                                       if sv['type'] == sv_type and sv['in_repeat']])
            detected_counts.append(detected_type_repeats)
        
        x = np.arange(len(sv_types))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, simulated_counts, width, 
                      label='Simulated (Truth)', color='#3498db', alpha=0.8)
        bars2 = ax.bar(x + width/2, detected_counts, width,
                      label='Detected (PASS)', color='#2ecc71', alpha=0.8)
        
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                           f'{int(height)}', ha='center', va='bottom', 
                           fontsize=10, fontweight='bold')
        
        for i, sv_type in enumerate(sv_types):
            if simulated_counts[i] > 0:
                rate = (detected_counts[i] / simulated_counts[i]) * 100
                ax.text(i, max(simulated_counts[i], detected_counts[i]) + 2,
                       f'{rate:.0f}%', ha='center', va='bottom', 
                       fontsize=11, fontweight='bold', color='red')
        
        ax.set_xlabel('SV Type', fontsize=14, fontweight='bold')
        ax.set_ylabel('Number of SVs in Repeat Regions', fontsize=14, fontweight='bold')
        ax.set_title('Simulated vs Detected SVs in Repeat Regions\n(Truth vs PASS calls in repeat regions by type)', 
                    fontsize=16, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(sv_types, fontsize=12)
        ax.legend(loc='upper right', fontsize=12)
        ax.grid(True, alpha=0.3, axis='y')
        
        total_sim = sum(simulated_counts)
        total_det = sum(detected_counts)
        if total_sim > 0:
            overall_rate = (total_det / total_sim) * 100
            ax.text(0.02, 0.98, f"Overall: {total_det}/{total_sim} ({overall_rate:.1f}%)", 
                   transform=ax.transAxes, fontsize=12, fontweight='bold', va='top')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/7_sv_repeat_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 7_sv_repeat_comparison.png")
        return True
    except Exception as e:
        print(f"Error generating plot 7: {e}")
        return False

def plot_8_repeat_types_distribution(truth_svs, enhanced_svs, output_dir):
    """Plot 8: Repeat types distribution"""
    try:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # Truth set: SVs in repeat regions
        truth_repeat_svs = [sv for sv in truth_svs if sv['in_repeat']]
        truth_repeat_counts = defaultdict(int)
        
        for sv in truth_repeat_svs:
            if sv['repeat_type']:
                repeat_type = sv['repeat_type'].lower()
                if 'homopolymer' in repeat_type:
                    truth_repeat_counts['Homopolymer'] += 1
                elif 'tandem' in repeat_type:
                    truth_repeat_counts['Tandem Repeat'] += 1
                elif 'low_complexity' in repeat_type or 'low' in repeat_type:
                    truth_repeat_counts['Low Complexity'] += 1
                else:
                    truth_repeat_counts['Other'] += 1
            else:
                truth_repeat_counts['Unspecified Repeat'] += 1
        
        # Enhanced: PASS SVs in repeat regions  
        enhanced_pass = [sv for sv in enhanced_svs if sv['filter'] == 'PASS']
        enhanced_repeat_svs = [sv for sv in enhanced_pass if sv['in_repeat']]
        enhanced_repeat_counts = defaultdict(int)
        
        for sv in enhanced_repeat_svs:
            if sv['repeat_type']:
                repeat_type = sv['repeat_type'].lower()
                if 'homopolymer' in repeat_type:
                    enhanced_repeat_counts['Homopolymer'] += 1
                elif 'tandem' in repeat_type:
                    enhanced_repeat_counts['Tandem Repeat'] += 1
                elif 'low_complexity' in repeat_type or 'low' in repeat_type:
                    enhanced_repeat_counts['Low Complexity'] += 1
                else:
                    enhanced_repeat_counts['Other'] += 1
            else:
                enhanced_repeat_counts['Unspecified Repeat'] += 1
        
        all_repeat_types = list(set(list(truth_repeat_counts.keys()) + list(enhanced_repeat_counts.keys())))
        colors = ['#e74c3c', '#3498db', '#f39c12', '#2ecc71', '#9b59b6']
        
        # Plot 1: Truth Set
        if truth_repeat_counts:
            truth_types = [rt for rt in all_repeat_types if truth_repeat_counts[rt] > 0]
            truth_counts = [truth_repeat_counts[rt] for rt in truth_types]
            
            bars1 = ax1.bar(truth_types, truth_counts, 
                           color=colors[:len(truth_types)], alpha=0.8)
            
            for bar, count in zip(bars1, truth_counts):
                if count > 0:
                    ax1.text(bar.get_x() + bar.get_width()/2., count + max(truth_counts)*0.02, 
                           str(count), ha='center', va='bottom', fontsize=11, fontweight='bold')
            
            ax1.set_xlabel('Repeat Region Type', fontsize=12, fontweight='bold')
            ax1.set_ylabel('Number of Simulated SVs', fontsize=12, fontweight='bold')
            ax1.set_title(f'Truth Set: SVs in Repeat Regions\n({len(truth_repeat_svs)} total repeat-associated variants)', 
                         fontsize=14, fontweight='bold')
            ax1.grid(True, alpha=0.3, axis='y')
            ax1.tick_params(axis='x', rotation=45, labelsize=10)
            
        else:
            ax1.text(0.5, 0.5, 'No repeat-associated SVs\nin truth set', 
                   ha='center', va='center', fontsize=14,
                   bbox=dict(boxstyle='round,pad=1', facecolor='lightgray', alpha=0.8))
            ax1.set_xlim(0, 1)
            ax1.set_ylim(0, 1)
        
        # Plot 2: Enhanced SVs
        if enhanced_repeat_counts:
            enhanced_types = [rt for rt in all_repeat_types if enhanced_repeat_counts[rt] > 0]
            enhanced_counts = [enhanced_repeat_counts[rt] for rt in enhanced_types]
            
            bars2 = ax2.bar(enhanced_types, enhanced_counts, 
                           color=colors[:len(enhanced_types)], alpha=0.8)
            
            for bar, count in zip(bars2, enhanced_counts):
                if count > 0:
                    ax2.text(bar.get_x() + bar.get_width()/2., count + max(enhanced_counts)*0.02, 
                           str(count), ha='center', va='bottom', fontsize=11, fontweight='bold')
            
            ax2.set_xlabel('Repeat Region Type', fontsize=12, fontweight='bold')
            ax2.set_ylabel('Number of Detected PASS SVs', fontsize=12, fontweight='bold')
            ax2.set_title(f'Enhanced Caller: PASS SVs in Repeat Regions\n({len(enhanced_repeat_svs)} total detected repeat variants)', 
                         fontsize=14, fontweight='bold')
            ax2.grid(True, alpha=0.3, axis='y')
            ax2.tick_params(axis='x', rotation=45, labelsize=10)
            
        else:
            ax2.text(0.5, 0.5, 'No repeat-associated SVs\ndetected as PASS', 
                   ha='center', va='center', fontsize=14,
                   bbox=dict(boxstyle='round,pad=1', facecolor='lightgray', alpha=0.8))
            ax2.set_xlim(0, 1)
            ax2.set_ylim(0, 1)
        
        if truth_repeat_svs and enhanced_repeat_svs:
            detection_rate = (len(enhanced_repeat_svs) / len(truth_repeat_svs)) * 100
            fig.suptitle(f'Repeat Region SV Detection: Truth vs Enhanced Caller\nOverall repeat detection rate: {detection_rate:.1f}% ({len(enhanced_repeat_svs)}/{len(truth_repeat_svs)})', 
                        fontsize=16, fontweight='bold')
        else:
            fig.suptitle('Repeat Region SV Detection: Truth vs Enhanced Caller\n(Simulated vs detected variants in repeat regions)', 
                        fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/8_repeat_region_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 8_repeat_region_comparison.png")
        return True
        
    except Exception as e:
        print(f"Error generating plot 8: {e}")
        return False

def plot_9_caller_overlap_venn(vcf_files, output_dir):
    """Plot 9: Caller overlap using Venn diagram"""
    try:
        # Try to install venn if not available
        try:
            import venn
        except ImportError:
            print("Installing venn package for 4-way diagram...")
            import subprocess
            subprocess.check_call([sys.executable, "-m", "pip", "install", "venn"])
            import venn
        
        fig = plt.figure(figsize=(14, 10))
        
        caller_svs = {}
        for caller, vcf_file in vcf_files.items():
            if os.path.exists(vcf_file):
                svs = parse_vcf(vcf_file)
                pass_svs = set()
                for sv in svs:
                    if sv['filter'] == 'PASS':
                        key = f"{sv['chrom']}:{sv['pos']//100}:{sv['type']}"
                        pass_svs.add(key)
                if pass_svs:
                    caller_svs[caller] = pass_svs
                    print(f"  {caller}: {len(pass_svs)} PASS calls")
        
        callers = list(caller_svs.keys())
        
        if len(callers) < 2:
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, 'Need at least 2 callers with PASS calls', 
                   ha='center', va='center', fontsize=14)
            ax.axis('off')
        
        elif len(callers) == 4:
            ax = fig.add_subplot(111)
            labels = {c.upper(): caller_svs[c] for c in callers}
            v = venn.venn(labels, ax=ax, fontsize=11)
            ax.set_title('4-Way Caller Overlap\n(PASS calls only)', 
                        fontsize=16, fontweight='bold')
            
            all_svs = set.union(*caller_svs.values())
            all_intersection = set.intersection(*caller_svs.values())
            
            stats_text = [
                f"Total unique SVs: {len(all_svs)}",
                f"Found by all 4 callers: {len(all_intersection)}",
                "",
                "Individual caller counts:"
            ]
            for caller in callers:
                stats_text.append(f"  {caller.upper()}: {len(caller_svs[caller])}")
            
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
            ax.text(0.02, 0.98, '\n'.join(stats_text), 
                   transform=ax.transAxes, fontsize=10, 
                   verticalalignment='top', bbox=props)
        
        elif len(callers) == 3:
            from matplotlib_venn import venn3, venn3_circles
            ax = fig.add_subplot(111)
            
            venn3(list(caller_svs.values()), 
                 set_labels=[c.upper() for c in callers], ax=ax)
            venn3_circles(list(caller_svs.values()), ax=ax, linewidth=1)
            
            ax.set_title('3-Way Caller Overlap', fontsize=16, fontweight='bold')
        
        elif len(callers) == 2:
            from matplotlib_venn import venn2, venn2_circles
            ax = fig.add_subplot(111)
            
            venn2(list(caller_svs.values()), 
                 set_labels=[c.upper() for c in callers], ax=ax)
            venn2_circles(list(caller_svs.values()), ax=ax, linewidth=1)
            
            ax.set_title('2-Way Caller Overlap', fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/9_caller_overlap.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 9_caller_overlap.png")
        return True
        
    except Exception as e:
        print(f"Error generating overlap plot: {e}")
        return False

def plot_10_complex_svs(truth_svs, enhanced_svs, output_dir):
    """Plot 10: Complex SVs analysis"""
    try:
        fig = plt.figure(figsize=(16, 10))
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2])
        
        truth_complex = [sv for sv in truth_svs if sv['is_complex'] or sv['type'] == 'CPX']
        enhanced_complex = [sv for sv in enhanced_svs if (sv['is_complex'] or sv['type'] == 'CPX') and sv['filter'] == 'PASS']
        
        # Bar comparison
        ax1 = fig.add_subplot(gs[0, 0])
        categories = ['Truth', 'Enhanced (PASS)']
        counts = [len(truth_complex), len(enhanced_complex)]
        colors = ['#3498db', '#2ecc71']
        
        bars = ax1.bar(categories, counts, color=colors, alpha=0.9, width=0.6)
        for bar, count in zip(bars, counts):
            ax1.text(bar.get_x() + bar.get_width()/2., count + 0.5,
                    str(count), ha='center', va='bottom', fontsize=12, fontweight='bold')
        
        if len(truth_complex) > 0:
            detection_rate = len(enhanced_complex) / len(truth_complex) * 100
            ax1.text(0.5, -0.15, f'Detection Rate: {detection_rate:.1f}%',
                    transform=ax1.transAxes, ha='center', fontsize=11, style='italic')
        
        ax1.set_ylabel('Number of Complex SVs', fontsize=12, fontweight='bold')
        ax1.set_title('Complex SV Detection', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='y')
        
        # Pattern distribution
        ax2 = fig.add_subplot(gs[0, 1])
        pattern_counts = defaultdict(int)
        for sv in enhanced_complex:
            if sv['complex_components']:
                pattern = ','.join(sorted(sv['complex_components']))
                pattern_counts[pattern] += 1
        
        if pattern_counts:
            patterns = list(pattern_counts.keys())[:20]
            counts = [pattern_counts[p] for p in patterns]
            
            bars = ax2.barh(range(len(patterns)), counts, color='#e74c3c', alpha=0.9)
            ax2.set_yticks(range(len(patterns)))
            ax2.set_yticklabels(patterns, fontsize=10)
            ax2.set_xlabel('Count', fontsize=12, fontweight='bold')
            ax2.grid(True, alpha=0.3, axis='x')
            
            for bar, count in zip(bars, counts):
                ax2.text(count + 0.1, bar.get_y() + bar.get_height()/2.,
                        str(count), ha='left', va='center', fontsize=10, fontweight='bold')
        else:
            ax2.text(0.5, 0.5, 'No complex patterns detected', ha='center', va='center', fontsize=12)
        
        ax2.set_title('Complex SV Patterns Detected', fontsize=14, fontweight='bold')
        
        # Summary table
        ax3 = fig.add_subplot(gs[1, :])
        ax3.axis('off')
        
        table_data = [['Pattern', 'Count', 'Example Components']]
        pattern_details = defaultdict(list)
        for sv in enhanced_complex:
            if sv['complex_components']:
                pattern = ','.join(sorted(sv['complex_components']))
                pattern_details[pattern].append(sv['size'])
        
        for pattern, sizes in sorted(pattern_details.items(), key=lambda x: len(x[1]), reverse=True)[:8]:
            count = len(sizes)
            components = pattern.split(',')
            example = ' + '.join(components[:3]) + ('...' if len(components) > 3 else '')
            table_data.append([pattern, str(count), example])
        
        if len(table_data) > 1:
            table = ax3.table(cellText=table_data, loc='center', cellLoc='left',
                             colWidths=[0.4, 0.1, 0.5])
            table.auto_set_font_size(False)
            table.set_fontsize(11)
            table.scale(1.2, 2)
            
            for i in range(len(table_data[0])):
                table[(0, i)].set_facecolor('#3498db')
                table[(0, i)].set_text_props(weight='bold', color='white')
        
        ax3.set_title('Complex SV Pattern Details', fontsize=14, fontweight='bold', pad=20)
        
        plt.suptitle('Complex SV Analysis\n(Detection and characterization of complex structural variants)', 
                     fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/10_complex_svs_analysis.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 10_complex_svs_analysis.png")
        return True
    except Exception as e:
        print(f"Error generating plot 10: {e}")
        return False

def plot_11_sv_size_distributions(truth_svs, enhanced_svs, output_dir):
    """Plot 11: SV size distributions"""
    try:
        fig, axes = plt.subplots(2, 4, figsize=(20, 12))
        axes = axes.flatten()
        
        # Enhanced PASS calls only
        enhanced_pass = [sv for sv in enhanced_svs if sv['filter'] == 'PASS']
        
        sv_types = ['DEL', 'INS', 'DUP_TANDEM', 'DUP_INTERSPERSED', 'INV', 'BND', 'CPX']
        truth_color = '#3498db'
        detected_colors = ['#e74c3c', '#2ecc71', '#f39c12', '#9b59b6', '#e67e22', '#1abc9c', '#34495e']
        
        for idx, sv_type in enumerate(sv_types):
            ax = axes[idx]
            
            if sv_type == 'DUP_TANDEM':
                truth_sizes = [sv['size'] for sv in truth_svs 
                              if sv['type'] == 'DUP' and sv.get('dup_type') == 'TANDEM' and sv['size'] > 0]
                detected_sizes = [sv['size'] for sv in enhanced_pass 
                                if sv['type'] == 'DUP' and sv.get('dup_subtype') == 'TANDEM' and sv['size'] > 0]
                plot_title = 'DUP TANDEM'
                
            elif sv_type == 'DUP_INTERSPERSED':
                truth_sizes = [sv['size'] for sv in truth_svs 
                              if sv['type'] == 'DUP' and sv.get('dup_type') == 'INTERSPERSED' and sv['size'] > 0]
                detected_sizes = [sv['size'] for sv in enhanced_pass 
                                if sv['type'] == 'DUP' and sv.get('dup_subtype') == 'INTERSPERSED' and sv['size'] > 0]
                plot_title = 'DUP INTERSPERSED'
                
            elif sv_type == 'BND':
                truth_bnds = [sv for sv in truth_svs if sv['type'] == 'BND']
                detected_bnds = [sv for sv in enhanced_pass if sv['type'] == 'BND']
                
                categories = ['Truth\n(Simulated)', 'Detected\n(PASS)']
                counts = [len(truth_bnds), len(detected_bnds)]
                
                bars = ax.bar(categories, counts, color=[truth_color, detected_colors[idx]], alpha=0.8)
                
                for bar, count in zip(bars, counts):
                    if count > 0:
                        ax.text(bar.get_x() + bar.get_width()/2., count + max(counts)*0.02,
                               str(count), ha='center', va='bottom', fontsize=12, fontweight='bold')
                
                ax.set_ylabel('Count of Breakends')
                ax.set_title('BND Count Distribution\n(BNDs are breakpoints, not size-based)', fontsize=11, fontweight='bold')
                ax.grid(True, alpha=0.3, axis='y')
                continue
                
            else:
                base_type = sv_type
                truth_sizes = [sv['size'] for sv in truth_svs 
                              if sv['type'] == base_type and sv['size'] > 0]
                detected_sizes = [sv['size'] for sv in enhanced_pass 
                                if sv['type'] == base_type and sv['size'] > 0]
                plot_title = base_type
            
            if sv_type != 'BND':
                if truth_sizes or detected_sizes:
                    all_sizes = truth_sizes + detected_sizes
                    if all_sizes and max(all_sizes) > min(all_sizes):
                        try:
                            bins = np.logspace(np.log10(max(50, min(all_sizes))), 
                                             np.log10(max(all_sizes)), 25)
                            ax.set_xscale('log')
                        except (ValueError, ZeroDivisionError):
                            bins = np.linspace(min(all_sizes), max(all_sizes), 25)
                        
                        if truth_sizes:
                            ax.hist(truth_sizes, bins=bins, alpha=0.7, 
                                   label=f'Truth (n={len(truth_sizes)})', 
                                   color=truth_color, density=True)
                        if detected_sizes:
                            ax.hist(detected_sizes, bins=bins, alpha=0.7, 
                                   label=f'Detected (n={len(detected_sizes)})', 
                                   color=detected_colors[idx], density=True)
                        
                        ax.set_xlabel('Size (bp)')
                        ax.set_ylabel('Density')
                        ax.legend(fontsize=9)
                        ax.grid(True, alpha=0.3)
                        
                        if truth_sizes and len(truth_sizes) > 1:
                            truth_median = np.median(truth_sizes)
                            ax.axvline(truth_median, color=truth_color, linestyle='--', alpha=0.8)
                        if detected_sizes and len(detected_sizes) > 1:
                            detected_median = np.median(detected_sizes)
                            ax.axvline(detected_median, color=detected_colors[idx], linestyle='--', alpha=0.8)
                    else:
                        ax.bar(['Truth', 'Detected'], [len(truth_sizes), len(detected_sizes)], 
                              color=[truth_color, detected_colors[idx]], alpha=0.7)
                        ax.set_ylabel('Count')
                else:
                    ax.text(0.5, 0.5, f'No {plot_title} data\n≥50bp', ha='center', va='center', 
                           transform=ax.transAxes, fontsize=12)
                
                ax.set_title(f'{plot_title} Size Distribution', fontsize=11, fontweight='bold')
        
        if len(axes) > len(sv_types):
            axes[-1].axis('off')
        
        plt.suptitle('SV Size Distributions: Truth vs Detected PASS Calls\n(DUP subtypes shown separately, BND shown as counts)', 
                     fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/11_sv_size_distributions.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 11_sv_size_distributions.png")
        return True
    except Exception as e:
        print(f"Error generating plot 11: {e}")
        return False

def plot_12_corrected_genotype_distribution(truth_svs, enhanced_svs, output_dir):
    """Plot 12: Corrected genotype distribution"""
    try:
        # Filter enhanced to PASS only
        enhanced_pass = [sv for sv in enhanced_svs if sv['filter'] == 'PASS']
        
        print(f"\nData Validation:")
        print(f"Truth SVs: {len(truth_svs)}")
        print(f"Enhanced total: {len(enhanced_svs)}")
        print(f"Enhanced PASS: {len(enhanced_pass)}")
        
        # Count genotypes
        truth_gt_counts = Counter([sv['gt_category'] for sv in truth_svs])
        enhanced_gt_counts = Counter([sv['gt_category'] for sv in enhanced_pass])
        
        print(f"\nTruth genotype counts: {dict(truth_gt_counts)}")
        print(f"Enhanced PASS genotype counts: {dict(enhanced_gt_counts)}")
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        gt_categories = ['Heterozygous (0/1)', 'Homozygous (1/1)', 'Reference (0/0)', 'Unknown']
        colors = ['#2ecc71', '#e74c3c', '#95a5a6', '#34495e']
        
        # Truth set pie chart
        truth_sizes = [truth_gt_counts.get(gt, 0) for gt in gt_categories]
        truth_labels = []
        truth_colors_filtered = []
        truth_sizes_filtered = []
        
        for i, (gt, count) in enumerate(zip(gt_categories, truth_sizes)):
            if count > 0:
                truth_labels.append(f'{gt}\n({count})')
                truth_colors_filtered.append(colors[i])
                truth_sizes_filtered.append(count)
        
        if truth_sizes_filtered:
            wedges1, texts1, autotexts1 = ax1.pie(truth_sizes_filtered, labels=truth_labels, 
                                                colors=truth_colors_filtered, autopct='%1.1f%%', 
                                                startangle=90, textprops={'fontsize': 11})
            for autotext in autotexts1:
                autotext.set_fontweight('bold')
                autotext.set_color('white')
        else:
            ax1.text(0.5, 0.5, 'No genotype data', ha='center', va='center', fontsize=12)
        
        ax1.set_title(f'Truth Set Genotypes\n({len(truth_svs)} total SVs)', fontsize=14, fontweight='bold')
        
        # Enhanced calls pie chart
        enhanced_sizes = [enhanced_gt_counts.get(gt, 0) for gt in gt_categories]
        enhanced_labels = []
        enhanced_colors_filtered = []
        enhanced_sizes_filtered = []
        
        for i, (gt, count) in enumerate(zip(gt_categories, enhanced_sizes)):
            if count > 0:
                enhanced_labels.append(f'{gt}\n({count})')
                enhanced_colors_filtered.append(colors[i])
                enhanced_sizes_filtered.append(count)
        
        if enhanced_sizes_filtered:
            wedges2, texts2, autotexts2 = ax2.pie(enhanced_sizes_filtered, labels=enhanced_labels,
                                                colors=enhanced_colors_filtered, autopct='%1.1f%%',
                                                startangle=90, textprops={'fontsize': 11})
            for autotext in autotexts2:
                autotext.set_fontweight('bold')
                autotext.set_color('white')
        else:
            ax2.text(0.5, 0.5, 'No genotype data', ha='center', va='center', fontsize=12)
        
        ax2.set_title(f'Enhanced PASS Calls\n({len(enhanced_pass)} PASS SVs)', fontsize=14, fontweight='bold')
        
        # Add summary statistics
        truth_total = sum(truth_sizes)
        enhanced_total = sum(enhanced_sizes)
        
        het_rate_truth = truth_gt_counts.get("Heterozygous (0/1)", 0) / max(1, truth_total) * 100
        het_rate_enhanced = enhanced_gt_counts.get("Heterozygous (0/1)", 0) / max(1, enhanced_total) * 100
        hom_rate_truth = truth_gt_counts.get("Homozygous (1/1)", 0) / max(1, truth_total) * 100
        hom_rate_enhanced = enhanced_gt_counts.get("Homozygous (1/1)", 0) / max(1, enhanced_total) * 100
        
        summary_text = (f'Truth: {truth_total} SVs | Enhanced: {enhanced_total} PASS SVs\n'
                       f'Het Rate - Truth: {het_rate_truth:.1f}% | Enhanced: {het_rate_enhanced:.1f}%\n'
                       f'Hom Rate - Truth: {hom_rate_truth:.1f}% | Enhanced: {hom_rate_enhanced:.1f}%')
        
        fig.text(0.5, 0.02, summary_text, ha='center', fontsize=11, style='italic')
        
        plt.suptitle('Corrected Genotype Distribution\n(Truth Set vs Enhanced PASS calls - Data Validated)', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.15)
        
        plt.savefig(f"{output_dir}/12_corrected_genotype_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: 12_corrected_genotype_distribution.png")
        
        # Print detailed comparison table
        print(f"\nDetailed Genotype Comparison:")
        print(f"{'Category':<20} {'Truth Count':<12} {'Truth %':<10} {'Enhanced Count':<15} {'Enhanced %':<12}")
        print("-" * 75)
        
        for gt in gt_categories:
            truth_count = truth_gt_counts.get(gt, 0)
            enhanced_count = enhanced_gt_counts.get(gt, 0)
            truth_pct = truth_count / max(1, truth_total) * 100
            enhanced_pct = enhanced_count / max(1, enhanced_total) * 100
            
            print(f"{gt:<20} {truth_count:<12} {truth_pct:<10.1f} {enhanced_count:<15} {enhanced_pct:<12.1f}")
        
        return True
    except Exception as e:
        print(f"Error generating plot 12: {e}")
        return False

def plot_7_svim_dup_subtypes_pie(output_dir):
    """Create pie chart of SVIM DUP subtypes like plot 5"""
    try:
        svim_vcf = f"{output_dir}/svim_calls.vcf"
        
        if not os.path.exists(svim_vcf):
            print("SVIM VCF not found, skipping DUP subtype plot")
            return False
        
        print("Analyzing SVIM DUP subtypes...")
        svs = parse_vcf(svim_vcf)
        svim_pass_dups = [sv for sv in svs if sv['type'] == 'DUP' and sv['filter'] == 'PASS']
        
        print(f"Found {len(svim_pass_dups)} PASS DUPs in SVIM")
        
        if not svim_pass_dups:
            print("No SVIM PASS DUPs found for subtype analysis")
            return False
        
        # Count subtypes
        subtype_counts = Counter()
        for dup in svim_pass_dups:
            subtype = dup.get('dup_subtype', 'Unknown')
            if subtype:
                subtype_counts[subtype] += 1
            else:
                subtype_counts['Unknown'] += 1
        
        print(f"SVIM DUP subtypes found: {dict(subtype_counts)}")
        
        # Create pie chart
        fig, ax = plt.subplots(figsize=(10, 8))
        
        if subtype_counts and sum(subtype_counts.values()) > 0:
            labels = []
            sizes = []
            colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6']
            
            for i, (subtype, count) in enumerate(subtype_counts.most_common()):
                labels.append(f'{subtype}\n({count})')
                sizes.append(count)
            
            wedges, texts, autotexts = ax.pie(sizes, labels=labels, 
                                             colors=colors[:len(sizes)],
                                             autopct='%1.1f%%', startangle=90,
                                             textprops={'fontsize': 11})
            
            for autotext in autotexts:
                autotext.set_fontweight('bold')
                autotext.set_color('white')
                autotext.set_fontsize(12)
        else:
            ax.text(0.5, 0.5, 'No DUP subtypes\nfound in SVIM', 
                   ha='center', va='center', fontsize=14)
        
        ax.set_title(f'SVIM Duplication Subtypes\n(Total PASS DUPs: {len(svim_pass_dups)})', 
                     fontsize=16, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/figures/svim_dup_subtypes_pie.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Generated: svim_dup_subtypes_pie.png")
        return True
        
    except Exception as e:
        print(f"Error generating SVIM DUP subtypes plot: {e}")
        return False

def main():
    """Main execution function"""
    if len(sys.argv) < 2:
        print("Usage: python3 complete_sv_visualization.py <output_dir> [eval_dir]")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    eval_dir = sys.argv[2] if len(sys.argv) > 2 else f"{output_dir}/eval_results"
    
    # Create figures directory
    figures_dir = f"{output_dir}/figures"
    ensure_dir(figures_dir)
    
    print("\n" + "="*70)
    print("COMPLETE SV ANALYSIS VISUALIZATION SUITE")
    print("="*70)
    
    # Find VCF files
    vcf_files = find_vcf_files(output_dir)
    
    if not vcf_files:
        print("No VCF files found!")
        return
    
    # Parse VCFs
    truth_svs = []
    enhanced_svs = []
    
    if 'truth' in vcf_files:
        print("\nParsing truth VCF...")
        truth_svs = parse_vcf(vcf_files['truth'])
        truth_in_repeats = len([sv for sv in truth_svs if sv['in_repeat']])
        print(f"  Found {len(truth_svs)} truth SVs")
        if truth_svs:
            print(f"  {truth_in_repeats} in repeat regions ({truth_in_repeats/len(truth_svs)*100:.1f}%)")
    
    if 'enhanced' in vcf_files:
        print("\nParsing enhanced VCF...")
        enhanced_svs = parse_vcf(vcf_files['enhanced'])
        enhanced_pass = [sv for sv in enhanced_svs if sv['filter'] == 'PASS']
        enhanced_in_repeats = len([sv for sv in enhanced_pass if sv['in_repeat']])
        print(f"  Found {len(enhanced_svs)} enhanced SVs ({len(enhanced_pass)} PASS)")
        if enhanced_pass:
            print(f"  {enhanced_in_repeats} PASS calls in repeat regions ({enhanced_in_repeats/len(enhanced_pass)*100:.1f}%)")
    
    # Load evaluation and runtime results
    print("\nLoading results...")
    results = load_truvari_results(eval_dir)
    runtime_data, memory_data = load_runtime_data(output_dir)
    
    print(f"Found Truvari results for: {list(results.keys())}")
    print(f"Found runtime data for: {list(runtime_data.keys())}")
    print(f"Found memory data for: {list(memory_data.keys())}")
    
    print("\nGenerating all plots...")
    print("-" * 50)
    
    success_count = 0
    
    # Generate all plots
    if plot_1_performance_metrics(results, figures_dir):
        success_count += 1
    
    if plot_2_runtime_comparison(runtime_data, memory_data, figures_dir):
        success_count += 1
    
    if plot_3_tp_fp_fn_stacked(results, figures_dir):
        success_count += 1

    if plot_7_svim_dup_subtypes_pie(output_dir):
        success_count += 1
    
    if truth_svs and enhanced_svs:
        if plot_4_sv_types_comparison(truth_svs, enhanced_svs, figures_dir):
            success_count += 1
        if plot_5_duplication_pie_charts(truth_svs, enhanced_svs, figures_dir):
            success_count += 1
        if plot_6_repeat_svs_overview(truth_svs, enhanced_svs, figures_dir):
            success_count += 1
        if plot_8_repeat_types_distribution(truth_svs, enhanced_svs, figures_dir):
            success_count += 1
        if plot_10_complex_svs(truth_svs, enhanced_svs, figures_dir):
            success_count += 1
        if plot_11_sv_size_distributions(truth_svs, enhanced_svs, figures_dir):
            success_count += 1
        if plot_12_corrected_genotype_distribution(truth_svs, enhanced_svs, figures_dir):
            success_count += 1
    
    # VCF files for overlap plot
    vcf_files_for_venn = {
        'enhanced': vcf_files.get('enhanced', ''),
        'sniffles2': vcf_files.get('sniffles2', ''),
        'cutesv': vcf_files.get('cutesv', ''),
        'svim': vcf_files.get('svim', '')
    }
    
    # Filter to existing files
    vcf_files_for_venn = {k: v for k, v in vcf_files_for_venn.items() if v and os.path.exists(v)}
    
    if len(vcf_files_for_venn) >= 2:
        if plot_9_caller_overlap_venn(vcf_files_for_venn, figures_dir):
            success_count += 1
    
    print("\n" + "="*70)
    print("VISUALIZATION COMPLETE")
    print("="*70)
    print(f"\nSuccessfully generated {success_count} plots")
    print(f"All plots saved to: {figures_dir}/")
    
    print("\nGenerated plots:")
    print("  1. Performance Metrics - Precision, Recall, F1 from Truvari evaluation")
    print("  2. Runtime Comparison - Time and memory usage")
    print("  3. TP/FP/FN Analysis - True/False Positives and False Negatives")
    print("  4. SV Types Comparison - Simulated vs Detected by type")
    print("  5. Duplication Types - Tandem vs Interspersed duplications")
    print("  6. Repeat Regions Overview - SVs in repeat regions")
    print("  7. SVIM DUP subtypes plot")
    print("  8. Repeat Type Distribution - Specific repeat type analysis")
    print("  9. Caller Overlap - Venn diagram of caller agreement")
    print("  10. Complex SVs - Complex structural variant analysis")
    print("  11. Size Distributions - SV size distributions by type")
    print("  12. Corrected Genotype Distribution - Fixed genotype analysis")

if __name__ == "__main__":
    main()

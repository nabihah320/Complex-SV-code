#!/usr/bin/env python3
"""
Extract Truvari results from JSON summary
"""
import json
import sys

def extract_results(summary_json, caller, csv_file):
    """Extract and format Truvari results"""
    with open(summary_json) as f:
        data = json.load(f)
        
    precision = float(data.get('precision', 0)) if data.get('precision') is not None else 0
    recall = float(data.get('recall', 0)) if data.get('recall') is not None else 0
    f1 = float(data.get('f1', 0)) if data.get('f1') is not None else 0
    tp = int(data.get('TP-base', data.get('TP', 0)))
    fp = int(data.get('FP', 0))
    fn = int(data.get('FN', 0))
    
    print(f'  Precision: {precision:.3f}')
    print(f'  Recall:    {recall:.3f}')
    print(f'  F1 Score:  {f1:.3f}')
    print(f'  TP: {tp}, FP: {fp}, FN: {fn}')
    
    with open(csv_file, 'a') as csv_out:
        csv_out.write(f'{caller},{precision:.3f},{recall:.3f},{f1:.3f},{tp},{fp},{fn}\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 extract_truvari_results.py <summary.json> <caller_name> <output.csv>")
        sys.exit(1)
    
    extract_results(sys.argv[1], sys.argv[2], sys.argv[3])

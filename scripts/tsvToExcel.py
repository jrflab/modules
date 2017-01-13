#!/usr/bin/env python
"""
Write tsv to excel sheet
"""
import argparse
import pandas as pd
import os
from openpyxl import load_workbook


def write_to_excel(tsv_file, excel_file, sheet_name, column_names, delimiter, overwrite):
    df = pd.read_csv(tsv_file, sep=delimiter)
    if not overwrite and os.path.isfile(excel_file):
        book = load_workbook(excel_file)
        writer = pd.ExcelWriter(excel_file, engine='openpyxl')
        writer.book = book
        writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
    else:
        writer = pd.ExcelWriter(excel_file)
    if column_names:
        df.to_excel(writer, sheet_name, columns=column_names, index=False)
    else:
        df.to_excel(writer, sheet_name, index=False)
    if not overwrite and os.path.isfile(excel_file):
        writer.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tsv_file", type=str, help="TSV")
    parser.add_argument("excel_file", type=str, help="Excel output")
    parser.add_argument("sheet_name", type=str, help="Sheet name")
    parser.add_argument("-c", "--column_names", type=str, default=None, help="Which columns to write (comma separated)")
    parser.add_argument("-d", "--delimiter", type=str, default="\t", help="Set delimiter")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing excel")
    args = parser.parse_args()
    if args.column_names:
        column_names = args.column_names.split(",")
    else:
        column_names = None
    sheet_name = (args.sheet_name[:25] + '..') if len(args.sheet_name) > 25 else args.sheet_name
    write_to_excel(args.tsv_file, args.excel_file, sheet_name, column_names, args.delimiter, args.overwrite)

if __name__ == "__main__":
    main()

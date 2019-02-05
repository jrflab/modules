#!/usr/bin/env python
import csv
from xlsxwriter.workbook import Workbook

tsv_file = 'summary/tsv/cravat_summary.tsv'
xlsx_file = 'summary/cravat_summary.xlsx'

workbook = Workbook(xlsx_file)
worksheet = workbook.add_worksheet()

tsv_reader = csv.reader(open(tsv_file, 'rb'), delimiter='\t')

for row, data in enumerate(tsv_reader):
    worksheet.write_row(row, 0, data)

workbook.close()

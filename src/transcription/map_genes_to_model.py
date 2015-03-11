import argparse
import csv

parser = argparse.ArgumentParser(description='Map genes to model.')
parser.add_argument('input', metavar='I', type=str, nargs='+',
                    help='Input files')

args = parser.parse_args()

for i in args.input:
  with open(i) as f:
    reader = csv.reader(f)
    
    lengthcol = 20
    for row in reader:
      print(', '.join([row[0], row[lengthcol]]))
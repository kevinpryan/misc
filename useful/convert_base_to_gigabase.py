#!/usr/bin/python
import csv
filename = "count_bases.txt"

# read in file as dictionary - every second line is filename, bases
with open(filename) as f:
    counts = {int(name.strip()):score.strip() for score, name in zip(f, f)}

# I want Sample to be key and bases to be value
inv_map = {v: k for k, v in counts.items()}

gigabases_dict = {}
gigabases_list = []

# create list of dictionaries containing gigabases for each sample/file
for x in inv_map:
    gigabases_dict["Sample"] = x
    gigabases_dict["Gigabases"] = str(round(inv_map[x]/1e9, 2))
    gigabases_list.append(gigabases_dict)
    gigabases_dict = {}

# write to file
csv_columns = ['Sample','Gigabases']
with open('gigabases.csv', 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
    writer.writeheader()
    writer.writerows(gigabases_list)

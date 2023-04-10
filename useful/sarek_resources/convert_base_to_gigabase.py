#!/usr/bin/python
import csv
filename = "count_bases.txt"

my_dict = {}
with open(filename) as f:
    lines = f.readlines()
    for i in range(0, len(lines), 2):
        try:
            sample_name = lines[i].strip()
            number_of_bases = int(lines[i+1])
            my_dict[sample_name] = number_of_bases
        except IndexError:
            print(f"Skipping last line: {lines[i].strip()}")

gigabases_dict = {}
gigabases_list = []

# create list of dictionaries containing gigabases for each sample/file
for x in my_dict:
    gigabases_dict["Sample"] = x
    gigabases_dict["Gigabases"] = str(round(my_dict[x]/1e9, 3))
    gigabases_list.append(gigabases_dict)
    gigabases_dict = {}

# write to file
csv_columns = ['Sample','Gigabases']
with open('gigabases.csv', 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
    writer.writeheader()
    writer.writerows(gigabases_list)

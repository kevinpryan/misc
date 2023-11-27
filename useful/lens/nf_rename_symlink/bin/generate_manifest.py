#!/usr/bin/env python3

import os
import sys
import errno
import argparse
import csv
def parse_args(args=None):
    Description = "Generate manifest for LENS"
    Epilog = "Example usage: python generate_manifest.py <SAMPLESHEET> <FILE_OUT> <SEQUENCING_METHOD> <-d DATASET>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("SAMPLESHEET", help="Input samplesheet file, output from rename_files_general.R")
    parser.add_argument("FILE_OUT", help="Output file.")
    parser.add_argument(
        "-d",
        "--dataset",
        dest="DATASET",
        default="caf_lens",
        help="dataset name",
    )
    return parser.parse_args(args)

def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)

def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception

def samplesheet_to_manifest(
    samplesheet_file,
    fileout,
    dataset,
):
    sample_mapping_dict = {}
    with open(samplesheet_file, "r", encoding="utf-8-sig") as fin:
        HEADER_MANIFEST = ["Patient_Name","Run_Name","Dataset","File_Prefix","Sequencing_Method","Normal"]
        HEADER_SAMPLESHEET = ["file_original", "patient", "abnormal_normal", "file_lens", "Sequencing_Method"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        # loop through each line, get manifest header for each
        # Dataset is the same for all files
        # Sequencing_Method is same as that in samplesheet
        for line in fin:
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]
                if len(lspl) != len(HEADER_SAMPLESHEET):
                    print_error("samplesheet not the required length", "Line", line)

                file_original, patient, abnormal_normal, file_lens, Sequencing_Method = lspl[: len(HEADER_SAMPLESHEET)]
            if int(patient) < 10:
                Patient_Name = "Pt-0" + str(patient)
            elif int(patient) >= 10:
                Patient_Name = "Pt-" + str(patient)
            else:
                print_error("Invalid patient name, must be an integer", "Line", line)
            Run_Name = file_lens.split("_")[0]
            File_Prefix = Run_Name
            if abnormal_normal == "a":
                Normal = "FALSE"
            elif abnormal_normal == "n":
                Normal = "TRUE"
            else:
                print_error("abnormal_normal must be either a (abnormal) or n (normal)", "Line", line)
            # put sample_info together
            sample_info = [Patient_Name, Run_Name, dataset, File_Prefix, Sequencing_Method, Normal]
            # add sample info to sample mapping dict
            if Run_Name not in sample_mapping_dict:
                sample_mapping_dict[Run_Name] = sample_info
    # write manifest to file
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(fileout)
        # make outfile
        make_dir(out_dir)
        with open(fileout, "w", newline = '') as fout:
            # not sure about this line
            writer = csv.writer(fout)
            writer.writerow(HEADER_MANIFEST)

            for key, value in sample_mapping_dict.items():
                writer.writerow(value)

def main(args=None):
    args = parse_args(args)
    samplesheet_to_manifest(
        samplesheet_file=args.SAMPLESHEET,
        fileout=args.FILE_OUT,
        dataset = args.DATASET
        )
if __name__ == "__main__":
    sys.exit(main())

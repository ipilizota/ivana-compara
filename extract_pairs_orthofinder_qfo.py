#!/usr/bin/env python3
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Write OrthoFinder pairwise predictions to a file compatible with Orthology Benchmark Service.

Typical usage examples::

    $ python extract_orthologues_qfo_orthofinder.py \
    --predictions /path/to/OrthoFinder/Results_XXXXX/Orthologues \
    --out_file /path/to/output/file.txt

"""

import argparse
import os
from typing import List
import warnings


def extract_orthologous_pairs(input_file: str) -> List[tuple]:
    """Extracts pairs of putative orthologues from a single OrthoFinder output file.

    Args:
        input_file: Path to the orthologues spreadsheet produced by OrthoFinder.

    Returns:
        A list of tuples where each tuple contains UniProtKB accessions for a pair of putative orthologues.

    Warns:
        UserWarning: If `input_file` cannot be opened or if a line in `input_file` is corrupted.

    """
    orthologues = []

    try:
        with open(input_file) as file_handler:
            for line in file_handler:
                if "Orthogroup" in line:
                    continue

                try:
                    temp = line.replace("\n", "").split("\t")
                    reference_genes = temp[1].replace(" ", "").split(",")
                    target_genes = temp[2].replace(" ", "").split(",")
                except IndexError:
                    warnings.warn(f"There is an issue with file '{input_file}' in line '{line}'.")
                    continue

                for reference in reference_genes:
                    reference_id = reference.split("|")[1]
                    for target in target_genes:
                        orthologues.append((reference_id, target.split("|")[1]))

    except EnvironmentError:
        warnings.warn(f"Could not open file '{input_file}'.")

    return orthologues


def write_orthologous_pairs(orthologous_pairs: List[tuple], output_file: str) -> None:
    """Appends pairs of orthologues to the end of a specified file.

    Args:
        orthologous_pairs: A list of tuples representing pairs of putative orthologues.
        output_file: Path to the output file containing pairs of putative orthologues.

    """
    with open(output_file, "a") as file_handler:
        for pair in orthologous_pairs:
            file_handler.write(pair[0] + "\t" + pair[1] + "\n")


def process_orthofinder_predictions(predictions: str, out_file: str) -> None:
    """Extracts OrthoFinder pairwise predictions and writes them to a single output file.

    Predictions are made on a QfO dataset and the code relies on the formatting of QfO fasta files.
    Output file is compatible with Orthology Benchmark Service.

    Args:
        predictions: Path to OrthoFinder's 'Orthologues' directory.
        output_file: Path to the output file containing pairs of putative orthologues.

    """
    predictions_subdirs = [f.path for f in os.scandir(predictions) if f.is_dir()]

    for dir in predictions_subdirs:
        files = [f.path for f in os.scandir(dir) if f.is_file() and f.path.endswith(".tsv")]
        for file in files:
            orthologs = extract_orthologous_pairs(file)
            write_orthologous_pairs(orthologs, out_file)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", required=True, type=str,
                        help="Path to OrthoFinder's 'Orthologues' directory")
    parser.add_argument("--out_file", required=True, type=str, help="Output file")

    args = parser.parse_args()

    process_orthofinder_predictions(args.predictions, args.out_file)

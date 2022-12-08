#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os
from Bio import SeqIO
import json
from abc import ABCMeta, abstractmethod
from script.settings import logger
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from pyfaidx import Fasta


class AMPBase(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def from_string(self): pass

    @abstractmethod
    def from_args(self): pass

    @abstractmethod
    def validate_inputs(self): pass

    @abstractmethod
    def run(self): pass

    @abstractmethod
    def create_databases(self): pass

    @abstractmethod
    def run_blast(self): pass

    @abstractmethod
    def filter_process(self): pass

    @abstractmethod
    def output(self): pass


class BaseModel(object):
    def extract_nth_bar(self, bunch_str, n):
        start = n + 3
        end = n + 4
        subset_split_str = bunch_str.split('|')[start: end]
        temporary_str = "|".join(subset_split_str)

        result = temporary_str[temporary_str.find(':') + 2:]
        result = result.rstrip()

        if result.isdigit():
            return int(result)

        else:
            try:
                return float(result)
            except ValueError:
                return result

    def extract_nth_hash(self, bunch_str, n):
        if "#" not in bunch_str:
            return 0
        else:
            arr = bunch_str.split("#")
            if n >= len(arr):
                return ""
            else:
                if n == 1 or n == 2:
                    return int(arr[n])

                elif n == 3:
                    if int(arr[n]) == 1:
                        return "+"
                    else:
                        return "-"
                else:
                    return arr[n]

    def find_num_dash(self, subject, index):
        dash_count = 0
        output = []

        for i in range(len(subject)):
            if subject[i] == '-':
                dash_count += 1
            else:
                output.append(subject[i])
            if len(output) == index:
                break

        return dash_count

    def get_submitted_protein_sequence(self, seq_filepath):
        submitted_proteins_dict = {}

        if os.stat(seq_filepath).st_size != 0:
            for record in SeqIO.parse(seq_filepath, 'fasta'):
                submitted_proteins_dict[record.id] = str(record.seq)

        return submitted_proteins_dict

    def get_orf_dna_sequence(self, input_file, input_type):
        filename = os.path.basename(input_file)
        predicted_genes_dict = {}

        if input_type in ["contig"]:
            contig_filepath = os.path.join(self.working_directory,
                                           filename + ".temp.contigToORF.fsa")
            if os.stat(contig_filepath).st_size != 0:
                for record in SeqIO.parse(contig_filepath, 'fasta'):
                    predicted_genes_dict[record.id] = str(record.seq)

        elif input_type in ["read"]:
            read_filepath = os.path.join(self.working_directory,
                                         filename + ".temp.read.fsa")
            if os.stat(read_filepath).st_size != 0:
                for record in SeqIO.parse(read_filepath, 'fasta'):
                    predicted_genes_dict[record.id] = str(record.seq)
        else:
            raise ValueError("input_type invalid \
                    (must be 'contig' or 'read'): {}".format(input_type))

        pjson = json.dumps(predicted_genes_dict)

        predicted_filepath = os.path.join(self.working_directory,
                                          filename + ".temp.predictedGenes.json")
        with open(predicted_filepath, 'w') as wf:
            wf.write(pjson)

        return predicted_genes_dict

    def get_orf_protein_sequence(self, input_file, input_type):

        filename = os.path.basename(input_file)
        predicted_genes_dict = {}

        if input_type in ["contig"]:
            contig_filepath = os.path.join(self.working_directory,
                                           filename + ".temp.contig.fsa")
            if os.stat(contig_filepath).st_size != 0:
                for record in SeqIO.parse(contig_filepath, 'fasta'):
                    predicted_genes_dict[record.id] = str(record.seq)

        elif input_type in ["read"]:
            read_filepath = os.path.join(self.working_directory,
                                         filename + ".temp.read.fsa")
            if os.stat(read_filepath).st_size != 0:
                for record in SeqIO.parse(read_filepath, 'fasta'):
                    predicted_genes_dict[record.id] = str(record.seq)
        else:
            raise ValueError("input_type invalid \
                    (must be 'contig' or 'read'): {}".format(input_type))

        # write json for all predicted file
        pjson = json.dumps(predicted_genes_dict)

        predicted_filepath = os.path.join(self.working_directory,
                                          filename + ".temp.predictedGenes.protein.json")
        with open(predicted_filepath, 'w') as wf:
            wf.write(pjson)

        return predicted_genes_dict

    def results(self, blast_results, query_id, perfect, strict, loose, exclude_nudge=True):
        nudged = False
        if len(perfect) == 0 and len(strict) == 0 and len(loose) > 0:
            if exclude_nudge is True:
                nudged, loose = self.nudge_loose_to_strict(loose)

            if nudged is True and self.loose is False:
                blast_results[query_id] = loose
            elif self.loose is True:
                blast_results[query_id] = loose

        elif len(perfect) == 0:
            if len(strict) > 0:
                # TODO:: add try catch here
                try:
                    nudged, strict = self.nudge_strict_to_perfect(strict)
                    blast_results[query_id] = strict
                except Exception as e:
                    logger.error(e)
        else:
            if len(perfect) > 0:
                blast_results[query_id] = perfect

        return blast_results

    def nudge_strict_to_perfect(self, strict):
        nudged = False
        for s in strict:
            if int(strict[s]["perc_identity"]) == 100 and strict[s]["type_match"] == "Strict" and strict[s][
                "model_type_id"] in [40292] and self.input_type == 'contig':
                reference = strict[s]["sequence_from_broadstreet"]
                query = strict[s]["orf_prot_sequence"]
                orf_dna_sequence_ = ""
                orf_end_ = ""
                orf_start_ = ""
                _partial_protein = ""
                partial_protein = ""
                if len(query) < len(reference) and query in reference:
                    length_nucleotides = (len(reference) - len(strict[s]["match"])) * 3

                    partial_bases, orf_start_, orf_end_ = self.get_part_sequence(
                        self.input_sequence, strict[s]["orf_from"],
                        strict[s]["orf_start"], strict[s]["orf_end"],
                        length_nucleotides, strict[s]["orf_strand"],
                        strict[s]["ARO_name"]
                    )

                    if len(partial_bases) % 3 == 0:
                        if strict[s]["orf_strand"] == "-":
                            strict[s]["orf_end_possible"] = orf_end_ + 1
                            strict[s]["orf_start_possible"] = strict[s]["orf_start"]
                            orf_dna_sequence_ = str(Seq(partial_bases, generic_dna).reverse_complement()) + strict[s][
                                "orf_dna_sequence"]
                            partial_protein = str(
                                Seq(partial_bases, generic_dna).reverse_complement().translate(table=11))
                        else:
                            strict[s]["orf_start_possible"] = orf_start_
                            strict[s]["orf_end_possible"] = int(strict[s]["orf_end"]) + 1
                            orf_dna_sequence_ = partial_bases + strict[s]["orf_dna_sequence"]
                            partial_protein = str(Seq(partial_bases, generic_dna).translate(table=11)).strip("*")
                        if len(partial_protein) > 0:
                            _partial_protein = partial_protein[0]
                            if partial_protein[0] in ["L", "M", "I", "V"]:
                                _partial_protein = "M" + partial_protein[1:]

                        combine = _partial_protein + strict[s]["match"]

                        if combine == strict[s]["sequence_from_broadstreet"]:
                            nudged = True
                            strict[s]["type_match"] = "Perfect"
                            strict[s]["nudged"] = nudged
                            strict[s]["note"] = "possible complete gene, missing n-terminus"
                            strict[s]["missed_bases_by_prodigal"] = partial_bases
                            strict[s]["orf_dna_sequence_possible"] = orf_dna_sequence_
                            strict[s]["orf_prot_sequence_possible"] = combine

                    else:
                        pass

                elif len(query) > len(reference) and reference in query:
                    if strict[s]["partial"] == "0":
                        start_codon = strict[s]["dna_sequence_from_broadstreet"][:3]
                        stop_codon = strict[s]["dna_sequence_from_broadstreet"][-3:]
                        if start_codon in ["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"] and stop_codon in ["TAA",
                                                                                                               "TAG",
                                                                                                               "TGA"]:
                            strict[s]["note"] = "possible complete gene, reference contained within open reading frame"
                            strict[s]["orf_prot_sequence_possible"] = ""

                            if strict[s]["orf_strand"] == "-":
                                strict[s]["orf_start_possible"] = strict[s]["orf_start"]
                                strict[s]["orf_end_possible"] = int(strict[s]["orf_start"]) + len(
                                    strict[s]["dna_sequence_from_broadstreet"]) - 1
                            else:
                                strict[s]["orf_start_possible"] = int(strict[s]["orf_end"]) - len(
                                    strict[s]["dna_sequence_from_broadstreet"]) + 1
                                strict[s]["orf_end_possible"] = int(strict[s]["orf_end"])

                            # pull nucleotides from query or submitted sequence
                            partial_bases, orf_start_, orf_end_ = self.get_part_sequence(
                                self.input_sequence, strict[s]["orf_from"],
                                strict[s]["orf_start_possible"], strict[s]["orf_end_possible"],
                                0, strict[s]["orf_strand"],
                                strict[s]["ARO_name"]
                            )

                            if strict[s]["orf_strand"] == "-":
                                strict[s]["orf_dna_sequence_possible"] = str(
                                    Seq(partial_bases, generic_dna).reverse_complement())
                            else:
                                strict[s]["orf_dna_sequence_possible"] = partial_bases

                            if len(strict[s]["orf_dna_sequence_possible"]) % 3 == 0:
                                orf_prot_sequence_possible = str(
                                    Seq(strict[s]["orf_dna_sequence_possible"], generic_dna).translate(table=11)).strip(
                                    "*")
                                strict[s]["orf_prot_sequence_possible"] = orf_prot_sequence_possible
                                if orf_prot_sequence_possible == strict[s]["sequence_from_broadstreet"]:
                                    strict[s]["type_match"] = "Perfect"
                                    nudged = True
                            else:
                                logger.warning(
                                    "incorrect open reading frame for coordinate: {}-{} on strand {} for {}".format
                                    (strict[s]["orf_start_possible"], strict[s]["orf_end_possible"],
                                     strict[s]["orf_strand"], strict[s]["orf_from"]))

                            strict[s]["nudged"] = nudged

                        else:
                            pass
                    else:
                        pass
            elif 95 <= int(strict[s]["perc_identity"]) < 100 and strict[s]["type_match"] == "Strict" and strict[s][
                "model_type_id"] in [40292] and self.input_type == 'contig':
                pass

        return nudged, strict

    def get_part_sequence(self, fasta_file, header, start, stop, nterminus, strand, name):
        header = header[:header.rfind("_")]
        genes = False
        try:
            genes = Fasta(fasta_file, sequence_always_upper=False, read_long_names=False, one_based_attributes=True)
        except Exception as e:
            logger.error(e)
        if genes:
            if strand == "-":
                _start = stop + 1
                _stop = stop + nterminus
                if nterminus == 0:
                    return str(genes.get_spliced_seq(header, [[start, stop]])), start, stop
                else:
                    return str(genes.get_spliced_seq(header, [[_start, _stop]])), _start, _stop
            elif strand == "+":
                _start = start - nterminus
                _stop = start - 1
                if _start <= 0:
                    _start = 1
                if _stop <= 0:
                    _stop = 1
                if nterminus == 0:
                    return str(genes.get_spliced_seq(header, [[start, stop]])), start, stop
                else:
                    return str(genes.get_spliced_seq(header, [[_start, _stop]])), _start, _stop

    def nudge_loose_to_strict(self, loose):
        nudged = False
        for i in loose:
            if 95 <= int(loose[i]["perc_identity"]) <= 100:
                loose[i]["type_match"] = "Strict"
                loose[i]["nudged"] = True
                nudged = True
                loose[i]["note"] = "loose hit with at least 95 percent identity pushed strict"

        return nudged, loose




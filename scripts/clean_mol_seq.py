import re
import pickle
import compress_json
import molseq

class CleanMolSeq(object):
    ID_RE = re.compile(">\s*(.+)")
    GAP_RE = re.compile("[-\s]+")
    WS_RE = re.compile("[\s]+")

    @staticmethod
    def stream_fasta(in_stream, remove_gaps=False):
        seq_id = None
        seq = []
        for line in in_stream:
            line = line.strip()
            if line:
                if line.startswith(">"):
                    if seq_id:
                        yield molseq.MolSeq(seq_id, "".join(seq))
                    seq_id = CleanMolSeq.ID_RE.search(line).group(1)
                    seq = []
                else:
                    if remove_gaps:
                        line = re.sub(CleanMolSeq.GAP_RE, "", line)
                    else:
                        line = re.sub(CleanMolSeq.WS_RE, "", line)
                    seq.append(line)
        if seq_id:
            yield molseq.MolSeq(seq_id, "".join(seq))

    @staticmethod
    def read_protein_fasta_file(infile, outfile, min_length):
        genomes = set()
        ignored_length = 0
        ignored_irr = 0
        ignored_name = 0
        total = 0
        f0 = open(infile)
        f1 = open(outfile, 'w')
        for mol_seq in CleanMolSeq.stream_fasta(f0, True):
            total += 1
            length = mol_seq.get_length()
            if length >= min_length:
                # Only need to call 'count_irregular_chars_aa' when dealing with protein fasta
                irr = mol_seq.count_irregular_chars_aa()
                if irr < 1:
                    name = mol_seq.get_seq_id()
                    s = name.split('|')
                    if len(s) > 2:
                        genome_acc = s[1]
                        protein_acc = s[2]
                        f1.write(genome_acc + '\t' + protein_acc + '\n')
                        genomes.add(genome_acc)
                    else:
                        ignored_name += 1
                else:
                    ignored_irr += 1
            else:
                ignored_length += 1
        f0.close()
        f1.close()
        print('Protein: Min Length               : ' + str(min_length))
        print('Protein: Total                    : ' + str(total))
        print('Protein: Ignored Length           : ' + str(ignored_length))
        print('Protein: Irreg Chars              : ' + str(ignored_irr))
        print('Protein: Ignored Ill-Formated Name: ' + str(ignored_name))
        print('Protein: Returned                 : ' + str(len(genomes)))

        return genomes

    @staticmethod
    def extract_from_protein_fasta_file(infile, outfile, keep_proteins_genome_acc):
        f0 = open(infile)
        f1 = open(outfile, 'w')
        total = 0
        for mol_seq in CleanMolSeq.stream_fasta(f0, True):
            name = mol_seq.get_seq_id()
            s = name.split('|')
            genome_acc = s[1]
            if genome_acc in keep_proteins_genome_acc:
                total += 1
                f1.write(mol_seq.to_fasta_wrapped(80))
                f1.write('\n')
        f0.close()
        f1.close()
        print('Protein: Sequences stored         : ' + str(total))

    @staticmethod
    def clean_mol_seqs(infile, outfile, min_length, min_ratio, verbose, qc_storage_outfile, qc_storage_infile = None, protein_file=None, protein_outfile=None,
                       protein_outfile_fasta=None, protein_min_length=0):
        genomes = None
        if protein_file:
            genomes = CleanMolSeq.read_protein_fasta_file(protein_file, protein_outfile, protein_min_length)
        f0 = open(infile)
        f1 = open(outfile, 'w')
        total = 0
        ignored_irr_chars = 0
        ignored_length = 0
        ignored_name = 0
        ignored_no_protein = 0
        ignored_id = 0
        kept = 0
        keep_proteins_genome_acc = set()
        keep_genome_acc = set()
        irreg_count = dict()

        if qc_storage_infile != None:
            storage = compress_json.load(qc_storage_infile)
        else:
            storage = {}

        for mol_seq in CleanMolSeq.stream_fasta(f0, True):
            name_lwr = mol_seq.get_seq_id().lower()
            acc = mol_seq.get_seq_id().split('|')[1]
            if storage.get(acc):
                continue
            seq = mol_seq.get_seq()
            reg = mol_seq.count_regular_chars_na()
            length = mol_seq.get_length()
            r = reg / length
            storage[acc] = {}
            storage[acc]['sequence'] = seq
            storage[acc]['length'] = length
            storage[acc]['ambiguity'] = 1 - r
            total += 1
            if len(name_lwr.split("|")) == 6:
                if '|severe_acute_respiratory_syndrome_related_coronavirus' in name_lwr or '2019_ncov' in name_lwr or 'hcov_19' in name_lwr or 'sars_cov_2' in name_lwr or 'sars_cov2' in name_lwr:
                    if length >= min_length:
                        if r >= min_ratio:
                            storage[acc]['QC'] = 'Pass'
                            keep_genome_acc.add(acc)
                            kept += 1
                            f1.write(mol_seq.to_fasta_wrapped(80))
                            f1.write('\n')
                        else:
                            storage[acc]['QC'] = 'Fail: Ambiguity'
                            ignored_irr_chars += 1
                    else:
                        storage[acc]['QC'] = 'Fail: Length'
                        ignored_length += 1
                else:
                    storage[acc]['QC'] = 'Fail: Name'
                    ignored_name += 1
            else:
                storage[acc]['QC'] = 'Fail: Metadata'
                ignored_id += 1

        f0.close()
        f1.close()
        if genomes:
            CleanMolSeq.extract_from_protein_fasta_file(protein_file, protein_outfile_fasta, keep_proteins_genome_acc)
        
        print('Stringency ratio             : ' + str(min_ratio))
        print('Total Inputed Sequences      : ' + str(total))
        print('Ignored, Missing Metadata    : ' + str(ignored_id))
        print('Ignored, Wrong Name          : ' + str(ignored_name))
        print('Ignored, Bad Length          : ' + str(ignored_length))
        print('Ignored, Irreg Chars         : ' + str(ignored_irr_chars))
        print('Sequences Kept               : ' + str(kept))
        print('\n')

        compress_json.dump(storage, qc_storage_outfile)

        return keep_genome_acc


if __name__ == "__main__":
     CleanMolSeq.clean_mol_seqs(
       'SARS2_Genome_Jan_2022.fasta',
       'SARS2_Genome_Jan_2022_clean.fasta',
       29400,
       0.9999,
       True,
       'SARS2_S_2-9-22.fasta',
       'SARS2_S_2-9-22_clean.txt',
       'SARS2_S_2-9-22_clean.fasta',
       1250)


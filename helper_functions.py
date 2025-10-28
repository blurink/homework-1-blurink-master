from typing import Tuple, Generator, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def codons(seq: str) -> Generator[str, None, None]:
    """Walk along the string, three nucleotides at a time. Cut off excess."""
    for i in range(0, len(seq) - 2, 3):
        yield seq[i:i + 3]


def extract_gt_orfs(record, start_codons, stop_codons, validate_cds=True, verbose=False):
    """Extract the ground truth ORFs as indicated by the NCBI annotator in the
    gene coding regions (CDS regins) of the genome.

    Parameters
    ----------
    record: SeqRecord
    start_codons: List[str]
    stop_codons: List[str]
    validate_cds: bool
        Filter out NCBI provided ORFs that do not fit our ORF criteria.
    verbose: bool

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    cds_regions = [f for f in record.features if f.type == "CDS"]

    orfs = []
    for region in cds_regions:
        loc = region.location
        seq = record.seq[loc.start:loc.end]
        if loc.strand == -1:
            seq = seq.reverse_complement()

        if not validate_cds:
            orfs.append((loc.strand, int(loc.start), int(loc.end)))
            continue

        try:
            assert seq[:3] in start_codons, "Start codon not found!"
            assert seq[-3:] in stop_codons, "Stop codon not found!"
            # Make sure there are no stop codons in the middle of the sequence
            for codon in codons(seq[3:-3]):
                assert (
                    codon not in stop_codons
                ), f"Stop codon {codon} found in the middle of the sequence!"

            # The CDS looks fine, add it to the ORFs
            orfs.append((loc.strand, int(loc.start), int(loc.end)))

        except AssertionError as ex:
            if verbose:
                print(
                    "Skipped CDS at region [%d - %d] on strand %d"
                    % (int(loc.start), int(loc.end), loc.strand)
                )
                print("\t", str(ex))

    # Some ORFs in paramecium have lenghts not divisible by 3. Remove these
    orfs = [orf for orf in orfs if (orf[2] - orf[1]) % 3 == 0]

    return orfs


def find_orfs(sequence : str, start_codons, stop_codons):
    """Find possible ORF candidates in a single reading frame.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int]]
        tuples of form (start_loc, stop_loc)

    """
    
    start_loc, stop_loc = [], []
    flag = False
    i = 0
    while i < len(sequence) - 2:
        if not flag:
            if sequence[i : i + 3] in start_codons:
                start_loc.append(i)
                flag = True
        else:
            if sequence[i : i + 3] in stop_codons:
                stop_loc.append(i + 3)
                flag = False
        i = i + 3
        
    return [(start_loc[i], stop_loc[i]) for i in range(len(stop_loc))]

def find_all_orfs(sequence : str, start_codons, stop_codons):
    """Find ALL the possible ORF candidates in the sequence using all six
    reading frames.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    orfs = []

    for i in range(3): orfs.extend([(1, j + i, k + i) for (j, k) in find_orfs(sequence[i:], start_codons, stop_codons)])

    nucleotide_map = {"A" : "T", "T" : "A", "G" : "C", "C" : "G"}
    strand_complement = "".join([nucleotide_map[x] for x in sequence])

    n = len(sequence)
    for i in range(3): orfs.extend([(-1, n - (k + i), n - (j + i)) for (j, k) in find_orfs(strand_complement[::-1][i:], start_codons, stop_codons)])
    return orfs

def translate_to_protein(seq):
    """Translate a nucleotide sequence into a protein sequence.

    Parameters
    ----------
    seq: str

    Returns
    -------
    str
        The translated protein sequence.

    """
    amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    codons =    [["GCT", "GCC", "GCA", "GCG"],
                ["TGT", "TGC"],
                ["GAT", "GAC"],
                ["GAA", "GAG"],
                ["TTT", "TTC"],
                ["GGT", "GGC", "GGA", "GGG"],
                ["CAT", "CAC"],
                ["ATT", "ATC", "ATA"],
                ["AAA", "AAG"],
                ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
                ["ATG"],
                ["AAT", "AAC"],
                ["CCT", "CCC", "CCA", "CCG"],
                ["CAA", "CAG"],
                ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
                ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
                ["ACT", "ACC", "ACA", "ACG"],
                ["GTT", "GTC", "GTA", "GTG"],
                ["TGG"],
                ["TAT", "TAC"]]
    amino_acid_sequence = ""
    for i in range(0, len(seq), 3):
        codon = seq[i : (i + 3)]
        for j in range(len(codons)):
            if codon in codons[j]: amino_acid_sequence += amino_acids[j]
    return amino_acid_sequence


def find_all_orfs_nested(sequence, start_codons, stop_codons):
    """Bonus problem: Find ALL the possible ORF candidates in the sequence using
    the updated definition of ORFs.

    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]

    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.

    """
    raise NotImplementedError()

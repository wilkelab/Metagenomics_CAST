import re
from collections import namedtuple
from typing import List, Tuple

import more_itertools
import parasail
from Bio.Seq import Seq


GAP_OPEN_PENALTY = 12
GAP_EXTEND_PENALTY = 12


cigar_regex = re.compile(r"(\d+)=")

Coordinates = List[Tuple[int, int]]

AlignmentResult = namedtuple('AlignmentResult',
                             ['score',
                              'match_count',
                              'haystack_start',
                              'haystack_end',
                              'aligned_haystack_sequence',
                              'aligned_needle_sequence',
                              'comparison_string',
                              'inverted'])


def _count_cigar_matches(string: str) -> int:
    """ parasail provides a CIGAR string to encode the alignment. We parse this to determine the number of exact matches. """
    matches = cigar_regex.findall(string)
    return sum([int(match) for match in matches])


class DNASlice(object):
    """ Represents a subsequence in a larger DNA, while retaining the global coordinates of the smaller slice. """
    def __init__(self, sequence: Seq, start: int, end: int):
        self.full_sequence = sequence
        self.start = start
        self.end = end

    def get_slice_from_local_coordinates(self, start: int, end: int):
        """ When we hand out the sequence for local alignment, we'll get coordinates relative to that slice.
        This method takes those coordinates directly and updates the global location. """
        return DNASlice(self.full_sequence, self.start + start, self.start + end)

    @property
    def sequence(self):
        return self.full_sequence[self.start: self.end]


class RepeatPair(object):
    def __init__(self, slice1: DNASlice, slice2: DNASlice, ar: AlignmentResult):
        self.slice1 = slice1
        self.slice2 = slice2
        self.ar = ar


def align(needle: str, haystack: str, inverted: bool) -> AlignmentResult:
    """ Run the local pairwise alignment of two strings and return alignment data. """
    # perform the alignment
    result = parasail.sw_trace(needle, haystack, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, parasail.nuc44)

    # extract pertinent data from the alignment result
    cigar_text = result.cigar.decode.decode("utf8")
    match_count = _count_cigar_matches(cigar_text)
    haystack_sequence = result.traceback.ref.replace("-", "")
    haystack_end = result.end_ref + 1
    haystack_start = haystack_end - len(haystack_sequence)
    return AlignmentResult(score=result.score,
                           match_count=match_count,
                           haystack_start=haystack_start,
                           haystack_end=haystack_end,
                           aligned_haystack_sequence=result.traceback.ref,
                           aligned_needle_sequence=result.traceback.query,
                           comparison_string=result.traceback.comp,
                           inverted=inverted)


def scan_for_repeats(source: DNASlice, target: DNASlice, inverted: bool, min_length: int, max_length: int, min_homology: float) -> List[RepeatPair]:
    repeats = []

    for length in range(min_length, max_length):
        for source_start_coord, nucleotides in enumerate(more_itertools.windowed(source.sequence, length)):

            if nucleotides[-1] is None:
                break
            # prepare the source sequence for alignment
            sequence = Seq("".join(nucleotides))
            sequence = str(sequence.reverse_complement()) if inverted else str(sequence)

            # align the sequences
            alignment_result = align(sequence, str(target.sequence), inverted)

            # assess the quality of this alignment

            if len(alignment_result.aligned_haystack_sequence) < length:
                continue
            homology = alignment_result.match_count / len(alignment_result.aligned_haystack_sequence)
            if homology < min_homology:
                continue

            # document the sequence and coordinates of the repeat pair
            source_end_coord = source_start_coord + length
            source_slice = source.get_slice_from_local_coordinates(source_start_coord, source_end_coord)
            target_slice = target.get_slice_from_local_coordinates(alignment_result.haystack_start, alignment_result.haystack_end)
            repeats.append(RepeatPair(source_slice, target_slice, alignment_result))

    return repeats

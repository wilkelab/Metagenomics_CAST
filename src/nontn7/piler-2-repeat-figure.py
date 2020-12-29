# Shows how well our modified piler-cr version detects arrays with two repeat-spacers

# re-run piler
# see if the first two repeat-spacers are at the same location and have the same repeat and spacer sequence
# count how many were found successfully

from collections import defaultdict
import random
import re
import sys

import jellyfish
from operon_analyzer import load, spacers

regex = re.compile(r'^Copies: (\d\d)')
random.seed(43)

REPEAT_COUNT = int(sys.argv[1])
BUFFER_SIZE = 100


def random_dna(length: int) -> str:
    return "".join([random.choice('ACGT') for _ in range(length)])


def create_array(repeat_length: int, repeat_spacer_count: int, spacer_length) -> str:
    repeat = random_dna(repeat_length)
    spacers = [random_dna(spacer_length) for _ in range(repeat_spacer_count)]
    array = "".join([f"{repeat}{spacer}" for spacer in spacers])
    print(repeat)
    print(array)
    return f"{random_dna(BUFFER_SIZE)}{array}{random_dna(BUFFER_SIZE)}"


counts = defaultdict(int)
for repeat_length in range(30, 50):
    for spacer_length in range(28, 36):
        for count in range(2):
            array_seq = create_array(repeat_length, REPEAT_COUNT, spacer_length)
            array_spacers = spacers._get_operon_spacers(0, len(array_seq), array_seq)
            print(array_spacers)
            counts[len(array_spacers)] += 1

print(counts)

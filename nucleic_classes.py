class RNA(str):
    def __init__(self, rna_seq):
        self.sequence = rna_seq.upper()
        if not all(symb in 'ACUG' for symb in self.sequence):
            raise Exception('Sequence contains inappropriate base.')

    def gc(self):
        gc = self.sequence.count('G') + self.sequence.count('C')
        try:
            return gc * 100 / len(self)
        except ZeroDivisionError:
            raise ZeroDivisionError('Unable to count GC-content for empty sequence.')

    def reverse_complement(self):
        rna_tab = str.maketrans("ACGU", "UGCA")
        return RNA(self.sequence.translate(rna_tab)[::-1])


class DNA(RNA):
    def __init__(self, dna_seq):
        self.sequence = dna_seq.upper()
        if not all(symb in 'ACTG' for symb in self.sequence):
            raise ZeroDivisionError('Sequence contains inappropriate base.')

    def transcribe(self):
        return RNA(self.sequence.replace('T', 'U'))

    def reverse_complement(self):
        dna_tab = str.maketrans("ACTG", "TGAC")
        return DNA(self.sequence.translate(dna_tab)[::-1])

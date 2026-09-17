"""THE CONTRACT: a reported beta is per copy of the declared `effect_allele`.

Every reader must count the allele it names. That is the whole agreement, and
it is what makes results comparable across input formats without anyone having
to know which column of which file the dosage came from.

It is deliberately NOT "every format counts ALT" or "every format counts A2".
The formats name their alleles differently -- PLINK1 `.bim` has A1/A2 columns,
PGEN has REF/ALT, BGEN carries a list -- and forcing one spelling on all of
them would mean flipping the sign of results a user already has. What has to
hold is the weaker, checkable property: dosage and declaration agree.

WHY THIS FILE EXISTS. On 2026-09-15 the BED path was found to produce dosages
summing to exactly 2.0 against the PGEN path on the same cohort -- opposite
alleles -- and the first fix attempted was to flip BED's decode table. That was
wrong: BED counts A2 and *declares* A2, verified against the official PLINK
example, so it already honoured the contract. Flipping the table alone left it
counting A1 while announcing A2, which is the one genuinely broken state. A
test of the invariant would have caught that in seconds; reading the decode
table did not.

The consequence a user must understand is downstream of the contract, not a
violation of it: `.bed` and `.pgen` built from the same data name DIFFERENT
effect alleles, so betas differ in sign between them. The `effect_allele`
column in the output is what resolves that, and it is always written.
"""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, "src")

import numpy as np

from torchgwas.bed import PlinkBedGenotype, resolve_plink_triplet


def write_bed_with_alleles(prefix: Path, calls, a1, a2):
    """A one-variant BED whose genotypes are given as A2 counts.

    `calls` is the A2 dosage per sample; the PLINK codes are
    00=A1/A1 (A2 dosage 0), 10=het (1), 11=A2/A2 (2), 01=missing.
    """
    bed, bim, fam = resolve_plink_triplet(prefix)
    samples = len(calls)
    fam.write_text("".join(f"F{i} I{i} 0 0 0 -9\n" for i in range(samples)))
    bim.write_text(f"1 rs1 0 1 {a1} {a2}\n")
    code_for = {0: 0b00, 1: 0b10, 2: 0b11}
    payload = bytearray(b"\x6c\x1b\x01")
    for start in range(0, samples, 4):
        byte = 0
        for offset in range(4):
            sample = start + offset
            code = 0b01 if sample >= samples else code_for[calls[sample]]
            byte |= code << (2 * offset)
        payload.append(byte)
    bed.write_bytes(bytes(payload))
    return bed


class EffectAlleleContractTestCase(unittest.TestCase):
    def test_bed_counts_the_allele_it_declares(self):
        # A2 = "T". Samples are 0, 1, 2 copies of T.
        with tempfile.TemporaryDirectory() as tmp:
            bed = write_bed_with_alleles(Path(tmp) / "c", [0, 1, 2, 0],
                                         a1="G", a2="T")
            source = PlinkBedGenotype(bed)
            declared = source.effect_alleles[0]
            dosage = np.asarray(source.read_chunk(0, 1)).reshape(-1)[:4]
            self.assertEqual(declared, "T",
                             "BED must declare the .bim column-6 allele")
            np.testing.assert_allclose(dosage, [0, 1, 2, 0],
                                       err_msg="dosage must COUNT the declared "
                                               "effect allele, not the other one")

    def test_bed_declaration_follows_the_bim_not_a_constant(self):
        # Swap the columns: the declared allele must swap with them, and the
        # dosage must follow, or the reader is counting a fixed position rather
        # than honouring the file.
        with tempfile.TemporaryDirectory() as tmp:
            bed = write_bed_with_alleles(Path(tmp) / "d", [0, 1, 2, 0],
                                         a1="T", a2="G")
            source = PlinkBedGenotype(bed)
            self.assertEqual(source.effect_alleles[0], "G")
            self.assertEqual(source.other_alleles[0], "T")

    def test_official_plink_example_counts_its_declared_allele(self):
        # cog-genomics.org/plink/1.9/formats#bed. bim "G A" so A2 = A; the PED
        # is GG, AA, missing, AA, AA, AA, which is 0, 2, missing, 2, 2, 2
        # copies of A. This is the authoritative check that the declaration and
        # the count agree, on bytes nobody here wrote.
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "official"
            Path(f"{prefix}.bed").write_bytes(
                bytes.fromhex("6c1b01dc0fe70f6b01"))
            Path(f"{prefix}.bim").write_text(
                "1 snp1 0 1 G A\n1 snp2 0 2 1 2\n1 snp3 0 3 A C\n")
            Path(f"{prefix}.fam").write_text(
                "".join(f"0 I{i} 0 0 0 -9\n" for i in range(6)))
            source = PlinkBedGenotype(prefix)
            self.assertEqual(source.effect_alleles[0], "A")
            np.testing.assert_allclose(
                np.asarray(source.read_chunk(0, 1))[:, 0],
                [0, 2, np.nan, 2, 2, 2], equal_nan=True)

    def test_every_reader_exposes_the_declaration_at_all(self):
        # A reader that does not name its effect allele cannot honour the
        # contract, because nothing downstream can state what the beta is per
        # copy of. This is the cheap structural half of the check.
        from torchgwas import bed as bed_module
        from torchgwas import pgen as pgen_module

        for module in (bed_module, pgen_module):
            source = module.__name__.rsplit(".", 1)[-1]
            names = [n for n in dir(module) if "effect_allele" in n.lower()]
            has_attr = any("effect_allele" in line
                           for line in open(module.__file__, encoding="utf-8"))
            self.assertTrue(has_attr or names,
                            f"{source} never mentions effect_alleles")


if __name__ == "__main__":
    unittest.main()

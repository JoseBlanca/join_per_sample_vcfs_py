import tempfile
from pathlib import Path

from join_vcfs.vcf_joining import _create_vcf_infos, _generate_var_bins

VCF1 = b"""#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001
20\t1\t.\tG\tA\t20\tPASS\t.\tGT\t0|0
20\t2\t.\tT\tA\t20\tPASS\t.\tGT\t.|0
20\t3\t.\tA\tG\t20\tPASS\t.\tGT\t.|0
20\t4\t.\tT\t.\t20\tPASS\t.\tGT\t0|0
20\t5\t.\tG\tA\t20\tPASS\t.\tGT\t0/1
20\t6\t.\tG\tA\t20\tPASS\t.\tGT\t0/1"""

VCF2 = b"""#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00002
20\t1\t.\tG\tA\t20\tPASS\t.\tGT\t0|0
20\t2\t.\tT\tA\t20\tPASS\t.\tGT\t.|0
20\t3\t.\tA\tG\t20\tPASS\t.\tGT\t.|0
20\t4\t.\tT\t.\t20\tPASS\t.\tGT\t0|0
20\t5\t.\tG\tA\t20\tPASS\t.\tGT\t0/1
20\t6\t.\tG\tA\t20\tPASS\t.\tGT\t0/1"""

VCF3 = b"""#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00003
20\t1\t.\tGATCGAT\tA\t20\tPASS\t.\tGT\t0|0"""


def get_var_bin_spans(vcf_paths, sorted_chromosomes):
    remaining_chromosomes = sorted_chromosomes[:]
    vcf_infos = _create_vcf_infos(vcf_paths)
    bin_spans = [
        bin.span for bin in _generate_var_bins(vcf_infos, remaining_chromosomes)
    ]

    return bin_spans


def test_simple_binning():
    with (
        tempfile.NamedTemporaryFile() as tmp1,
        tempfile.NamedTemporaryFile() as tmp2,
    ):
        tmp1.write(VCF1)
        tmp1.flush()
        tmp1_path = Path(tmp1.name)
        tmp2.write(VCF2)
        tmp2.flush()
        tmp2_path = Path(tmp2.name)

        bin_spans = get_var_bin_spans(
            [tmp1_path, tmp2_path], sorted_chromosomes=["1", "20"]
        )
        print(bin_spans)
        assert bin_spans == [
            ("20", 1, 1),
            ("20", 2, 2),
            ("20", 3, 3),
            ("20", 4, 4),
            ("20", 5, 5),
            ("20", 6, 6),
        ]


def test_binning_with_deletion_spanning_several_snps():
    with (
        tempfile.NamedTemporaryFile() as tmp1,
        tempfile.NamedTemporaryFile() as tmp2,
        tempfile.NamedTemporaryFile() as tmp3,
    ):
        tmp1.write(VCF1)
        tmp1.flush()
        tmp1_path = Path(tmp1.name)
        tmp3.write(VCF3)
        tmp3.flush()
        tmp3_path = Path(tmp3.name)

        bin_spans = get_var_bin_spans(
            [tmp1_path, tmp3_path], sorted_chromosomes=["1", "20"]
        )
        print(bin_spans)
        assert bin_spans == [
            ("20", 1, 7),
        ]

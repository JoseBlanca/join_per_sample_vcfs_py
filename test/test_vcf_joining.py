import tempfile
from pathlib import Path

import pytest

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

VCF4 = b"""#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00004
1\t4\t.\tGATCGAT\tA\t20\tPASS\t.\tGT\t0|0
20\t4\t.\tGATCGAT\tA\t20\tPASS\t.\tGT\t0|0
21\t4\t.\tGATCGAT\tA\t20\tPASS\t.\tGT\t0|0
"""

VCF5 = b"""#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00005
20\t3\t.\tG\tA\t20\tPASS\t.\tGT\t0|0
20\t8\t.\tG\tA\t20\tPASS\t.\tGT\t0|0
20\t20\t.\tG\tA\t20\tPASS\t.\tGT\t0|0"""


def get_var_bins(vcf_paths, sorted_chromosomes):
    remaining_chromosomes = sorted_chromosomes[:]
    vcf_infos = _create_vcf_infos(vcf_paths)
    var_bins = list(_generate_var_bins(vcf_infos, remaining_chromosomes))
    bin_spans = [bin.span for bin in var_bins]
    bin_vars = [bin.vars for bin in var_bins]

    return bin_spans, bin_vars


def write_in_temp_file(temp_file, contents):
    temp_file.write(contents)
    temp_file.flush()
    tmp_path = Path(temp_file.name)
    return tmp_path


def test_simple_binning():
    with (
        tempfile.NamedTemporaryFile() as tmp1,
        tempfile.NamedTemporaryFile() as tmp2,
    ):
        tmp1_path = write_in_temp_file(tmp1, VCF1)
        tmp2_path = write_in_temp_file(tmp2, VCF2)

        bin_spans = get_var_bins(
            [tmp1_path, tmp2_path], sorted_chromosomes=["1", "20"]
        )[0]
        assert bin_spans == [
            ("20", 1, 1),
            ("20", 2, 2),
            ("20", 3, 3),
            ("20", 4, 4),
            ("20", 5, 5),
            ("20", 6, 6),
        ]

    with (
        tempfile.NamedTemporaryFile() as tmp1,
    ):
        tmp1_path = write_in_temp_file(tmp1, VCF1)

        bin_spans = get_var_bins([tmp1_path], sorted_chromosomes=["20"])[0]
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
        tempfile.NamedTemporaryFile() as tmp3,
    ):
        tmp1_path = write_in_temp_file(tmp1, VCF1)
        tmp3_path = write_in_temp_file(tmp3, VCF3)

        bin_spans = get_var_bins(
            [tmp1_path, tmp3_path], sorted_chromosomes=["1", "20"]
        )[0]
        assert bin_spans == [
            ("20", 1, 7),
        ]

    with (
        tempfile.NamedTemporaryFile() as tmp1,
        tempfile.NamedTemporaryFile() as tmp3,
        tempfile.NamedTemporaryFile() as tmp4,
        tempfile.NamedTemporaryFile() as tmp5,
    ):
        tmp1_path = write_in_temp_file(tmp1, VCF1)
        tmp3_path = write_in_temp_file(tmp3, VCF3)
        tmp4_path = write_in_temp_file(tmp4, VCF4)
        tmp5_path = write_in_temp_file(tmp5, VCF5)

        bin_spans = get_var_bins(
            [tmp1_path, tmp3_path, tmp4_path, tmp5_path],
            sorted_chromosomes=["1", "20", "21"],
        )[0]
        assert bin_spans == [
            ("1", 4, 10),
            ("20", 1, 10),
            ("20", 20, 20),
            ("21", 4, 10),
        ]


VCF_WITH_WRONG_CHROMOSOME_ORDER = b"""#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00005
1\t3\t.\tG\tA\t20\tPASS\t.\tGT\t0|0
2\t8\t.\tG\tA\t20\tPASS\t.\tGT\t0|0
1\t20\t.\tG\tA\t20\tPASS\t.\tGT\t0|0"""


def test_wrong_chrom_order():
    with (
        tempfile.NamedTemporaryFile() as tmp1,
    ):
        tmp1_path = write_in_temp_file(tmp1, VCF_WITH_WRONG_CHROMOSOME_ORDER)

        with pytest.raises(RuntimeError):
            get_var_bins([tmp1_path], sorted_chromosomes=["1", "2"])


VCF_WITH_WRONG_ORDER = b"""#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00005
1\t3\t.\tG\tA\t20\tPASS\t.\tGT\t0|0
1\t4\t.\tG\tA\t20\tPASS\t.\tGT\t0|0
1\t5\t.\tG\tA\t20\tPASS\t.\tGT\t0|0
1\t8\t.\tG\tA\t20\tPASS\t.\tGT\t0|0
1\t2\t.\tG\tA\t20\tPASS\t.\tGT\t0|0"""


def test_wrong_order():
    with (
        tempfile.NamedTemporaryFile() as tmp1,
    ):
        tmp1_path = write_in_temp_file(tmp1, VCF_WITH_WRONG_ORDER)

        with pytest.raises(RuntimeError):
            get_var_bins([tmp1_path], sorted_chromosomes=["1"])


# TODO
#
# ------
#  -  -
#
# -------
#  -   -----
#  -- ---  -
#
# -------   different chrom
#   -       -
#
# -------
#    -      no more vars
#
#

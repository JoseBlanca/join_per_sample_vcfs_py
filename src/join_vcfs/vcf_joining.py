from pathlib import Path
from typing import Iterator
from collections import defaultdict, namedtuple

from more_itertools import peekable

from join_vcfs.vcf_parser import parse_vcf


def _fill_incoming_snps(incomping_snps, snp_iterators, iterator_ids):
    for id_ in iterator_ids:
        if id_ not in incomping_snps:
            try:
                incomping_snps[id_] = next(snp_iterators[id_]["vars_iter"])
            except StopIteration:
                pass


def _overlaps(var_span, bin_span):
    if var_span[0] != bin_span[0]:
        return False
    var_start, var_end = var_span[1], var_span[2]
    bin_start, bin_end = bin_span[1], bin_span[2]
    if var_start < bin_start:
        msg = "Implementation error, a var has a start lower than the current var bin"
        msg += f"var: {var_span[0]}:{var_start}-{var_end} bin_span: {bin_span[0]}:{bin_start}-{bin_end}"
        raise RuntimeError(msg)
    if var_start <= bin_end:
        return True
    else:
        return False


def _calculate_var_span(var):
    alleles_len = max(len(allele) for allele in var["alleles"])
    return var["chrom"], var["pos"], var["pos"] + alleles_len - 1


def _add_var_to_bin(vars_in_bin, var, var_span, iterator_id, bin_span):
    vars_in_bin[iterator_id].append(var)
    var_span = _calculate_var_span(var)
    span_has_been_elongated = False
    if var_span[1] > bin_span[1]:
        bin_span[1] = var_span[1]
        span_has_been_elongated = True
    return span_has_been_elongated


NO_VARS_LEFT = 1
NO_VARS_IN_CURRENT_CHROM = 2


def _get_first_span(vcf_infos, current_chrom):
    first_span = None
    no_vars_left = True
    no_vars_in_current_chrom = True

    for vcf_info in vcf_infos.values():
        try:
            next_var = vcf_info["vars_iter"].peek()
            no_vars_left = False
        except StopIteration:
            continue
        if next_var["chrom"] != current_chrom:
            continue
        no_vars_in_current_chrom = False
        next_var_span = _calculate_var_span(next_var)
        if first_span is None:
            first_span = next_var_span
        else:
            if next_var_span[1] < first_span[1]:
                first_span = next_var_span
    return first_span, no_vars_left, no_vars_in_current_chrom


VarBin = namedtuple("VarBin", ["vars", "span"])


def _create_vars_bin(vcf_infos, current_chrom):
    first_span, no_vars_left, no_vars_in_current_chrom = _get_first_span(
        vcf_infos, current_chrom
    )

    if no_vars_in_current_chrom:
        return NO_VARS_IN_CURRENT_CHROM
    elif no_vars_left:
        return NO_VARS_LEFT

    vars_in_bin = defaultdict(list)
    span_has_been_elongated = False
    bin_span = first_span
    while True:
        for vcf_id, vcf_info in vcf_infos.items():
            try:
                next_var = vcf_info["vars_iter"].peek()
            except StopIteration:
                continue
            next_var_span = _calculate_var_span(next_var)
            if _overlaps(next_var_span, bin_span):
                try:
                    next_var = next(vcf_info["vars_iter"])
                except StopIteration:
                    msg = "Implementation error, we have previously peeked the var iterator and we made sure that a var was coming"
                    raise RuntimeError(msg)
                span_has_been_elongated_this_time = _add_var_to_bin(
                    vars_in_bin, next_var, next_var_span, vcf_id, bin_span
                )
                if span_has_been_elongated_this_time:
                    span_has_been_elongated = True
        if not span_has_been_elongated:
            break
    return VarBin(vars_in_bin, bin_span)


def _generate_var_bins(
    vcf_infos: dict[int, dict],
    remaining_chromosomes: list[str],
):
    remaining_chromosomes = remaining_chromosomes[:]

    current_chrom = remaining_chromosomes.pop(0)
    while current_chrom is not None:
        res = _create_vars_bin(vcf_infos, current_chrom)
        if res == NO_VARS_IN_CURRENT_CHROM:
            if remaining_chromosomes:
                current_chrom = remaining_chromosomes.pop(0)
            else:
                current_chrom = None
            continue
        elif res == NO_VARS_LEFT:
            break
        else:
            # This is a vars_bin
            yield res


def _create_vcf_infos(vcf_paths):
    vcf_paths = [Path(path) for path in vcf_paths]
    parsing_results = [parse_vcf(path) for path in vcf_paths]

    vcf_infos = {}
    samples_seen = set()
    for idx, result in enumerate(parsing_results):
        metadata = result["metadata"]
        this_samples = list(map(str, metadata["samples"]))

        overlaping_samples = samples_seen.intersection(this_samples)
        if overlaping_samples:
            raise RuntimeError(
                "Some samples are found in different VCFs: ",
                ",".join(overlaping_samples),
            )
        samples_seen.update(this_samples)

        vcf_info = {
            "vars_iter": peekable(result["vars"]),
            "samples": metadata["samples"],
            "ploidy": metadata["ploidy"],
            "fhand": result["fhand"],
        }
        vcf_infos[idx] = vcf_info
    return vcf_infos


def join_vcfs(vcf_paths: list[Path], ordered_chromosomes: list) -> Iterator:
    if not ordered_chromosomes:
        raise ValueError("Al least one chromosome should be given")
    remaining_chromosomes = ordered_chromosomes[:]
    vcf_infos = _create_vcf_infos(vcf_paths)

    var_bins = _generate_var_bins(vcf_infos, remaining_chromosomes)

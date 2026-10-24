"""Read-only, standard-library validation of an explicitly supplied table bundle."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path


TABLES = ("toxicant_cpg_evidence_curated", "controlled_id_vocabulary_crosswalk",
          "exposure_definition_qualifiers", "candidate_analysis_guards")
MISSING = {"not_reported", "not_applicable"}
DICTIONARY_COLUMNS = ["table_name", "column_name", "data_type", "definition", "allowed_values", "missing_code", "requirement"]
PLATFORM_GROUPS = {
    "Illumina Infinium HumanMethylation450 BeadChip": "450K",
    "Illumina Infinium HumanMethylation450": "450K",
    "Illumina HumanMethylation450": "450K",
    "Illumina Infinium MethylationEPIC BeadChip": "EPIC",
    "Illumina Infinium MethylationEPIC": "EPIC",
    "Illumina HumanMethylationEPIC": "EPIC",
}


def digest_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames or []
        if not columns or len(set(columns)) != len(columns):
            raise ValueError("empty or duplicated CSV header")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError("CSV row width differs from its header")
    return columns, rows


def numeric(raw: str) -> Decimal | None:
    if raw in MISSING:
        return None
    try:
        value = Decimal(raw)
    except InvalidOperation as error:
        raise ValueError("invalid numeric representation") from error
    if not value.is_finite():
        raise ValueError("non-finite numeric value")
    return value


def nonnegative_integer(raw: str) -> int:
    value = numeric(raw)
    if value is None or value < 0 or value != value.to_integral_value():
        raise ValueError("invalid nonnegative integer")
    return int(value)


def indexed(rows: list[dict[str, str]], keys: list[str]) -> dict[tuple[str, ...], dict[str, str]]:
    result = {}
    for row in rows:
        key = tuple(row[column] for column in keys)
        if not all(key) or key in result:
            raise ValueError("empty or duplicated primary key")
        result[key] = row
    return result


def check_evidence(rows: list[dict[str, str]]) -> dict:
    zero_count, non_cg_count, positive_se = 0, 0, 0
    for row in rows:
        if any(not value for value in row.values()):
            raise ValueError("empty required evidence value")
        effect, standard_error, p_value, sample_size = (
            numeric(row[field]) for field in ("effect_estimate", "standard_error", "p_value", "sample_size"))
        if p_value is not None and not Decimal(0) <= p_value <= Decimal(1):
            raise ValueError("p-value outside zero-to-one range")
        if standard_error is not None and standard_error < 0:
            raise ValueError("negative standard error")
        if sample_size is not None and (sample_size <= 0 or sample_size != sample_size.to_integral_value()):
            raise ValueError("sample size is not a positive integer")
        if effect is not None and effect != 0 and row["effect_direction"] != ("positive" if effect > 0 else "negative"):
            raise ValueError("effect direction differs from estimate sign")
        zero_count += p_value == 0
        non_cg_count += not bool(re.fullmatch(r"cg[0-9]{8}", row["cpg_id"]))
        positive_se += standard_error is not None and standard_error > 0
    return {"source_zero_pvalues": zero_count, "non_cg_identifiers": non_cg_count,
            "positive_reported_standard_errors": positive_se}


def check_guards(evidence: dict, rows: list[dict[str, str]]) -> None:
    guarded = indexed(rows, ["evidence_id"])
    if set(guarded) != set(evidence):
        raise ValueError("guard coverage differs from evidence")
    for key, guard in guarded.items():
        source = evidence[key]
        for field, source_field in (("pmid", "pmid"), ("cpg_id", "cpg_id"),
                                    ("declared_platform", "methylation_platform"),
                                    ("source_p_value", "p_value"), ("source_standard_error", "standard_error")):
            if guard[field] != source[source_field]:
                raise ValueError("guard source value differs from evidence")
        cg_form = bool(re.fullmatch(r"cg[0-9]{8}", source["cpg_id"]))
        mapped = guard["manifest_mapping_status"] == "mapped"
        if mapped and not cg_form:
            raise ValueError("non-cg identifier marked as CpG-only mapped")
        if guard["cpg_only_lookup_candidate"] != str(cg_form and mapped).lower():
            raise ValueError("CpG-only guard differs from reported mapping")
        p_value, standard_error = numeric(source["p_value"]), numeric(source["standard_error"])
        representation = ("controlled_missing_value" if p_value is None else
                          "numeric_zero_exact_precision_unknown" if p_value == 0 else "positive_numeric_value_as_reported")
        if (guard["p_value_representation"] != representation
                or guard["p_value_transform_candidate"] != str(p_value is not None and p_value > 0).lower()
                or guard["inverse_variance_weight_candidate"] != str(standard_error is not None and standard_error > 0).lower()):
            raise ValueError("numeric operation guard differs from source representation")
        if guard["platform_group"] not in {"450K", "EPIC", "unknown"}:
            raise ValueError("unknown platform-group enumeration")
        if guard["platform_group"] != PLATFORM_GROUPS.get(source["methylation_platform"], "unknown"):
            raise ValueError("platform group differs from declared platform")
        if guard["platform_group"] == "unknown" and mapped:
            raise ValueError("unknown platform marked as mapped")


def check_companions(evidence: dict, crosswalk: list[dict[str, str]], qualifiers: list[dict[str, str]], contract: dict) -> None:
    counts = Counter(row["toxicant_id"] for row in evidence.values())
    exposures = indexed(crosswalk, ["release_toxicant_id"])
    if {key[0] for key in exposures} != set(counts):
        raise ValueError("crosswalk exposure coverage differs")
    first_total, second_total = 0, 0
    for key, row in exposures.items():
        first, second = nonnegative_integer(row["first_round_rows"]), nonnegative_integer(row["second_round_rows"])
        if first + second != counts[key[0]] or row["external_identifier_transfer_permitted"] != "false":
            raise ValueError("crosswalk counts or identifier-transfer boundary differs")
        first_total += first
        second_total += second
    rounds = contract["tables"][TABLES[0]].get("provenance_by_round", {})
    if rounds and (first_total != rounds["first_round"]["records"] or second_total != rounds["second_round"]["records"]):
        raise ValueError("crosswalk round totals differ from contract")
    required = {key for key, row in evidence.items() if row["toxicant_id"] in {"urinary_total_arsenic", "air_pm_coarse"}}
    qualified = indexed(qualifiers, ["evidence_id"])
    if set(qualified) != required:
        raise ValueError("qualifier coverage differs from limited exposure scope")
    allowed = contract["tables"][TABLES[2]]["analyte_definition_values"]
    for key, row in qualified.items():
        source = evidence[key]
        if (row["release_toxicant_id"] != source["toxicant_id"] or row["pmid"] != source["pmid"]
                or row["source_locator"] != source["source_locator"] or row["analyte_definition"] not in allowed):
            raise ValueError("qualifier identity or definition differs")
        if source["toxicant_id"] == "air_pm_coarse" and row["analyte_definition"] != "not_applicable":
            raise ValueError("coarse-PM row has an arsenic analyte definition")
        if source["toxicant_id"] == "urinary_total_arsenic" and row["analyte_definition"] == "not_applicable":
            raise ValueError("urinary arsenic row lacks an analyte definition")


def validate(paths: dict[str, Path]) -> dict:
    if set(paths) != {*TABLES, "contract", "dictionary"}:
        raise ValueError("exactly four tables, one contract and one dictionary are required")
    before = {role: digest_file(path) for role, path in paths.items()}
    contract = json.loads(paths["contract"].read_text(encoding="utf-8"))
    if set(contract["tables"]) != set(TABLES):
        raise ValueError("unsupported contract table set")
    datasets = {table: read_csv(paths[table]) for table in TABLES}
    for table, (columns, rows) in datasets.items():
        specification = contract["tables"][table]
        if columns != specification["required_columns"] or len(rows) != specification["rows"]:
            raise ValueError("table header or row count differs from contract")
        keys = specification["primary_key"]
        indexed(rows, [keys] if isinstance(keys, str) else keys)
    evidence = indexed(datasets[TABLES[0]][1], ["evidence_id"])
    checks = check_evidence(datasets[TABLES[0]][1])
    check_guards(evidence, datasets[TABLES[3]][1])
    check_companions(evidence, datasets[TABLES[1]][1], datasets[TABLES[2]][1], contract)
    columns, definitions = read_csv(paths["dictionary"])
    if columns != DICTIONARY_COLUMNS:
        raise ValueError("dictionary header differs")
    expected = {(table, column) for table in TABLES for column in contract["tables"][table]["required_columns"]}
    actual = indexed(definitions, ["table_name", "column_name"])
    if set(actual) != expected or any(not row["definition"] or not row["data_type"] or row["requirement"] != "required" for row in definitions):
        raise ValueError("dictionary does not define exactly every required table field")
    if any(digest_file(path) != before[role] for role, path in paths.items()):
        raise ValueError("input changed during validation")
    return {
        "status": "table_bundle_structure_and_consistency_validated_not_release_clearance",
        "table_rows": {table: len(rows) for table, (_, rows) in datasets.items()},
        "dictionary_fields": len(definitions), "numeric_checks": checks,
        "source_values_modified": False, "source_extraction_rebuilt": False,
        "original_review_authority_revalidated": False, "manifest_lookup_reperformed": False,
        "per_record_round_assignment_revalidated": False,
        "source_rights_assessed": False, "publication_permitted": False,
        "inputs": {role: {"filename": path.name, "sha256": before[role]} for role, path in paths.items()},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    for name in ("evidence", "crosswalk", "qualifiers", "guards", "contract", "dictionary"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    paths = dict(zip(TABLES, (args.evidence, args.crosswalk, args.qualifiers, args.guards)))
    paths.update(contract=args.contract, dictionary=args.dictionary)
    print(json.dumps(validate(paths), indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
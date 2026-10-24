import copy
import csv
import json
import sys
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts import validate_candidate_tables as validator


class CandidateTableValidatorTests(unittest.TestCase):
    def fixture(self):
        evidence_columns = ["evidence_id", "toxicant_id", "pmid", "cpg_id", "source_locator", "effect_estimate",
                            "effect_direction", "standard_error", "p_value", "sample_size", "methylation_platform"]
        first = dict(zip(evidence_columns, ["row1", "urinary_total_arsenic", "1", "cg00000001", "table::row1",
                        "0.1", "positive", "0.2", "1e-400", "20", "Illumina Infinium HumanMethylation450 BeadChip"]))
        second = {**first, "evidence_id": "row2", "cpg_id": "ch.synthetic", "source_locator": "table::row2", "p_value": "0"}
        crosswalk = [{"release_toxicant_id": "urinary_total_arsenic", "first_round_rows": "1", "second_round_rows": "1",
                      "external_identifier_transfer_permitted": "false"}]
        qualifiers = [{"evidence_id": row["evidence_id"], "release_toxicant_id": row["toxicant_id"], "pmid": row["pmid"],
                       "source_locator": row["source_locator"], "analyte_definition": "total_arsenic_gfaas"} for row in (first, second)]
        guards = []
        for row in (first, second):
            cg_form = row is first
            guards.append({"evidence_id": row["evidence_id"], "pmid": row["pmid"], "cpg_id": row["cpg_id"],
                "declared_platform": row["methylation_platform"], "source_p_value": row["p_value"],
                "source_standard_error": row["standard_error"], "platform_group": "450K",
                "manifest_mapping_status": "mapped" if cg_form else "unmapped",
                "cpg_only_lookup_candidate": str(cg_form).lower(), "p_value_transform_candidate": str(cg_form).lower(),
                "p_value_representation": "positive_numeric_value_as_reported" if cg_form else "numeric_zero_exact_precision_unknown",
                "inverse_variance_weight_candidate": "true"})
        datasets = dict(zip(validator.TABLES, ([first, second], crosswalk, qualifiers, guards)))
        contract = {"tables": {table: {"required_columns": list(rows[0]), "rows": len(rows),
                     "primary_key": ["release_toxicant_id" if table == validator.TABLES[1] else "evidence_id"]}
                    for table, rows in datasets.items()}}
        contract["tables"][validator.TABLES[0]]["provenance_by_round"] = {"first_round": {"records": 1}, "second_round": {"records": 1}}
        contract["tables"][validator.TABLES[2]]["analyte_definition_values"] = ["total_arsenic_gfaas", "sum_iAs_MMA_DMA", "not_applicable"]
        dictionary = [{"table_name": table, "column_name": column, "data_type": "string", "definition": "Synthetic field",
                       "allowed_values": "", "missing_code": "", "requirement": "required"}
                      for table, rows in datasets.items() for column in rows[0]]
        return datasets, contract, dictionary

    def write_fixture(self, directory, datasets, contract, dictionary):
        paths = {}
        for table, rows in {**datasets, "dictionary": dictionary}.items():
            path = Path(directory) / (table + ".csv")
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            paths[table] = path
        paths["contract"] = Path(directory) / "contract.json"
        paths["contract"].write_text(json.dumps(contract), encoding="utf-8")
        return paths

    def test_synthetic_bundle_preserves_tiny_positive_and_zero_boundaries(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = self.write_fixture(directory, *self.fixture())
            report = validator.validate(paths)
        self.assertEqual(validator.numeric("1e-400"), Decimal("1e-400"))
        self.assertEqual(report["numeric_checks"]["source_zero_pvalues"], 1)
        self.assertEqual(report["numeric_checks"]["non_cg_identifiers"], 1)
        self.assertFalse(report["publication_permitted"])
        self.assertFalse(report["original_review_authority_revalidated"])

    def test_rejects_zero_transform_and_source_value_changes(self):
        for field, value in (("p_value_transform_candidate", "true"), ("source_p_value", "1e-400"),
                             ("cpg_only_lookup_candidate", "true")):
            datasets, contract, dictionary = self.fixture()
            datasets[validator.TABLES[3]][1][field] = value
            with tempfile.TemporaryDirectory() as directory:
                paths = self.write_fixture(directory, datasets, contract, dictionary)
                with self.assertRaises(ValueError):
                    validator.validate(paths)

    def test_rejects_duplicate_key_negative_error_and_bad_dictionary(self):
        for case in ("duplicate", "negative_error", "missing_dictionary"):
            datasets, contract, dictionary = self.fixture()
            if case == "duplicate":
                datasets[validator.TABLES[0]][1] = copy.deepcopy(datasets[validator.TABLES[0]][0])
            elif case == "negative_error":
                datasets[validator.TABLES[0]][0]["standard_error"] = "-0.1"
            else:
                dictionary.pop()
            with tempfile.TemporaryDirectory() as directory:
                paths = self.write_fixture(directory, datasets, contract, dictionary)
                with self.assertRaises(ValueError):
                    validator.validate(paths)

    def test_csv_reader_rejects_extra_cells_and_nonfinite_numbers(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.csv"
            path.write_text("field\nvalue,extra\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "width"):
                validator.read_csv(path)
        for value in ("NaN", "Infinity", "-Infinity"):
            with self.assertRaises(ValueError):
                validator.numeric(value)

    def test_unknown_declared_platform_cannot_default_to_450k(self):
        datasets, contract, dictionary = self.fixture()
        datasets[validator.TABLES[0]][0]["methylation_platform"] = "unrecognized_platform"
        datasets[validator.TABLES[3]][0]["declared_platform"] = "unrecognized_platform"
        with tempfile.TemporaryDirectory() as directory:
            paths = self.write_fixture(directory, datasets, contract, dictionary)
            with self.assertRaisesRegex(ValueError, "platform group differs"):
                validator.validate(paths)


if __name__ == "__main__":
    unittest.main()
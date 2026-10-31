# Provenance-aware toxicant-methylation evidence resource

Version 1.0.1 is based on the frozen 2026-10-01 scientific snapshot. Public repository: https://github.com/shizhf999/toxicant-methylation-evidence. This independent data and validation resource is not tied to a particular journal. GitHub releases do not themselves mint a DOI; cite an archive DOI only after its record has been published and verified.

## Contents and counting units

The association table contains 16,681 rows and 26 fields from 12 publications, tracked through 13 source-artifact entries. It covers 22 exposure concepts and 16,053 distinct probe identifiers. A row is a published association, not an independent participant, publication or biological discovery. The crosswalk has 22 rows, definition qualifiers 109 rows, and operation guards 16,681 rows. The reader dictionary defines the 58 fields across these four tables. The unchanged dictionary is also retained as data_dictionary.csv for exact-file license identification.

The separate reference/clock_probe_membership.csv contains 937 clock-probe memberships representing 889 distinct identifiers; it supplies no coefficients or scoring implementation. Do not join clock membership as independent association evidence.

## Validation

Use Python 3.10 or newer; the validator and synthetic tests use only the standard library. From this folder run:

```console
python scripts/validate_candidate_tables.py --evidence data/association_evidence.csv --crosswalk data/exposure_crosswalk.csv --qualifiers data/exposure_qualifiers.csv --guards data/operation_guards.csv --contract config/table_contract.json --dictionary data_dictionary_reader.csv
python -m unittest discover -s tests -p test_validate_candidate_tables.py
```

Checks cover dimensions, identifiers, joins, exact numeric representations and operation prerequisites. A passing check is not independent re-extraction, biological validation or a sharing authorization. file_manifest.csv records hashes of the actual bundled files, not the original input hashes.

## Interpretation and metadata preparation

Numerical strings, row order and all 25 non-notes fields are unchanged. Source-review notes retain scientific qualifications while replacing an internal mapping filename/checksum and coordinator-letter references with self-contained descriptions. Retrieval wording is clarified without asserting a publisher reuse grant. Reviewer identifiers are workflow slots, not personal identities or proof that every final field was independently supplied twice. Earlier source-review statements describe their review stage; the final combined study identifiers are coordinator-issued. Scientific coordinator notes are retained.

Join guards one-to-one on evidence_id; join the crosswalk many-to-one from toxicant_id to release_toxicant_id. Qualifiers cover only 109 scoped rows and are joined by evidence_id, with PMID and exposure key checked; absence does not establish an unqualified measurement. Do not transfer CAS or ChEBI identifiers merely from these vocabulary links.

Retain 326 source-stored zero p-values of unknown precision, without guessed bounds or replacements. Do not impute standard errors; 5,248 rows contain positive reported standard errors. Keep ch.6.19775382R outside cg-only lookup. The 29 suspended records remain excluded. A selected published list is not a full tested-probe universe; absent probes are not null findings. No enrichment, pooling, clock scoring or participant-level fitting is performed or authorized by these prerequisites.

source_manifest.csv provides publication identifiers, source URLs, artifact names and registered hashes without private local paths. Acquisition dates and versions not established in the source index remain not_reported. PMID30055357 retains a historical source-file registration gap for 13 rows; a locally checked candidate PDF is not silently substituted for the registered source.

## Scope, licenses and limitations

No original publication PDFs, XML, supplementary workbooks, publisher figures, manuscript files, manuscript-specific analysis outputs, correspondence, private handoffs, signed records, credentials, remote-transfer helpers, installed libraries or binaries are bundled.

Original compilation and documentation use CC BY 4.0 only within the authors' rights; unchanged validator/tests use MIT. See LICENSE.md and LICENSES/file_licenses.json. No publisher permission grant is claimed. The exact-file historical notice and original license-scoped bytes are retained. The reader dictionary indicates modifications; the separate scientific snapshot is unchanged.

## Citation and reproducibility

See CITATION.cff and .zenodo.json for creator and version metadata. The validator and five synthetic tests use Python 3.10+ and the standard library only; no dependency installation is needed. Run both commands above from the extracted repository root, not an unrelated workspace. No participant-level data, coefficients, enrichment, pooling or clock scoring are supplied.

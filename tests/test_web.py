"""Cross-check browser outputs against the installed command-line workflow."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import zipfile

from eukasmchecker import cli

ROOT = Path(__file__).resolve().parents[1]


class BrowserParityTests(unittest.TestCase):
    def test_browser_outputs_match_cli(self):
        node = shutil.which("node")
        self.assertIsNotNone(node, "Node.js >=20 is needed for browser-edition tests")
        subprocess.run(["python3", str(ROOT / "web/build.py")], check=True)
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            subprocess.run([node, str(ROOT / "tests/web_workflow.cjs"), str(work)], check=True)
            # ECM input is normalized to the CLI's established fixed column names.
            with (work / "classification.tsv").open() as source:
                content = source.read().replace("id\tbp\tkingdom", "contig_id\tlength_bp\tlabel_euk", 1)
            (work / "cli_ecm.tsv").write_text(content)
            prepared = work / "sample_A.eac.json"
            cli.prepare(argparse.Namespace(sample="sample_A", fasta=str(work / "assembly.fasta"),
                ecm=str(work / "cli_ecm.tsv"), coverage=str(work / "coverage.tsv"),
                seq_id="identifier", seq_cov="mean_reads", out=str(prepared),
                scratch=work / "scratch", log_dir=None))
            manifest = cli.build_review(prepared, work / "review")
            meta = manifest["samples"][0]
            for case in json.loads((work / "cases.json").read_text()):
                with self.subTest(case=case["name"]):
                    decision = dict(sample_id="sample_A", fingerprint=meta["fingerprint"],
                        eac_path=meta["eac_path"], eac_sha256=meta["eac_sha256"],
                        reviewed=True, retain_sample=case["retain"], min_length_bp=500,
                        min_coverage=case["min"], max_coverage=case["max"], labels=case["labels"], excluded_seq_ids=case.get("excluded", []),advanced_rules=case.get("rules", []),
                        expected_kept_contigs=case["count"], expected_kept_bp=case["bp"],
                        mode="manual", note="", reviewed_at="synthetic")
                    document = dict(schema_version=9, dataset_id=manifest["dataset_id"], samples=[decision])
                    path = work / "decisions.json"
                    path.write_text(json.dumps(document))
                    output = work / ("cli_" + case["name"])
                    cli.apply(argparse.Namespace(decisions=str(path),out=str(output),
                        scratch=work / "scratch",log_dir=None,threads=1))
                    folder = output / "sample_A"
                    self.assertEqual((folder / "sample_A.fasta").read_bytes(), (work / (case["name"] + ".fasta")).read_bytes())
                    with (folder / "sample_A.sequences.tsv").open() as left, (work / (case["name"] + ".tsv")).open() as right:
                        a,b=list(csv.reader(left,delimiter="\t")),list(csv.reader(right,delimiter="\t"))
                    self.assertEqual(a[0],b[0])
                    self.assertEqual(len(a),len(b))
                    for x,y in zip(a[1:],b[1:]):
                        self.assertEqual(x[:2]+x[3:],y[:2]+y[3:])
                        self.assertEqual(float(x[2]),float(y[2]))
                    with zipfile.ZipFile(work / (case["name"] + ".zip")) as archive:
                        self.assertIsNone(archive.testzip())
                        self.assertEqual(archive.read("sample_A/sample_A.fasta"),(folder / "sample_A.fasta").read_bytes())
                        self.assertEqual(archive.read("sample_A/sample_A.sequences.tsv"),(work / (case["name"] + ".tsv")).read_bytes())

            with zipfile.ZipFile(work / "worker-manual.zip") as archive:
                saved = json.loads(archive.read("decisions.json"))
                self.assertEqual(saved["schema_version"], 8)
                self.assertEqual(saved["samples"][0]["excluded_seq_ids"], ["b"])
                self.assertEqual(archive.read("sample_A/sample_A.fasta"), (work / "manual.fasta").read_bytes())

            with zipfile.ZipFile(work / "worker-advanced.zip") as archive:
                saved = json.loads(archive.read("decisions.json"))
                self.assertEqual(saved["schema_version"], 9)
                self.assertEqual(len(saved["samples"][0]["advanced_rules"]), 2)
                self.assertEqual(archive.read("sample_A/sample_A.fasta"), (work / "advanced_multi.fasta").read_bytes())

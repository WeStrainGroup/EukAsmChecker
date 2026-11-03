#!/usr/bin/env python3
"""Synthetic CLI, output and local-save regression checks."""
import argparse
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
import threading
import unittest
from unittest.mock import patch
from urllib.request import Request, urlopen
from urllib.error import HTTPError
import zipfile

from eukasmchecker import cli as g
import sys


def fixture(root, sample="sample_A"):
    base = root / sample
    base.mkdir()
    (base / "assembly.fasta").write_bytes(b">a original header\r\n" + b"A" * 10001 + b"\r\n>b description\n" + b"C" * 500 + b"\n>c\n" + b"G" * 499 + b"\n")
    (base / "coverage.tsv").write_text("sequence\treads_cov\textra\nc\t100.5\tx\nb\t10\tx\na\t100\tx\n")
    (base / "ecm.tsv").write_text("contig_id\tlength_bp\tlabel_euk\nb\t500\tother\nc\t499\tNA\na\t10001\teukaryota\n")
    return argparse.Namespace(sample=sample, fasta=str(base / "assembly.fasta"),
        ecm=str(base / "ecm.tsv"), coverage=str(base / "coverage.tsv"), seq_id="sequence", seq_cov="reads_cov",
        out=str(base / (sample + ".eac.json")), scratch=root / "scratch", log_dir=None)


def decisions(manifest):
    entries = []
    for meta in manifest["samples"]:
        entries.append(dict(sample_id=meta["sample_id"],fingerprint=meta["fingerprint"],
            eac_path=meta["eac_path"],eac_sha256=meta["eac_sha256"],retain_sample=True,
            min_length_bp=500,min_coverage=10,max_coverage=None,labels=["eukaryota","other","NA"],
            reviewed=True,reviewed_at="fixture",note="Synthetic review",mode="default",
            expected_kept_contigs=2,expected_kept_bp=10501))
    return dict(schema_version=7,dataset_id=manifest["dataset_id"],samples=entries)


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.args = fixture(self.root)
        g.prepare(self.args)
        self.folder = self.root / "review"
        self.manifest = g.build_review(Path(self.args.out), self.folder)
        self.document = decisions(self.manifest)

    def tearDown(self):
        self.temp.cleanup()

    def apply(self, name="filtered"):
        path = self.root / "decisions.json"
        g.dump(path,self.document)
        output=self.root / name
        g.apply(argparse.Namespace(decisions=str(path),out=str(output),scratch=self.root/"scratch",log_dir=None))
        return output

    def test_prepare_mapping_and_precision(self):
        obj, checksum=g.load_sample(self.args.out)
        self.assertEqual(obj["rows"], [["a",10001,100.0,"eukaryota"],["b",500,10.0,"other"],["c",499,100.5,"NA"]])
        self.assertEqual(obj["fasta"],str(Path(self.args.fasta).resolve()))
        self.assertEqual(checksum,self.manifest["samples"][0]["eac_sha256"])
        self.assertEqual(self.manifest["samples"][0]["reference"]["min_coverage"],10)

    def test_apply_exact_two_files_and_all_reasons(self):
        d=self.document["samples"][0]
        d.update(max_coverage=99,labels=["other"],expected_kept_contigs=1,expected_kept_bp=500)
        output=self.apply()/"sample_A"
        self.assertEqual(sorted(p.name for p in output.iterdir()),["sample_A.fasta","sample_A.sequences.tsv"])
        self.assertEqual((output/"sample_A.fasta").read_bytes(),b">b description\n"+b"C"*500+b"\n")
        rows=list(g.table(output/"sample_A.sequences.tsv"))
        self.assertEqual(list(rows[0]),g.STATS)
        self.assertEqual(rows[0]["filter_reason"],"high_coverage;excluded_label")
        self.assertEqual(rows[2]["filter_reason"],"short_length;high_coverage;excluded_label")
        self.assertEqual(rows[1]["filter_status"],"retained")
        self.assertEqual(rows[2]["coverage"],"100.5")
        self.assertEqual(list(output.parent.iterdir()),[output])

    def test_manual_exclusions_and_restore(self):
        self.document["schema_version"] = 8
        d = self.document["samples"][0]
        d.update(excluded_seq_ids=["a", "c"], expected_kept_contigs=1, expected_kept_bp=500)
        output = self.apply("manual") / "sample_A"
        self.assertEqual((output / "sample_A.fasta").read_bytes(), b">b description\n" + b"C" * 500 + b"\n")
        rows = list(g.table(output / "sample_A.sequences.tsv"))
        self.assertEqual(rows[0]["filter_reason"], "manual_exclusion")
        self.assertEqual(rows[2]["filter_reason"], "short_length;manual_exclusion")
        d.update(excluded_seq_ids=[], expected_kept_contigs=2, expected_kept_bp=10501)
        restored = self.apply("restored") / "sample_A/sample_A.fasta"
        self.assertEqual(restored.read_bytes(), b">a original header\r\n" + b"A" * 10001 + b"\r\n>b description\n" + b"C" * 500 + b"\n")

    def test_invalid_manual_exclusions(self):
        self.document["schema_version"] = 8
        d = self.document["samples"][0]
        for excluded in [None, "a", ["a", "a"], [1], [""]]:
            with self.subTest(excluded=excluded):
                d["excluded_seq_ids"] = excluded
                with self.assertRaisesRegex(ValueError, "manual exclusions"):
                    g.validate_decisions(self.manifest, self.document)
        del d["excluded_seq_ids"]
        with self.assertRaisesRegex(ValueError, "Missing manual"):
            g.validate_decisions(self.manifest, self.document)
        d["excluded_seq_ids"] = ["missing"]
        with self.assertRaisesRegex(ValueError, "Unknown manually excluded"):
            self.apply("unknown_manual")

    def test_manual_exclusions_require_matching_review_counts(self):
        self.document["schema_version"] = 8
        self.document["samples"][0]["excluded_seq_ids"] = ["b"]
        with self.assertRaisesRegex(ValueError, "Retained counts differ"):
            self.apply("stale_manual_counts")

    def test_advanced_rules_filter_restore_and_manual_override(self):
        self.document["schema_version"] = 9
        d = self.document["samples"][0]
        first = dict(id=1, enabled=True, labels=["other"], length=dict(operator=">=", value=500),
                     coverage=dict(operator="<=", value=10))
        second = dict(id=2, enabled=True, labels=["NA"], length=dict(operator="any", value=0),
                      coverage=dict(operator=">", value=100))
        d.update(advanced_rules=[first, second], excluded_seq_ids=[], expected_kept_contigs=1, expected_kept_bp=10001)
        output = self.apply("advanced") / "sample_A"
        self.assertEqual((output / "sample_A.fasta").read_bytes(), b">a original header\r\n" + b"A" * 10001 + b"\r\n")
        rows = list(g.table(output / "sample_A.sequences.tsv"))
        self.assertEqual(rows[1]["filter_reason"], "advanced_rule_1")
        self.assertEqual(rows[2]["filter_reason"], "short_length;advanced_rule_2")
        first["enabled"] = False
        d.update(expected_kept_contigs=2, expected_kept_bp=10501)
        self.apply("disabled")
        first["enabled"] = True
        d.update(excluded_seq_ids=["a"], expected_kept_contigs=0, expected_kept_bp=0)
        self.assertEqual((self.apply("manual_advanced") / "sample_A/sample_A.fasta").read_bytes(), b"")
        d["retain_sample"] = False
        output = self.apply("excluded_advanced") / "sample_A"
        self.assertTrue(all(row["filter_reason"] == "excluded_sample" for row in g.table(output / "sample_A.sequences.tsv")))

    def test_advanced_boundaries_overlap_and_precision(self):
        rule = dict(id=1, enabled=True, labels=["other"], length=dict(operator=">", value=500),
                    coverage=dict(operator="any", value=0))
        row = ["b", 500, 10.4, "other"]
        for operator, expected in [("any", True), (">", False), (">=", True), ("<", False), ("<=", True)]:
            rule["length"]["operator"] = operator
            self.assertEqual(g.rule_matches(row, rule), expected)
        numeric_rule = dict(rule, id=1.0, length=dict(operator=">", value=500.0))
        validated = g.validate_advanced_rules([numeric_rule], {"other"})
        self.assertIs(type(validated[0]["id"]), int)
        self.assertIs(type(validated[0]["length"]["value"]), int)
        rule["length"]["operator"] = "any"
        rule["coverage"] = dict(operator=">", value=10)
        self.assertTrue(g.rule_matches(row, rule))
        d = self.document["samples"][0]
        d["advanced_rules"] = [rule, dict(rule, id=2)]
        self.assertEqual(g.reason(row, d), "advanced_rule_1;advanced_rule_2")
        rule["labels"] = ["eukaryota"]
        self.assertFalse(g.rule_matches(row, rule))

    def test_invalid_advanced_rules_and_legacy_import(self):
        import copy
        self.document["schema_version"] = 9
        d = self.document["samples"][0]
        d["excluded_seq_ids"] = []
        with self.assertRaisesRegex(ValueError, "Missing advanced"):
            g.validate_decisions(self.manifest, self.document)
        rule = dict(id=1, enabled=True, labels=["other"], length=dict(operator=">", value=500),
                    coverage=dict(operator="any", value=0))
        invalid = [None, "rule", [None], [dict(rule, id=True)], [rule, rule], [dict(rule, enabled=1)],
                   [dict(rule, labels=[])], [dict(rule, labels=["missing"])], [dict(rule, labels=["other", "other"])],
                   [dict(rule, length=dict(operator="=", value=500))],
                   [dict(rule, coverage=dict(operator=">", value=10.5))],
                   [dict(rule, coverage=dict(operator="any", value=True))],
                   [dict(rule, coverage=dict(operator="any", value=2**53))]]
        for value in invalid:
            with self.subTest(value=value):
                d["advanced_rules"] = copy.deepcopy(value)
                with self.assertRaisesRegex(ValueError, "advanced"):
                    g.validate_decisions(self.manifest, self.document)
        for version in (7, 8):
            self.document["schema_version"] = version
            d.pop("advanced_rules", None)
            restored = g.validate_decisions(self.manifest, self.document)
            self.assertEqual(restored["sample_A"]["advanced_rules"], [])

    def test_advanced_rules_require_matching_counts(self):
        self.document["schema_version"] = 9
        self.document["samples"][0].update(excluded_seq_ids=[], advanced_rules=[
            dict(id=1, enabled=True, labels=["other"], length=dict(operator="any", value=0), coverage=dict(operator="any", value=0))])
        with self.assertRaisesRegex(ValueError, "Retained counts differ"):
            self.apply("advanced_stale_counts")

    def test_default_maximum_and_inclusive_bound(self):
        d=self.document["samples"][0]
        d.update(min_length_bp=0,min_coverage=0,expected_kept_contigs=3,expected_kept_bp=11000)
        output=self.apply()/"sample_A/sample_A.fasta"
        self.assertEqual(output.read_bytes(),Path(self.args.fasta).read_bytes())
        d.update(min_length_bp=500,min_coverage=10,max_coverage=100,expected_kept_contigs=2,expected_kept_bp=10501)
        self.apply("bounded")
        self.assertEqual(g.reason(["x",500,100,"other"],d),"pass")
        self.assertEqual(g.reason(["x",500,100.01,"other"],d),"high_coverage")

    def test_excluded_sample(self):
        self.document["samples"][0].update(retain_sample=False,min_coverage=None,labels=[],expected_kept_contigs=0,expected_kept_bp=0)
        output=self.apply()/"sample_A"
        self.assertEqual((output/"sample_A.fasta").read_bytes(),b"")
        self.assertEqual({r["filter_reason"] for r in g.table(output/"sample_A.sequences.tsv")},{"excluded_sample"})

    def test_original_tables_not_required(self):
        Path(self.args.ecm).unlink()
        Path(self.args.coverage).unlink()
        self.apply()

    def test_unreviewed_no_output(self):
        self.document["samples"][0]["reviewed"]=False
        with self.assertRaisesRegex(ValueError,"Unreviewed samples: sample_A"):
            self.apply()
        self.assertFalse((self.root/"filtered.applying").exists())

    def test_changed_prepared_file(self):
        with Path(self.args.out).open("a") as stream:stream.write(" ")
        with self.assertRaisesRegex(ValueError,"Prepared sample changed"):
            self.apply()

    def test_changed_fasta(self):
        with Path(self.args.fasta).open("ab") as stream:stream.write(b"A")
        with self.assertRaisesRegex(ValueError,"Original FASTA changed"):
            self.apply()
        self.assertFalse((self.root/"filtered").exists())

    def test_invalid_thresholds_and_counts(self):
        for value in [-1,1.5,True,9]:
            with self.subTest(value=value):
                self.document["samples"][0]["max_coverage"]=value
                with self.assertRaises(ValueError):g.validate_decisions(self.manifest,self.document)
        self.document["samples"][0].update(max_coverage=None,expected_kept_bp=1)
        with self.assertRaisesRegex(ValueError,"Retained counts differ"):
            self.apply()

    def test_no_overwrite(self):
        original=Path(self.args.out).read_bytes()
        with self.assertRaisesRegex(ValueError,"already exists"):g.prepare(self.args)
        self.assertEqual(Path(self.args.out).read_bytes(),original)
        self.apply()
        with self.assertRaisesRegex(ValueError,"already exists"):self.apply()

    def test_input_errors(self):
        cases=[("coverage.tsv","sequence\treads_cov\na\t10\na\t20\n","Duplicate coverage"),
               ("coverage.tsv","sequence\treads_cov\na\t10\n","Missing coverage"),
               ("coverage.tsv","sequence\treads_cov\nx\t10\n","Extra or invalid coverage"),
               ("coverage.tsv","sequence\treads_cov\na\t-1\n","Invalid coverage"),
               ("coverage.tsv","sequence\treads_cov\na\tnan\n","Invalid coverage"),
               ("coverage.tsv","id\tvalue\na\t10\n","Coverage must contain"),
               ("ecm.tsv","contig_id\tlength_bp\tlabel_euk\na\t1\teuk\n","ECM length differs"),
               ("ecm.tsv","contig_id\tlength_bp\tlabel_euk\na\t10001\teuk\na\t10001\teuk\n","Duplicate ECM"),
               ("ecm.tsv","contig_id\tlength_bp\tlabel_euk\na\t10001\teuk\n","Missing ECM")]
        for index,(filename,content,message) in enumerate(cases):
            with self.subTest(case=index):
                args=fixture(self.root,"bad"+str(index))
                Path(args.out).parent.joinpath(filename).write_text(content)
                with self.assertRaisesRegex(ValueError,message):g.prepare(args)
        args=fixture(self.root,"duplicate_fasta")
        with Path(args.fasta).open("ab") as stream:stream.write(b">a\nA\n")
        with self.assertRaisesRegex(ValueError,"Duplicate FASTA"):g.prepare(args)

    def test_recursive_sort_and_conflicts(self):
        args=fixture(self.root,"aaa")
        g.prepare(args)
        manifest=g.build_review(self.root,self.root/"sorted")
        self.assertEqual([m["sample_id"] for m in manifest["samples"]],["aaa","sample_A"])
        copy=self.root/"copy.eac.json"
        copy.write_bytes(Path(self.args.out).read_bytes())
        with self.assertRaisesRegex(ValueError,"Duplicate sample IDs") as error:
            g.build_review(self.root,self.root/"duplicate")
        self.assertIn(str(copy),str(error.exception))
        self.assertIn(self.args.out,str(error.exception))

    def test_zip_and_cli(self):
        output=self.root/"bundle.zip"
        result=subprocess.run([sys.executable,"-m","eukasmchecker","review",self.args.out,"--out",str(output)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        with zipfile.ZipFile(output) as bundle:
            self.assertEqual(sorted(bundle.namelist()),["data/000000.js","index.html"])
            html=bundle.read("index.html").decode()
            self.assertIn(self.args.out,html)
            self.assertIn("Save decisions",html)
            self.assertNotIn("/*PROJECT_MANIFEST*/null",html)
            self.assertNotIn('src="https://',html)
        with patch.object(g,"browser_available",return_value=False):
            with self.assertRaisesRegex(ValueError,"--out"):
                g.review(argparse.Namespace(input=self.args.out,out=None,scratch=self.root/"scratch",log_dir=None))

    def test_legacy_migration_validation(self):
        for version in (2,3,4,5,6):
            self.document["schema_version"]=version
            self.document["dataset_id"]="old-dataset"
            migrated=g.validate_decisions(self.manifest,self.document)
            self.assertIsNone(migrated["sample_A"]["max_coverage"])
        with self.assertRaisesRegex(ValueError,"Legacy decisions"):
            self.apply()

    def test_parallel_matches_serial(self):
        second = fixture(self.root, "sample_B")
        g.prepare(second)
        serial = g.build_review(self.root, self.root / "serial_review", threads=1)
        parallel = g.build_review(self.root, self.root / "parallel_review", threads=2)
        self.assertEqual(serial, parallel)
        for meta in serial["samples"]:
            relative = meta["data_file"]
            self.assertEqual((self.root / "serial_review" / relative).read_bytes(),
                             (self.root / "parallel_review" / relative).read_bytes())
        self.document = decisions(serial)
        first = self.apply("serial_filtered")
        path = self.root / "decisions.json"
        second_output = self.root / "parallel_filtered"
        g.apply(argparse.Namespace(decisions=str(path), out=str(second_output),
                                  scratch=self.root / "scratch", log_dir=None, threads=2))
        for file in first.rglob("*"):
            if file.is_file():
                self.assertEqual(file.read_bytes(), (second_output / file.relative_to(first)).read_bytes())

    def test_apply_beside_prepared_sample(self):
        prepared = Path(self.args.out)
        original = prepared.read_bytes()
        # The root already contains sample_A/sample_A.eac.json.
        path = self.root / "decisions.json"
        g.dump(path, self.document)
        g.apply(argparse.Namespace(decisions=str(path), out=str(self.root),
                                  scratch=self.root / "scratch", log_dir=None, threads=2))
        self.assertEqual(prepared.read_bytes(), original)
        self.assertTrue(prepared.with_name("sample_A.fasta").is_file())
        self.assertTrue(prepared.with_name("sample_A.sequences.tsv").is_file())
        self.assertFalse(self.root.with_name(self.root.name + ".applying").exists())

    def test_publication_failure_rolls_back(self):
        partial = self.root / "staged_results"
        output = self.root / "existing"
        sample = partial / "sample_A"
        sample.mkdir(parents=True)
        target = output / "sample_A"
        target.mkdir(parents=True)
        marker = target / "notes.txt"
        marker.write_text("preserve")
        for suffix in (".fasta", ".sequences.tsv"):
            (sample / ("sample_A" + suffix)).write_text("verified result")
        real_link = g.os.link
        count = 0
        def failing_link(source, destination):
            nonlocal count
            count += 1
            if count == 2:
                raise OSError("Synthetic publication failure")
            real_link(source, destination)
        with patch.object(g.os, "link", failing_link):
            with self.assertRaisesRegex(OSError, "Synthetic publication failure"):
                g.publish_results(partial, output, self.manifest["samples"], True)
        self.assertEqual(list(target.iterdir()), [marker])
        self.assertEqual(marker.read_text(), "preserve")
        self.assertEqual(len(list(sample.iterdir())), 2)
        g.publish_results(partial, output, self.manifest["samples"], True)
        self.assertTrue((target / "sample_A.fasta").exists())

    def test_output_source_mismatch(self):
        output = self.root / "other_prepared"
        folder = output / "sample_A"
        folder.mkdir(parents=True)
        (folder / "sample_A.eac.json").write_bytes(Path(self.args.out).read_bytes())
        with self.assertRaisesRegex(ValueError, "different prepared sample"):
            g.check_destinations(output, self.manifest["samples"])

    def test_local_server_save_conflict_and_restart(self):
        progress=self.root/"progress/decisions.json"
        instances=[]
        ready=threading.Event()
        class Server(g.HTTPServer):
            def server_activate(self):
                super().server_activate()
                instances.append(self)
        def browser(url,new):
            ready.set()
            return True
        errors=[]
        def run():
            try:g.serve_review(self.manifest,self.folder,progress)
            except Exception as exc:errors.append(exc)
        with patch.object(g,"HTTPServer",Server),patch.object(g.webbrowser,"open",browser):
            thread=threading.Thread(target=run, daemon=True)
            thread.start()
            self.assertTrue(ready.wait(60), errors)
            server=instances[0]
            base=f"http://127.0.0.1:{server.server_port}"
            token=self.manifest["review_server"]["token"]
            def request(document=None,revision="",auth=True):
                data=json.dumps(dict(decisions=document,revision=revision)).encode() if document else None
                headers={"X-EAC-Token":token if auth else "bad","Content-Type":"application/json"}
                return json.load(urlopen(Request(base+"/api/decisions",data=data,headers=headers),timeout=5))
            try:
                self.assertIsNone(request()["decisions"])
                with self.assertRaises(HTTPError) as error:request(auth=False)
                self.assertEqual(error.exception.code,403)
                saved=request(self.document)
                self.assertEqual(json.loads(progress.read_text()),self.document)
                with self.assertRaises(HTTPError) as error:request(self.document)
                self.assertEqual(error.exception.code,409)
                self.document["samples"][0].update(note="Server note",reviewed=False,expected_kept_contigs=None,expected_kept_bp=None)
                request(self.document,saved["revision"])
                self.assertEqual(request()["decisions"]["samples"][0]["note"],"Server note")
            finally:
                server.shutdown()
                thread.join(5)
        self.assertFalse(errors,errors)
        self.assertFalse((progress.parent/"review.lock").exists())
        with patch.object(g.webbrowser,"open",return_value=False):
            with self.assertRaisesRegex(ValueError,"Browser launch request failed"):
                g.serve_review(self.manifest,self.folder,progress)
        self.assertFalse((progress.parent/"review.lock").exists())


if __name__ == "__main__":
    unittest.main()

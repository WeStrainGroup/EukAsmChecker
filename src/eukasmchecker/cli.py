#!/usr/bin/env python3
"""Prepare, review and filter eukaryotic genome assemblies."""
import argparse
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import secrets
import shlex
import shutil
import subprocess
import sys
import tempfile
import webbrowser
import zipfile
from datetime import datetime, timezone
from http.server import HTTPServer, SimpleHTTPRequestHandler
from urllib.parse import urlsplit
from concurrent.futures import ProcessPoolExecutor
from . import __version__

VERSION = __version__
COLUMNS = ["seq_id", "length_bp", "coverage", "ecm_label"]
STATS = COLUMNS + ["filter_status", "filter_reason"]

def require(condition, message):
    if not condition:
        raise ValueError(message)


def now():
    return datetime.now(timezone.utc).isoformat()


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=True, allow_nan=False, indent=2) + "\n")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def table(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)),
                f"Missing or duplicate headers: {path}")
        for row in reader:
            require(None not in row and None not in row.values(), f"Malformed row: {path}")
            yield row


def number(value, field, integer=False):
    require(value is not None and str(value).strip() != "", f"Missing {field}")
    require(not isinstance(value, bool), f"Invalid {field}: {value}")
    v = float(value)
    require(math.isfinite(v) and v >= 0, f"Invalid {field}: {value}")
    if integer:
        require(v.is_integer() and v <= 9007199254740991, f"Invalid integer {field}: {value}")
        return int(v)
    return v


def fasta_records(path):
    header, lines, cid, length = None, [], None, 0
    with Path(path).open("rb") as f:
        for line in f:
            if line.startswith(b">"):
                if header is not None:
                    yield cid, length, header, lines
                header, lines, length = line, [], 0
                name = line[1:].split()
                require(name, f"Empty FASTA header: {path}")
                cid = name[0].decode("utf-8")
            else:
                require(header is not None, f"Sequence before FASTA header: {path}")
                seq = line.rstrip(b"\r\n")
                require(not re.search(rb"[ \t\v\f]", seq), f"Whitespace within sequence: {cid}")
                length += len(seq)
                lines.append(line)
        if header is not None:
            yield cid, length, header, lines


def check_fasta(path, rows, output=None, keep=None):
    expected = {r[0]: r[1] for r in rows}
    seen = set()
    for cid, length, header, lines in fasta_records(path):
        require(cid in expected and cid not in seen, f"FASTA extra/duplicate ID: {cid}")
        require(length == expected[cid], f"FASTA length mismatch: {cid}")
        seen.add(cid)
        if output is not None and cid in keep:
            output.write(header)
            output.writelines(lines)
    require(seen == set(expected), "FASTA missing sequence IDs.")


def reference(rows):
    depths = sorted(r[2] for r in rows if r[1] > 10000 and r[3] == "eukaryota")
    n = len(depths)
    median = depths[n // 2] / 2 + depths[(n - 1) // 2] / 2 if n else None
    return dict(n=n, cov_10kb=median,
                min_coverage=math.floor(max(median / 25, 10) + 0.5) if median is not None else None)


def validate_advanced_rules(rules, available):
    require(isinstance(rules, list), "Invalid advanced rules.")
    ids = set()
    for rule in rules:
        require(isinstance(rule, dict), "Invalid advanced rule.")
        rule_id = rule.get("id")
        require(type(rule_id) in (int, float) and 0 < rule_id <= 2**53 - 1 and int(rule_id) == rule_id and rule_id not in ids,
                "Invalid advanced rule ID.")
        rule["id"] = int(rule_id)
        ids.add(rule_id)
        require(isinstance(rule.get("enabled"), bool), "Invalid advanced rule state.")
        labels = rule.get("labels")
        require(isinstance(labels, list) and bool(labels) and all(isinstance(label, str) for label in labels)
                and len(labels) == len(set(labels)) and set(labels) <= available,
                "Invalid advanced rule labels.")
        for field in ("length", "coverage"):
            condition = rule.get(field)
            require(isinstance(condition, dict) and condition.get("operator") in ("any", ">", ">=", "<", "<="),
                    "Invalid advanced rule operator.")
            value = condition.get("value")
            require(type(value) in (int, float) and 0 <= value <= 2**53 - 1 and int(value) == value,
                    "Invalid advanced rule threshold.")
            condition["value"] = int(value)
    return rules


def condition_matches(value, condition):
    operator, threshold = condition["operator"], condition["value"]
    if operator == "any":
        return True
    if operator == ">":
        return value > threshold
    if operator == ">=":
        return value >= threshold
    if operator == "<":
        return value < threshold
    return value <= threshold


def rule_matches(row, rule):
    return (row[3] in rule["labels"] and condition_matches(row[1], rule["length"])
            and condition_matches(row[2], rule["coverage"]))


def reason(row, decision):
    if decision.get("retain_sample", True) is False:
        return "excluded_sample"
    reasons = []
    if row[1] < decision["min_length_bp"]:
        reasons.append("short_length")
    if row[2] < decision["min_coverage"]:
        reasons.append("low_coverage")
    if decision.get("max_coverage") is not None and row[2] > decision["max_coverage"]:
        reasons.append("high_coverage")
    if decision["labels"] != "all" and row[3] not in (decision["labels"] if isinstance(decision["labels"], list) else [decision["labels"]]):
        reasons.append("excluded_label")
    if row[0] in decision.get("excluded_seq_ids", ()):
        reasons.append("manual_exclusion")
    for rule in decision.get("advanced_rules", []):
        if rule["enabled"] and rule_matches(row, rule):
            reasons.append("advanced_rule_" + str(rule["id"]))
    return ";".join(reasons) or "pass"


def write_tsv(path, columns, rows):
    with Path(path).open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows([["NA" if v is None else ";".join(v) if isinstance(v, list) else v for v in row] for row in rows])



def sample_id(value):
    require(isinstance(value, str) and re.fullmatch(r"[A-Za-z0-9_.-]+", value)
            and value not in (".", ".."), f"Invalid sample ID: {value}")
    return value


def copy_input(source, destination):
    subprocess.run(["cp", "--", str(source), str(destination)], check=True)


def scratch(args, action):
    args.scratch.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix="eukasmchecker_" + action + "_", dir=args.scratch))


def publish_file(source, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".eac_", dir=output.parent) as temporary:
        staged = Path(temporary) / output.name
        copy_input(source, staged)
        require(sha(staged) == sha(source), f"Copyback mismatch: {output}")
        # Hard-link publication fails if another process created the destination.
        os.link(staged, output)


def log_result(args, action, **info):
    if args.log_dir is not None:
        args.log_dir.mkdir(parents=True, exist_ok=True)
        dump(args.log_dir / (datetime.now().strftime("%Y%m%d_%H%M%S_%f") + "_" + action + ".json"),
             dict(version=VERSION, action=action, time=now(), **info))


def read_inputs(files, id_column, cov_column):
    lengths = {}
    for cid, length, _, _ in fasta_records(files["fasta"]):
        require(cid not in lengths, f"Duplicate FASTA ID: {cid}")
        require(length > 0, f"Empty FASTA sequence: {cid}")
        lengths[cid] = length
    require(lengths, "FASTA contains no sequences.")
    coverage = {}
    for row in table(files["coverage"]):
        require(id_column in row and cov_column in row,
                f"Coverage must contain columns {id_column!r} and {cov_column!r}.")
        cid = row[id_column]
        require(cid in lengths, f"Extra or invalid coverage ID: {cid}")
        require(cid not in coverage, f"Duplicate coverage ID: {cid}")
        coverage[cid] = number(row[cov_column], f"coverage for {cid}")
    missing = lengths.keys() - coverage.keys()
    require(not missing, "Missing coverage IDs: " + ", ".join(sorted(missing)))
    labels = {}
    for row in table(files["ecm"]):
        require(all(k in row for k in ("contig_id", "length_bp", "label_euk")),
                "ECM requires contig_id, length_bp and label_euk columns.")
        cid = row["contig_id"]
        require(cid in lengths, f"Extra or invalid ECM ID: {cid}")
        require(cid not in labels, f"Duplicate ECM ID: {cid}")
        require(number(row["length_bp"], "ECM length", True) == lengths[cid],
                f"ECM length differs from FASTA: {cid}")
        label = row["label_euk"]
        require(label and not re.search(r"\s|;", label), f"Invalid ECM label: {cid}")
        labels[cid] = "eukaryota" if label == "euk" else label
    missing = lengths.keys() - labels.keys()
    require(not missing, "Missing ECM IDs: " + ", ".join(sorted(missing)))
    return [[cid, length, coverage[cid], labels[cid]] for cid, length in lengths.items()]


def prepare(args):
    sample = sample_id(args.sample)
    output = Path(args.out).resolve()
    require(output.name.endswith(".eac.json"), "prepare --out must end with .eac.json")
    require(not output.exists(), f"Output already exists: {output}")
    inputs = {k: Path(getattr(args, k)).resolve() for k in ("fasta", "ecm", "coverage")}
    task = scratch(args, "prepare")
    try:
        files = {}
        for kind, source in inputs.items():
            require(source.is_file(), f"Missing {kind}: {source}")
            files[kind] = task / kind
            copy_input(source, files[kind])
        rows = read_inputs(files, args.seq_id, args.seq_cov)
        hashes = {k: sha(v) for k, v in files.items()}
        obj = dict(format="EukAsmChecker sample", schema_version=1, tool_version=VERSION,
                   sample_id=sample, created_at=now(), fasta=str(inputs["fasta"]),
                   fasta_sha256=hashes["fasta"], input_sha256=hashes,
                   fingerprint=digest(hashes), columns=COLUMNS, rows=rows)
        dump(task / "sample.eac.json", obj)
        publish_file(task / "sample.eac.json", output)
        log_result(args, "prepare", sample_id=sample, output=str(output), sequences=len(rows))
    except Exception:
        print(f"Preparation failed. Scratch retained at {task}", file=sys.stderr)
        raise
    shutil.rmtree(task)
    print(f"Prepared {sample}: {output}", flush=True)


def load_sample(path, expected_sha=None):
    path = Path(path).resolve()
    raw = path.read_bytes()
    checksum = hashlib.sha256(raw).hexdigest()
    require(expected_sha is None or checksum == expected_sha, f"Prepared sample changed: {path}")
    obj = json.loads(raw)
    require(obj.get("format") == "EukAsmChecker sample" and obj.get("schema_version") == 1,
            f"Unsupported .eac.json: {path}")
    sample_id(obj.get("sample_id"))
    require(obj.get("columns") == COLUMNS, f"Invalid sequence columns: {path}")
    require(Path(obj["fasta"]).is_absolute(), f"FASTA path must be absolute: {path}")
    require(obj["fasta_sha256"] == obj["input_sha256"]["fasta"]
            and obj["fingerprint"] == digest(obj["input_sha256"]), f"Invalid source fingerprint: {path}")
    require(all(re.fullmatch(r"[a-f0-9]{64}", h) for h in obj["input_sha256"].values()),
            f"Invalid input checksum: {path}")
    rows = obj["rows"]
    require(isinstance(rows, list) and rows, f"No sequences: {path}")
    seen = set()
    for row in rows:
        require(isinstance(row, list) and len(row) == 4, f"Malformed sequence row: {path}")
        cid, length, cov, label = row
        require(isinstance(cid, str) and cid and not re.search(r"\s", cid) and cid not in seen,
                f"Invalid or duplicate sequence ID: {cid}")
        seen.add(cid)
        require(number(length, "sequence length", True) > 0, f"Empty sequence: {cid}")
        number(cov, "coverage")
        require(isinstance(length, int) and isinstance(cov, (int, float)), "Sequence measurements must be numeric.")
        require(isinstance(label, str) and label and not re.search(r"\s|;", label), "Invalid ECM label.")
    return obj, checksum


def metadata(obj, path, checksum, index):
    rows = obj["rows"]
    groups = {}
    for _, length, _, label in rows:
        group = groups.setdefault(label, dict(contigs=0, bp=0))
        group["contigs"] += 1
        group["bp"] += length
    return dict(sample_id=obj["sample_id"], fingerprint=obj["fingerprint"],
                eac_path=str(path), eac_sha256=checksum, index=index,
                data_file=f"data/{index:06d}.js", contigs=len(rows),
                bp=sum(r[1] for r in rows), groups=groups, reference=reference(rows))


def run_tasks(function, tasks, threads):
    if threads == 1:
        return [function(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=min(threads, len(tasks))) as pool:
        return list(pool.map(function, tasks))


def review_sample(task):
    path, destination, index = task
    staged = destination / f"sample_{index:06d}.eac.json"
    copy_input(path, staged)
    obj, checksum = load_sample(staged)
    staged.unlink()
    meta = metadata(obj, path, checksum, index)
    payload = dict(sample_id=obj["sample_id"], fingerprint=obj["fingerprint"], rows=obj["rows"])
    (destination / meta["data_file"]).write_text("window.receiveGenomeSample(" +
        json.dumps(payload, ensure_ascii=True, allow_nan=False, separators=(",", ":")) + ");\n")
    return meta


def build_review(source, destination, threads=1):
    paths = [source] if source.is_file() else sorted(source.rglob("*.eac.json"))
    require(paths, f"No .eac.json files: {source}")
    (destination / "data").mkdir(parents=True)
    tasks = [(path, destination, index) for index, path in enumerate(paths)]
    samples = run_tasks(review_sample, tasks, threads)
    conflicts = {}
    for meta in samples:
        conflicts.setdefault(meta["sample_id"], []).append(meta["eac_path"])
    duplicates = [name + ": " + ", ".join(paths) for name, paths in conflicts.items() if len(paths) > 1]
    require(not duplicates, "Duplicate sample IDs:\n" + "\n".join(duplicates))
    samples.sort(key=lambda m: m["sample_id"])
    identity = digest([[s["sample_id"], s["fingerprint"], s["eac_sha256"], s["eac_path"]] for s in samples])
    return dict(schema_version=3, tool_version=VERSION, dataset_id=identity,
                project_name=source.stem.replace(".eac", "") if source.is_file() else source.name,
                columns=COLUMNS, samples=samples)


def render_review(manifest, destination):
    template = Path(__file__).with_name("review.html").read_text()
    require(template.count("/*PROJECT_MANIFEST*/null") == 1, "Missing manifest template marker.")
    html = template.replace("/*PROJECT_MANIFEST*/null", json.dumps(manifest, ensure_ascii=True,
                             allow_nan=False).replace("<", "\\u003c"))
    (destination / "index.html").write_text(html)


def validate_decisions(manifest, decisions, complete=True):
    version = decisions.get("schema_version")
    require(version in (2, 3, 4, 5, 6, 7, 8, 9), "Unsupported decisions schema.")
    if version >= 7:
        require(decisions.get("dataset_id") == manifest["dataset_id"], "Decisions belong to another dataset.")
    entries = decisions.get("samples")
    require(isinstance(entries, list), "Missing decisions.samples.")
    names = [d.get("sample_id") for d in entries]
    expected = [s["sample_id"] for s in manifest["samples"]]
    require(len(names) == len(set(names)) and set(names) == set(expected), "Decision sample set differs from inputs.")
    result = {d["sample_id"]: copy.deepcopy(d) for d in entries}
    pending = [name for name in expected if result[name].get("reviewed") is not True]
    require(not complete or not pending, "Unreviewed samples: " + ", ".join(pending))
    available = {"eukaryota", "other", "NA"} | {k for m in manifest["samples"] for k in m["groups"]}
    for meta in manifest["samples"]:
        d = result[meta["sample_id"]]
        require(d.get("fingerprint") == meta["fingerprint"], f"Input mismatch: {meta['sample_id']}")
        if version >= 7:
            require(d.get("eac_path") == meta["eac_path"] and d.get("eac_sha256") == meta["eac_sha256"],
                    f"Prepared sample reference mismatch: {meta['sample_id']}")
        require(isinstance(d.get("reviewed"), bool), "reviewed must be boolean.")
        require(isinstance(d.get("note", ""), str), "Note must be text.")
        d["retain_sample"] = d.get("retain_sample") if version >= 5 else True
        require(isinstance(d["retain_sample"], bool), "retain_sample must be boolean.")
        labels = d.get("labels")
        if isinstance(labels, str):
            labels = sorted(available) if labels == "all" else [labels]
        require(isinstance(labels, list) and all(isinstance(k, str) for k in labels)
                and len(labels) == len(set(labels)) and set(labels) <= available, "Invalid ECM label selection.")
        d["labels"] = labels
        if version >= 8:
            require("excluded_seq_ids" in d, "Missing manual exclusions.")
        excluded = d.get("excluded_seq_ids", [])
        require(isinstance(excluded, list) and all(isinstance(value, str) and value for value in excluded)
                and len(excluded) == len(set(excluded)), "Invalid manual exclusions.")
        d["excluded_seq_ids"] = excluded
        require(version < 9 or "advanced_rules" in d, "Missing advanced rules.")
        d["advanced_rules"] = validate_advanced_rules(d.get("advanced_rules", []), available)
        for field in ("min_length_bp", "min_coverage"):
            if d.get(field) is not None:
                d[field] = number(d[field], field, field == "min_length_bp" or version >= 4)
        if version < 6:
            d["max_coverage"] = None
        else:
            require("max_coverage" in d, "Missing max_coverage.")
            if d["max_coverage"] is not None:
                d["max_coverage"] = number(d["max_coverage"], "max_coverage", True)
        if d["reviewed"]:
            if d["retain_sample"]:
                require(d.get("min_length_bp") is not None and d.get("min_coverage") is not None and labels,
                        f"Reviewed sample has incomplete filters: {meta['sample_id']}")
                require(d["max_coverage"] is None or d["max_coverage"] >= d["min_coverage"],
                        "Maximum coverage must be at least minimum coverage.")
            for field in ("expected_kept_contigs", "expected_kept_bp"):
                d[field] = number(d.get(field), field, True)
            require(d["retain_sample"] or (d["expected_kept_contigs"] == 0 and d["expected_kept_bp"] == 0),
                    "Excluded samples must have zero retained counts.")
    return result


def browser_available():
    return sys.platform in ("darwin", "win32") or bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def serve_review(manifest, folder, progress):
    progress.parent.mkdir(parents=True, exist_ok=True)
    lock = progress.parent / "review.lock"
    require(not lock.exists(), f"Another review session owns {lock}. Close it first. If it crashed, remove this lock after checking the process.")
    token = secrets.token_urlsafe(32)
    current = None
    if progress.exists():
        current = json.loads(progress.read_text())
        validate_decisions(manifest, current, complete=False)
    revision = sha(progress) if progress.exists() else ""
    lock.mkdir()

    class Handler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(folder), **kwargs)

        def log_message(self, *args):
            pass

        def allowed(self):
            return self.headers.get("Host") == f"127.0.0.1:{self.server.server_port}"

        def reply(self, code, data):
            body = json.dumps(data, ensure_ascii=True).encode()
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            if not self.allowed():
                self.reply(403, {"error": "Invalid host"})
            elif urlsplit(self.path).path == "/api/decisions":
                if self.headers.get("X-EAC-Token") != token:
                    self.reply(403, {"error": "Invalid session"})
                else:
                    self.reply(200, dict(decisions=current, revision=revision, path=str(progress)))
            elif urlsplit(self.path).path in ("/", "/index.html") or re.fullmatch(r"/data/\d{6}\.js", urlsplit(self.path).path):
                super().do_GET()
            else:
                self.reply(404, {"error": "Not found"})

        def do_POST(self):
            nonlocal current, revision
            origin = self.headers.get("Origin")
            if (not self.allowed() or self.path != "/api/decisions" or self.headers.get("X-EAC-Token") != token
                    or origin not in (None, f"http://127.0.0.1:{self.server.server_port}")):
                self.reply(403, {"error": "Invalid review session"})
                return
            try:
                size = int(self.headers.get("Content-Length", "0"))
                require(0 < size <= 64 * 1024 * 1024, "Invalid save request size.")
                request = json.loads(self.rfile.read(size))
                if request.get("revision") != revision:
                    self.reply(409, {"error": "Decisions changed in another tab. Reload before continuing."})
                    return
                document = request["decisions"]
                validate_decisions(manifest, document, complete=False)
                temporary = progress.with_suffix(".tmp")
                dump(temporary, document)
                temporary.replace(progress)
                revision, current = sha(progress), document
                self.reply(200, dict(revision=revision, path=str(progress)))
            except (ValueError, KeyError, TypeError, OSError) as exc:
                self.reply(400, {"error": str(exc)})

    try:
        with HTTPServer(("127.0.0.1", 0), Handler) as server:
            manifest["review_server"] = dict(token=token, path=str(progress))
            render_review(manifest, folder)
            url = f"http://127.0.0.1:{server.server_port}/index.html"
            try:
                requested = webbrowser.open(url, new=2)
            except webbrowser.Error:
                requested = False
            require(requested, "Browser launch request failed.")
            print(f"Browser launch requested. Open this address on this computer if needed: {url}", flush=True)
            print(f"Decisions are saved to {progress}\nKeep this command running. Press Ctrl+C to stop.", flush=True)
            try:
                server.serve_forever()
            except KeyboardInterrupt:
                pass
    finally:
        lock.rmdir()


def review(args):
    source = Path(args.input).resolve()
    require(source.exists() and (source.is_dir() or source.name.endswith(".eac.json")),
            "review expects a .eac.json file or directory.")
    fallback = "eukasmchecker review " + shlex.quote(str(source)) + " --out " + shlex.quote(str(source.parent / "review.zip"))
    if not args.out:
        require(browser_available(), "No graphical browser environment. Create an offline bundle instead:\n" + fallback)
    task = scratch(args, "review")
    try:
        manifest = build_review(source, task, getattr(args, "threads", 1))
        if args.out:
            output = Path(args.out).resolve()
            require(output.suffix.lower() == ".zip", "review --out must end with .zip")
            require(not output.exists(), f"Output already exists: {output}")
            render_review(manifest, task)
            archive = task / "review.zip"
            with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as bundle:
                bundle.write(task / "index.html", "index.html")
                for path in sorted((task / "data").iterdir()):
                    bundle.write(path, "data/" + path.name)
            publish_file(archive, output)
            print(f"Offline review: {output}\nExtract the entire archive and open index.html. Use Save decisions to download decisions.json.", flush=True)
            log_result(args, "review", output=str(output), samples=len(manifest["samples"]))
        else:
            parent = source.parent if source.is_file() else source
            progress = parent / ".eukasmchecker" / manifest["dataset_id"] / "decisions.json"
            try:
                serve_review(manifest, task, progress)
            except (ValueError, OSError) as exc:
                raise ValueError(str(exc) + "\nOffline alternative:\n" + fallback) from exc
    finally:
        shutil.rmtree(task)


def inspect_prepared(task):
    decision, scratch_root, index = task
    name = sample_id(decision.get("sample_id"))
    path = Path(decision["eac_path"])
    require(path.is_absolute(), "Prepared sample references must be absolute paths.")
    work = scratch_root / name
    work.mkdir()
    staged = work / "sample.eac.json"
    copy_input(path, staged)
    obj, checksum = load_sample(staged, decision["eac_sha256"])
    return metadata(obj, path, checksum, index)


def filter_sample(task):
    meta, decision, scratch_root, partial = task
    name = meta["sample_id"]
    work = scratch_root / name
    obj, _ = load_sample(work / "sample.eac.json", meta["eac_sha256"])
    fasta = work / "original.fasta"
    copy_input(obj["fasta"], fasta)
    require(sha(fasta) == obj["fasta_sha256"], f"Original FASTA changed: {name}")
    rows = obj["rows"]
    excluded = set(decision.get("excluded_seq_ids", []))
    require(excluded <= {row[0] for row in rows}, f"Unknown manually excluded sequence: {name}")
    decision = dict(decision, excluded_seq_ids=excluded)
    reasons = [reason(row, decision) for row in rows]
    keep = {row[0] for row, why in zip(rows, reasons) if why == "pass"}
    require(len(keep) == decision["expected_kept_contigs"] and
            sum(row[1] for row in rows if row[0] in keep) == decision["expected_kept_bp"],
            f"Retained counts differ from review: {name}")
    result = work / "result"
    result.mkdir()
    with (result / (name + ".fasta")).open("wb") as stream:
        check_fasta(fasta, rows, stream, keep)
    write_tsv(result / (name + ".sequences.tsv"), STATS,
              [row + ["retained" if why == "pass" else "removed", why]
               for row, why in zip(rows, reasons)])
    destination = partial / name
    destination.mkdir()
    for path in result.iterdir():
        copy_input(path, destination / path.name)
        require(sha(path) == sha(destination / path.name), f"Copyback mismatch: {name}")
    # Retain scratch until every sample has been published successfully.
    return dict(sample_id=name, retained=len(keep), total=len(rows))


def check_destinations(output, samples):
    require(not output.exists() or output.is_dir(), f"Output is not a directory: {output}")
    for meta in samples:
        name = meta["sample_id"]
        folder = output / name
        require(not folder.is_symlink(), f"Sample output directory is a symlink: {folder}")
        require(not folder.exists() or folder.is_dir(), f"Sample output is not a directory: {folder}")
        prepared = folder / (name + ".eac.json")
        if prepared.exists():
            require(prepared.resolve() == Path(meta["eac_path"]).resolve(),
                    f"Output folder contains a different prepared sample: {folder}")
        for suffix in (".fasta", ".sequences.tsv"):
            target = folder / (name + suffix)
            require(not target.exists() and not target.is_symlink(), f"Output already exists: {target}")


def publish_results(partial, output, samples, existing_output):
    check_destinations(output, samples)
    if not existing_output:
        require(not output.exists(), f"Output already exists: {output}")
        partial.rename(output)
        return
    published = []
    created_folders = []
    try:
        for meta in samples:
            name = meta["sample_id"]
            folder = output / name
            if not folder.exists():
                folder.mkdir()
                created_folders.append(folder)
            for suffix in (".fasta", ".sequences.tsv"):
                source = partial / name / (name + suffix)
                target = folder / source.name
                # Same-filesystem links publish verified results without overwriting files.
                os.link(source, target)
                published.append((source, target))
    except Exception:
        for source, target in reversed(published):
            if target.exists() and os.path.samestat(source.stat(), target.stat()):
                target.unlink()
        for folder in reversed(created_folders):
            if not any(folder.iterdir()):
                folder.rmdir()
        raise
    shutil.rmtree(partial)


def apply(args):
    decisions_path = Path(args.decisions).resolve()
    output = Path(args.out).resolve()
    partial = output.with_name(output.name + ".applying")
    require(not partial.exists(), f"Incomplete output already exists: {partial}")
    existing_output = output.exists()
    task = scratch(args, "apply")
    threads = getattr(args, "threads", 1)
    try:
        staged_decisions = task / "decisions.json"
        copy_input(decisions_path, staged_decisions)
        document = json.loads(staged_decisions.read_text())
        require(document.get("schema_version") in (7, 8, 9),
                "Legacy decisions need source references. Import them into a review of the same prepared samples, then Save decisions and retry.")
        entries = document.get("samples")
        require(isinstance(entries, list) and entries, "No sample decisions.")
        names = [sample_id(d.get("sample_id")) for d in entries]
        require(len(names) == len(set(names)), "Duplicate sample decisions.")
        pending = [d["sample_id"] for d in entries if d.get("reviewed") is not True]
        require(not pending, "Unreviewed samples: " + ", ".join(pending))
        tasks = [(decision, task, index) for index, decision in enumerate(entries)]
        samples = run_tasks(inspect_prepared, tasks, threads)
        manifest = dict(dataset_id=document["dataset_id"], samples=samples)
        choices = validate_decisions(manifest, document)
        check_destinations(output, samples)
        partial.mkdir(parents=True)
        tasks = [(meta, choices[meta["sample_id"]], task, partial) for meta in samples]
        results = run_tasks(filter_sample, tasks, threads)
        publish_results(partial, output, samples, existing_output)
        log_result(args, "apply", output=str(output), decisions=str(decisions_path),
                   decisions_sha256=sha(staged_decisions), samples=len(samples), threads=threads)
    except Exception:
        print(f"Filtering failed. Inspect partial results {partial} and scratch {task}", file=sys.stderr)
        raise
    shutil.rmtree(task)
    for result in results:
        print(f"{result['sample_id']}: retained {result['retained']} / {result['total']} sequences", flush=True)
    print(f"Results: {output}", flush=True)


def positive_integer(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("Use a positive integer.")
    return parsed


def main():
    parser = argparse.ArgumentParser(prog="eukasmchecker", description=__doc__)
    parser.add_argument("--version", action="version", version="EukAsmChecker " + VERSION)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare", help="Prepare one sample from FASTA, ECM and coverage tables.")
    for name in ("sample", "fasta", "ecm", "coverage", "seq-id", "seq-cov", "out"):
        p.add_argument("--" + name, required=True)
    r = sub.add_parser("review", help="Review prepared samples locally or create an offline ZIP.")
    r.add_argument("input")
    r.add_argument("--out", help="Write an offline ZIP instead of starting a local browser session.")
    a = sub.add_parser("apply", help="Apply the explicit decisions file to its referenced samples.")
    a.add_argument("decisions")
    a.add_argument("--out", required=True)
    for command in (r, a):
        command.add_argument("--threads", type=positive_integer, default=1,
                             help="Number of parallel sample workers (default: 1).")
    for command in (p, r, a):
        command.add_argument("--scratch", type=Path, default=Path(tempfile.gettempdir()),
                             help="Temporary work directory. HPC jobs must use their node-local /data directory.")
        command.add_argument("--log-dir", type=Path, help="Optional central run-log directory.")
    args = parser.parse_args()
    try:
        {"prepare": prepare, "review": review, "apply": apply}[args.command](args)
    except (ValueError, OSError, KeyError, TypeError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

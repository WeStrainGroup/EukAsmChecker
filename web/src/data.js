/* Browser input parsing. All sequence measurements retain their original precision. */
(function () {
  "use strict";
  const encoder = new TextEncoder();
  const decoder = new TextDecoder("utf-8", {fatal:true});
  const fail = message => {throw new Error(message);};
  const digest = text => sha256(text);
  async function* lines(file, hasher) {
    let pending = new Uint8Array(0), offset = 0;
    for await (const chunk of file.stream()) {
      if (hasher) hasher.update(chunk);
      const buffer = new Uint8Array(pending.length + chunk.length);
      buffer.set(pending); buffer.set(chunk, pending.length);
      let start = 0;
      for (let i=0;i<buffer.length;i++) {
        if (buffer[i]===10) {
          yield {bytes:buffer.subarray(start,i+1),start:offset+start,end:offset+i+1};
          start=i+1;
        }
      }
      offset+=start;
      pending=buffer.slice(start);
    }
    if (pending.length) yield {bytes:pending,start:offset,end:offset+pending.length};
  }
  function content(bytes) {
    let end=bytes.length;
    if (end && bytes[end-1]===10) end--;
    if (end && bytes[end-1]===13) end--;
    return bytes.subarray(0,end);
  }
  async function fasta(file) {
    const hasher=sha256.create(), records=[], ids=new Set();
    let current=null;
    for await (const line of lines(file,hasher)) {
      const value=content(line.bytes);
      if (value[0]===62) {
        if (current) {
          if (!current.length) fail("Empty FASTA sequence: "+current.id);
          current.end=line.start; records.push(current);
        }
        const header=decoder.decode(value.subarray(1));
        const id=header.trim().split(/\s+/)[0];
        if (!id || ids.has(id)) fail("Missing or duplicate FASTA ID: "+id);
        ids.add(id); current={id,length:0,start:line.start,end:file.size};
      } else {
        if (!current) fail("FASTA must begin with a > header.");
        for (const byte of value) if ([9,11,12,13,32].includes(byte)) fail("Whitespace in sequence: "+current.id);
        current.length+=value.length;
      }
    }
    if (!current || !current.length) fail("Empty FASTA or empty final sequence.");
    records.push(current);
    return {records,hash:hasher.hex()};
  }
  async function readTable(file, required, visit) {
    const hash=sha256.create();
    let header=null, count=0;
    for await (const line of lines(file,hash)) {
      const text=decoder.decode(content(line.bytes));
      if (!header) {
        header=text.replace(/^\uFEFF/,"").split("\t");
        if (new Set(header).size!==header.length || header.some(v=>!v)) fail("Invalid or duplicate column names in "+file.name);
        for (const column of required) if (!header.includes(column)) fail(file.name+": missing column "+column);
      } else if (text!=="") {
        const cells=text.split("\t");
        if (cells.length!==header.length) fail(file.name+": malformed TSV row "+(count+2));
        visit(Object.fromEntries(header.map((key,i)=>[key,cells[i]])));
        count++;
      }
    }
    if (!header || !count) fail("Empty table: "+file.name);
    return hash.hex();
  }
  function numeric(value, field) {
    const text=String(value).trim();
    if (!/^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$/.test(text)) fail("Invalid "+field+": "+value);
    const number=Number(text);
    if (!Number.isFinite(number) || number<0) fail("Invalid "+field+": "+value);
    return number;
  }
  async function prepare(sample, columns) {
    if (!/^[A-Za-z0-9][A-Za-z0-9_.-]*$/.test(sample.id) || sample.id==="." || sample.id==="..") fail("Use letters, numbers, dots, underscores or hyphens for sample IDs.");
    const parsed=await fasta(sample.fasta);
    const byId=new Map(parsed.records.map(record=>[record.id,record]));
    const coverage=new Map(), labels=new Map();
    const coverageHash=await readTable(sample.coverage,[columns.coverageId,columns.coverageValue],row=>{
      const id=row[columns.coverageId];
      if (!byId.has(id) || coverage.has(id)) fail("Extra or duplicate coverage ID: "+id);
      coverage.set(id,numeric(row[columns.coverageValue],"coverage"));
    });
    const required=[columns.ecmId,columns.ecmLabel];
    if (columns.ecmLength) required.push(columns.ecmLength);
    const ecmHash=await readTable(sample.ecm,required,row=>{
      const id=row[columns.ecmId], record=byId.get(id);
      if (!record || labels.has(id)) fail("Extra or duplicate ECM ID: "+id);
      if (columns.ecmLength && numeric(row[columns.ecmLength],"ECM length")!==record.length) fail("ECM/FASTA length differs: "+id);
      const value=row[columns.ecmLabel];
      if (!value || /[\s;]/.test(value)) fail("Invalid ECM label: "+id);
      labels.set(id,value==="euk"?"eukaryota":value);
    });
    if (coverage.size!==byId.size || labels.size!==byId.size) fail("FASTA, coverage and ECM must contain exactly the same IDs.");
    const rows=parsed.records.map(record=>[record.id,record.length,coverage.get(record.id),labels.get(record.id)]);
    const fingerprint=digest(JSON.stringify({coverage:coverageHash,ecm:ecmHash,fasta:parsed.hash}));
    const prepared={sample_id:sample.id,fingerprint,rows,ranges:parsed.records.map(record=>[record.start,record.end]),
      input_sha256:{fasta:parsed.hash,coverage:coverageHash,ecm:ecmHash},
      source_files:{fasta:sample.fasta.name,coverage:sample.coverage.name,ecm:sample.ecm.name}};
    const blob=new Blob([JSON.stringify(prepared)],{type:"application/json"});
    const eacHash=await blobHash(blob);
    const mapped=rows.map(row=>({id:row[0],length:row[1],coverage:row[2],label:row[3]}));
    const reference=EACCore.reference(mapped), groups={};
    for (const row of rows) {
      const group=groups[row[3]] ||= {contigs:0,bp:0};
      group.contigs++;group.bp+=row[1];
    }
    const meta={sample_id:sample.id,fingerprint,eac_path:"browser://"+fingerprint+"/"+sample.id+".eac.json",
      eac_sha256:eacHash,contigs:rows.length,bp:rows.reduce((sum,row)=>sum+row[1],0),groups,
      reference:{n:reference.n,cov_10kb:reference.median,min_coverage:reference.cutoff}};
    return {meta,blob};
  }
  async function blobHash(blob) {
    const hash=sha256.create();
    for await (const chunk of blob.stream()) hash.update(chunk);
    return hash.hex();
  }
  function validateDecision(decision,meta) {
    if (!decision || decision.sample_id!==meta.sample_id || decision.fingerprint!==meta.fingerprint ||
        decision.eac_path!==meta.eac_path || decision.eac_sha256!==meta.eac_sha256) fail("Decisions do not match "+meta.sample_id);
    if (decision.reviewed!==true || typeof decision.retain_sample!=="boolean") fail("Confirm review for "+meta.sample_id);
    if (!Array.isArray(decision.labels) || new Set(decision.labels).size!==decision.labels.length ||
        decision.labels.some(label=>typeof label!=="string")) fail("Invalid ECM labels.");
    if ("excluded_seq_ids" in decision && (!Array.isArray(decision.excluded_seq_ids) ||
        decision.excluded_seq_ids.some(id=>typeof id!=="string"||!id) ||
        new Set(decision.excluded_seq_ids).size!==decision.excluded_seq_ids.length)) fail("Invalid manual exclusions.");
    if (decision.retain_sample) {
      if (!Number.isSafeInteger(decision.min_length_bp) || decision.min_length_bp<0 ||
          !Number.isSafeInteger(decision.min_coverage) || decision.min_coverage<0 ||
          !decision.labels.length) fail("Invalid thresholds for "+meta.sample_id);
      if (decision.max_coverage!==null && (!Number.isSafeInteger(decision.max_coverage) || decision.max_coverage<decision.min_coverage)) fail("Invalid maximum coverage.");
    }
    if (!Number.isSafeInteger(decision.expected_kept_contigs) || decision.expected_kept_contigs<0 ||
        !Number.isSafeInteger(decision.expected_kept_bp) || decision.expected_kept_bp<0) fail("Invalid reviewed counts.");
  }
  async function apply(preparedBlob,original,decision,meta) {
    validateDecision(decision,meta);
    if (await blobHash(preparedBlob)!==meta.eac_sha256) fail("Prepared input changed.");
    const prepared=JSON.parse(await preparedBlob.text());
    if (await blobHash(original)!==prepared.input_sha256.fasta) fail("Original FASTA differs from the reviewed input.");
    const settings={retain_sample:decision.retain_sample,length:decision.min_length_bp,coverage:decision.min_coverage,
      maxCoverage:decision.max_coverage,labels:decision.labels,
      excludedIds:EACCore.exclusions(decision.excluded_seq_ids ?? [],prepared.rows.map(row=>row[0]))};
    const parts=[], table=["seq_id\tlength_bp\tcoverage\tecm_label\tfilter_status\tfilter_reason\n"];
    let kept=0,bp=0,rangeStart=null,rangeEnd=null;
    function flush() {if(rangeStart!==null) parts.push(original.slice(rangeStart,rangeEnd));rangeStart=null;}
    prepared.rows.forEach((row,index)=>{
      const why=EACCore.reasons({id:row[0],length:row[1],coverage:row[2],label:row[3]},settings);
      const retain=why.length===0;
      table.push([...row,retain?"retained":"removed",why.length?why.join(";"):"pass"].join("\t")+"\n");
      if (retain) {
        kept++;bp+=row[1];
        const [start,end]=prepared.ranges[index];
        if(rangeStart===null)rangeStart=start;
        rangeEnd=end;
      } else flush();
    });
    flush();
    if(kept!==decision.expected_kept_contigs || bp!==decision.expected_kept_bp) fail("Retained counts differ from review: "+meta.sample_id);
    return {fasta:new Blob(parts,{type:"text/plain"}),table:new Blob(table,{type:"text/tab-separated-values"}),kept,bp};
  }
  globalThis.EACData={lines,fasta,prepare,apply,blobHash,numeric,validateDecision};
})();

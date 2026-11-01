"use strict";
importScripts("vendor/sha256.min.js","core.js","data.js","zip.js");
let samples=new Map(),manifest=null;
async function preparedData(sample){
  if(sample.blob)return JSON.parse(await sample.blob.text());
  if(!sample.payload){
    const response=await fetch(sample.demoPath);
    if(!response.ok)throw new Error("Could not load the example profile.");
    const blob=await response.blob();
    if(await EACData.blobHash(blob)!==sample.meta.eac_sha256)throw new Error("Example checksum mismatch.");
    sample.payload=JSON.parse(await blob.text());
  }
  return sample.payload;
}
self.addEventListener("message",async event=>{
  const message=event.data;
  const progress=(stage,done,total,sample="")=>self.postMessage({id:message.id,type:"progress",stage,done,total,sample});
  try {
    if(message.type==="demo"){
      const response=await fetch("demo/manifest.json");
      if(!response.ok)throw new Error("Examples are unavailable.");
      manifest=await response.json();
      const legacyIndex=["WeGIF_334","WeGIF_339","WeGIF_653"].indexOf(message.sample);
      if(legacyIndex>=0)message.sample=manifest.samples[legacyIndex].sample_id;
      manifest.initial_sample=manifest.samples.some(row=>row.sample_id===message.sample)?message.sample:manifest.samples[0].sample_id;
      samples=new Map(manifest.samples.map(meta=>[meta.sample_id,{meta,demoPath:meta.data_file}]));
      self.postMessage({id:message.id,type:"result",value:manifest});
      return;
    }
    if(message.type==="prepare"){
      samples=new Map();manifest=null;
      const ids=new Set();
      for(const sample of message.samples){
        if(ids.has(sample.id))throw new Error("Duplicate sample ID: "+sample.id);
        ids.add(sample.id);
      }
      if(!ids.size)throw new Error("Add at least one sample.");
      let index=0;
      for(const sample of message.samples){
        progress("prepare",index,message.samples.length,sample.id);
        const result=await EACData.prepare(sample,message.columns);
        result.meta.index=index;result.meta.data_file="browser:"+sample.id;
        samples.set(sample.id,{...result,fasta:sample.fasta});
        index++;
      }
      const meta=[...samples.values()].map(sample=>sample.meta).sort((a,b)=>a.sample_id.localeCompare(b.sample_id,"en"));
      const identity=sha256(JSON.stringify(meta.map(sample=>[sample.sample_id,sample.fingerprint,sample.eac_sha256,sample.eac_path])));
      manifest={schema_version:3,tool_version:"1.0",browser_local:true,dataset_id:identity,
        project_name:message.projectName||"Assembly review",columns:["seq_id","length_bp","coverage","ecm_label"],samples:meta};
      self.postMessage({id:message.id,type:"result",value:manifest});
    }else if(message.type==="payload"){
      const sample=samples.get(message.sample);
      if(!sample)throw new Error("Unknown sample.");
      const prepared=await preparedData(sample);
      self.postMessage({id:message.id,type:"result",value:{sample_id:prepared.sample_id,fingerprint:prepared.fingerprint,rows:prepared.rows}});
    }else if(message.type==="export"){
      const document=message.decisions;
      if(!manifest || document.schema_version!==7 || document.dataset_id!==manifest.dataset_id ||
         document.samples.length!==samples.size || new Set(document.samples.map(row=>row.sample_id)).size!==samples.size)
        throw new Error("Decisions do not match this project.");
      const entries=[],summary=["sample_id\tretain_sample\tmin_length_bp\tmin_coverage\tmax_coverage\tlabels\tretained_sequences\tretained_bp\tnote\n"];
      const sources={};
      let index=0;
      for(const decision of document.samples){
        const sample=samples.get(decision.sample_id);
        if(!sample)throw new Error("Unknown sample in decisions.");
        progress("filter",index,samples.size,decision.sample_id);
        const name=decision.sample_id;
        let output;
        if(manifest.browser_demo){
          EACData.validateDecision(decision,sample.meta);
          const prepared=await preparedData(sample);
          const rows=prepared.rows.map(row=>({id:row[0],length:row[1],coverage:row[2],label:row[3]}));
          const state={retain_sample:decision.retain_sample,length:decision.min_length_bp,coverage:decision.min_coverage,maxCoverage:decision.max_coverage,labels:decision.labels};
          const result=EACCore.evaluate(rows,state);
          if(result.kept!==decision.expected_kept_contigs||result.keptBp!==decision.expected_kept_bp)throw new Error("Reviewed example counts differ.");
          const table=["seq_id\tlength_bp\tcoverage\tecm_label\tfilter_status\tfilter_reason\n"];
          prepared.rows.forEach((row,index)=>{const mark=result.annotated[index];table.push([...row,mark.keep?"retained":"removed",mark.reasons.join(";")||"pass"].join("\t")+"\n");});
          output={table:new Blob(table),kept:result.kept,bp:result.keptBp};
        }else{
          output=await EACData.apply(sample.blob,sample.fasta,decision,sample.meta);
          entries.push({name:name+"/"+name+".fasta",blob:output.fasta});
        }
        entries.push({name:name+"/"+name+".sequences.tsv",blob:output.table});
        summary.push([name,decision.retain_sample,decision.min_length_bp,decision.min_coverage,decision.max_coverage,
          decision.labels.join(";"),output.kept,output.bp,String(decision.note||"").replace(/[\t\r\n]+/g," ")].join("\t")+"\n");
        const prepared=await preparedData(sample);
        sources[name]=manifest.browser_demo?{example:sample.demoPath,fingerprint:sample.meta.fingerprint}:{files:prepared.source_files,sha256:prepared.input_sha256};
        index++;
      }
      const saved={...document,origin:manifest.browser_demo?"EukAsmChecker example":"EukAsmChecker browser",source_files:sources};
      entries.push({name:"decisions.json",blob:new Blob([JSON.stringify(saved,null,2)+"\n"])},
        {name:"summary.tsv",blob:new Blob(summary)});
      const zip=await EACZip.build(entries,(done,total)=>progress("archive",done,total));
      self.postMessage({id:message.id,type:"result",value:zip});
    }else throw new Error("Unknown operation.");
  }catch(error){
    if(message.type==="prepare"){samples.clear();manifest=null;}
    self.postMessage({id:message.id,type:"error",message:error.message||String(error)});
  }
});

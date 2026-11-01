"use strict";
const $=id=>document.getElementById(id);
for(const input of document.querySelectorAll("input"))input.value=input.type==="file"?"":input.defaultValue;
try{
  for(const key of Object.keys(localStorage))if(key.startsWith("genome-review-"))localStorage.removeItem(key);
}catch(e){}
window.addEventListener("pageshow",event=>{if(event.persisted)location.reload();});
const en=Object.fromEntries([...document.querySelectorAll("[data-text]")].map(element=>[element.dataset.text,element.textContent]));
const zh={
  description:"根据序列长度、reads coverage 和 EukContigMiner 标签，对真核基因组组装进行交互式审阅与过滤。", inputsTitle:"输入文件",
  demo:"示例", fastaButton:"FASTA 文件", coverageButton:"Coverage 表格", ecmButton:"EukContigMiner 表格",
  matching:"请检查下方样本名称和文件对应关系。", empty:"添加文件后，样本会显示在这里。", sample:"样本",
  columnsTitle:"对应列", columnsHelp:"选择对应列名，序列 ID 须与 FASTA 标题中的第一个字段一致。",
  coverageLegend:"Coverage 表格", sequenceId:"序列 ID", readsCoverage:"Reads coverage", ecmLegend:"EukContigMiner 表格", ecmLabel:"标签", ecmLength:"长度核对",
  projectName:"项目名称", prepare:"开始审阅 →", savedDecisions:"导入审阅记录 (.json)", cancel:"取消"
};
let language="en",files={fasta:[],coverage:[],ecm:[]},samples=[],worker=null,manifest=null,sequence=0,busy=false;
const pending=new Map();
let archiveURL=null,frameReady=false,prepareGeneration=0;
const say=(english,chinese)=>language==="zh"?chinese:english;
function setLanguage(){
  document.documentElement.lang=language==="en"?"en":"zh-CN";
  for(const element of document.querySelectorAll("[data-text]"))element.textContent=(language==="zh"?zh:en)[element.dataset.text]||en[element.dataset.text];
  $("language").textContent=language==="en"?"中文":"English";
  updateDecisionCaption();
  if(frameReady)$("reviewFrame").contentWindow.postMessage({type:"eac:language",language},location.origin);
}
function updateDecisionCaption(){
  const file=$("decisionFile").files[0];
  $("decisionCaption").textContent=file?file.name:say(en.savedDecisions,zh.savedDecisions);
  $("decisionCaption").parentElement.title=file?file.name:say("Select decisions.json from a ZIP downloaded with Save results. It restores thresholds, labels, sample choices and notes. Select the matching FASTA and tables again.","选择之前通过“保存结果”下载的 ZIP 中的 decisions.json，恢复阈值、标签、样本取舍和备注。还需重新选择对应的 FASTA 和表格。");
  $("clearDecisions").hidden=!file;
  $("clearDecisions").setAttribute("aria-label",say("Clear imported review","清除导入的审阅记录"));
}
async function readSavedDecisions(){
  const file=$("decisionFile").files[0];
  return file?JSON.parse(await file.text()):null;
}
$("decisionFile").onchange=updateDecisionCaption;
$("clearDecisions").onclick=()=>{$("decisionFile").value="";updateDecisionCaption();};
$("language").onclick=()=>{language=language==="en"?"zh":"en";setLanguage();};
function showError(message){$("error").textContent=message;$("error").hidden=!message;}
function path(file){return file.name;}
function sampleName(file){
  const name=file.name;
  let stem=name.replace(/\.(fasta|fna|fas|fa|tsv|txt)$/i,"").replace(/[._-](coverage|classification|ecm)$/i,"");
  return stem.replace(/[^A-Za-z0-9_.-]+/g,"_")||"sample";
}
function guessed(file,kind){
  const match=files[kind].findIndex(candidate=>sampleName(candidate)===sampleName(file));
  return match>=0?match:files.fasta.length===1&&files[kind].length===1?0:-1;
}
function selectFiles(kind,list){
  files[kind]=[...list];
  if(kind==="fasta")samples=files.fasta.map((file,index)=>({id:sampleName(file),fasta:index,coverage:guessed(file,"coverage"),ecm:guessed(file,"ecm")}));
  else samples.forEach(sample=>{sample[kind]=guessed(files.fasta[sample.fasta],kind);});
  renderSamples();refreshColumns();
}
for(const kind of ["fasta","coverage","ecm"])$(kind+"Files").onchange=event=>selectFiles(kind,event.target.files);
function renderSamples(){
  for(const kind of ["fasta","coverage","ecm"])$(kind+"Count").textContent=files[kind].length;
  $("empty").hidden=!!samples.length;$("sampleWrap").hidden=!samples.length;$("sampleRows").replaceChildren();
  samples.forEach(sample=>{
    const row=document.createElement("tr"),idCell=document.createElement("td"),input=document.createElement("input");
    input.value=sample.id;input.setAttribute("aria-label","Sample ID");input.oninput=()=>sample.id=input.value;
    idCell.append(input);row.append(idCell);
    const fastaCell=document.createElement("td");fastaCell.textContent=path(files.fasta[sample.fasta]);fastaCell.title=fastaCell.textContent;row.append(fastaCell);
    for(const kind of ["coverage","ecm"]){
      const cell=document.createElement("td"),select=document.createElement("select");select.setAttribute("aria-label",kind+" file for "+sample.id);
      select.add(new Option(say("Choose a table","选择表格"),"-1"));
      files[kind].forEach((file,index)=>select.add(new Option(path(file),String(index))));
      select.value=String(sample[kind]);select.onchange=()=>{sample[kind]=Number(select.value);refreshColumns();};
      cell.append(select);row.append(cell);
    }
    $("sampleRows").append(row);
  });
}
async function headers(file){
  if(!file)return[];
  const text=await file.slice(0,65536).text();
  return text.replace(/^\uFEFF/,"").split(/\r?\n/,1)[0].split("\t");
}
async function refreshColumns(){
  const generation=++prepareGeneration,selected=samples[0];
  const coverage=await headers(files.coverage[selected?.coverage]||files.coverage[0]);
  const ecm=await headers(files.ecm[selected?.ecm]||files.ecm[0]);
  if(generation!==prepareGeneration)return;
  const definitions={coverageId:[coverage,["scaffold_id","seq_id","sequence","contig_id","identifier","id"]],coverageValue:[coverage,["mean_depth","coverage","reads_cov","mean_coverage"]],
    ecmId:[ecm,["contig_id","seq_id","scaffold_id","id"]],ecmLabel:[ecm,["label_euk","ecm_label","label"]],ecmLength:[ecm,["length_bp","length"]]};
  for(const [id,[columns,preferred]] of Object.entries(definitions)){
    const element=$(id),previous=element.value;
    element.replaceChildren();
    if(id==="ecmLength")element.add(new Option(say("Use FASTA length","使用 FASTA 长度"),""));
    else element.add(new Option(say("Choose a column","选择列名"),""));
    for(const name of columns)element.add(new Option(name,name));
    element.value=columns.includes(previous)?previous:preferred.find(name=>columns.includes(name))|| "";
  }
}
function startWorker(){
  if(worker)worker.terminate();
  pending.clear();worker=new Worker("worker.js");
  worker.onmessage=event=>{
    const message=event.data,task=pending.get(message.id);
    if(message.type==="progress"){
      const labels={prepare:say("Preparing","准备"),filter:say("Filtering","筛选"),archive:say("Packaging","打包")};
      $("activityText").textContent=labels[message.stage]+" "+message.done+" / "+message.total+" "+message.sample;
      $("progress").value=message.total?message.done/message.total:0;
    }else if(task){pending.delete(message.id);message.type==="error"?task.reject(new Error(message.message)):task.resolve(message.value);}
  };
  worker.onerror=event=>{for(const task of pending.values())task.reject(new Error(event.message||"Worker failed."));pending.clear();};
}
function request(type,data={}){
  return new Promise((resolve,reject)=>{const id=++sequence;pending.set(id,{resolve,reject});worker.postMessage({id,type,...data});});
}
function setBusy(value){
  busy=value;$("activity").hidden=!value;
  $("demo").disabled=value;
  $("reviewFrame").inert=value;
  for(const element of $("setup").querySelectorAll("input,select,button"))element.disabled=value;
  if(value){showError("");$("downloadStatus").hidden=true;$("activityText").textContent=say("Working…","正在处理…");$("progress").value=0;}
}
$("cancel").onclick=()=>{
  if(worker)worker.terminate();worker=null;
  for(const task of pending.values())task.reject(new Error(say("Cancelled. Prepare the inputs again to continue.","已取消。重新生成审阅页面后可继续。")));
  pending.clear();manifest=null;frameReady=false;
  document.body.classList.remove("reviewing");$("reviewArea").hidden=true;$("setup").hidden=false;setBusy(false);
};
async function openReview(project,decisions=null){
  const response=await fetch("review-template.html");
  if(!response.ok)throw new Error("Could not load the review interface.");
  const template=await response.text();
  manifest=project;
  const html=template.replace("/*PROJECT_MANIFEST*/null",JSON.stringify({...manifest,initial_decisions:decisions}).replace(/</g,"\\u003c"));
  frameReady=false;
  $("reviewFrame").srcdoc=html;
  document.body.classList.add("reviewing");$("setup").hidden=true;$("reviewArea").hidden=false;window.scrollTo(0,0);
}
async function openExamples(sample="sample_1"){
  if(busy)return;
  try{
    setBusy(true);
    startWorker();
    await openReview(await request("demo",{sample}));
  }
  catch(error){showError(error.message);}
  finally{setBusy(false);}
}
$("prepare").onclick=async()=>{
  if(busy)return;
  try{
    if(!samples.length)throw new Error(say("Choose FASTA files first.","请先选择 FASTA 文件。"));
    const input=samples.map(sample=>{
      if(sample.coverage<0||sample.ecm<0)throw new Error(say("Assign both tables for ","请为样本指定两个表格 ")+sample.id);
      return{id:sample.id.trim(),fasta:files.fasta[sample.fasta],coverage:files.coverage[sample.coverage],ecm:files.ecm[sample.ecm]};
    });
    setBusy(true);
    const decisions=await readSavedDecisions();
    startWorker();
    const columns=Object.fromEntries(["coverageId","coverageValue","ecmId","ecmLabel","ecmLength"].map(id=>[id,$(id).value]));
    if(["coverageId","coverageValue","ecmId","ecmLabel"].some(id=>!columns[id]))throw new Error(say("Choose all required column names.","请选择所有必需的列名。"));
    const project=await request("prepare",{samples:input,columns,projectName:$("projectName").value.trim()});
    await openReview(project,decisions);
  }catch(error){showError(error.message);}
  finally{setBusy(false);}
};
window.addEventListener("message",async event=>{
  if(event.origin!==location.origin||event.source!==$("reviewFrame").contentWindow)return;
  const message=event.data;
  if(message.type==="eac:ready"){
    frameReady=true;
    $("reviewFrame").contentWindow.postMessage({type:"eac:language",language},location.origin);
  }else if(message.type==="eac:back"||message.type==="eac:restore-error"){
    document.body.classList.remove("reviewing");
    $("reviewArea").hidden=true;$("setup").hidden=false;window.scrollTo(0,0);
    if(message.type==="eac:restore-error"){frameReady=false;showError(message.message);}
  }else if(message.type==="eac:sample"){
    try{const payload=await request("payload",{sample:message.sample});event.source.postMessage({type:"eac:payload",payload},location.origin);}
    catch(error){showError(error.message);}
  }else if(message.type==="eac:export"){
    try{
      if(busy||!frameReady)return;
      setBusy(true);
      const archive=await request("export",{decisions:message.document});
      if(archiveURL)URL.revokeObjectURL(archiveURL);
      archiveURL=URL.createObjectURL(archive);
      const name=manifest.browser_demo?"example.zip":(manifest.project_name.replace(/[^A-Za-z0-9_.-]+/g,"_")||"assembly")+".filtered.zip";
      const link=document.createElement("a");link.href=archiveURL;link.download=name;link.textContent=say("Download ZIP","下载 ZIP");
      $("downloadStatus").replaceChildren(document.createTextNode(say("Results are ready. ","结果已生成。")),link);
      $("downloadStatus").hidden=false;link.click();
    }catch(error){showError(error.message);}
    finally{setBusy(false);}
  }
});
$("demo").onclick=()=>openExamples();
refreshColumns();
const requestedExample=new URLSearchParams(location.search).get("demo");
if(requestedExample!==null){
  const url=new URL(location.href);
  url.searchParams.delete("demo");
  history.replaceState(null,"",url);
  openExamples(requestedExample);
}
setLanguage();

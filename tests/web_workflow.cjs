const assert=require("node:assert/strict");
const fs=require("node:fs");
const path=require("node:path");
const crypto=require("node:crypto");
const vm = require("node:vm");
const root=path.resolve(__dirname,".."),out=process.argv[2];
global.sha256=require(path.join(root,"web/dist/vendor/sha256.min.js"));
require(path.join(root,"web/dist/core.js"));
require(path.join(root,"web/dist/data.js"));
require(path.join(root,"web/dist/zip.js"));
// Files are assigned by the chosen input type, with samples matched by name.
const appSource = fs.readFileSync(path.join(root, "web/src/app.js"), "utf8");
const controls = Object.fromEntries(["fastaFiles", "coverageFiles", "ecmFiles"].map(id => [id, {}]));
const inputContext = {
  files: {fasta: [], coverage: [], ecm: []}, samples: [],
  $: id => controls[id], renderSamples() {}, refreshColumns() {}
};
vm.createContext(inputContext);
vm.runInContext(appSource.slice(appSource.indexOf("function path(file)"), appSource.indexOf("function renderSamples()")), inputContext);
controls.fastaFiles.onchange({target: {files: [{name: "a.fasta"}, {name: "b.fasta"}]}});
controls.coverageFiles.onchange({target: {files: [{name: "b.coverage.tsv"}, {name: "a.coverage.tsv"}]}});
controls.ecmFiles.onchange({target: {files: [{name: "a.classification.tsv"}, {name: "b.classification.tsv"}]}});
assert.equal(inputContext.samples.length, 2);
for (const sample of inputContext.samples) {
  assert.equal(inputContext.files.coverage[sample.coverage].name, sample.id + ".coverage.tsv");
  assert.equal(inputContext.files.ecm[sample.ecm].name, sample.id + ".classification.tsv");
}

// Mixed-label caches must preserve filtering, seeded order and threshold-drag reuse.
const reviewSource = fs.readFileSync(path.join(root, "src/eukasmchecker/review.html"), "utf8");
// Saving/restoring must retain manual IDs and accept existing schema 7 reviews.
const advancedRule={id:1,enabled:true,labels:["other"],length:{operator:">",value:500},coverage:{operator:"any",value:0}};
const ruleLabels=["eukaryota","other","NA"];
assert.equal(EACCore.ruleMatches({label:"other",length:500,coverage:20},advancedRule),false);
assert.equal(EACCore.ruleMatches({label:"other",length:501,coverage:20},advancedRule),true);
for(const operator of ["any",">",">=","<","<="]){
 const expected={any:true,">":false,">=":true,"<":false,"<=":true}[operator];
 assert.equal(EACCore.ruleMatches({label:"other",length:500,coverage:20},{...advancedRule,length:{operator,value:500}}),expected);
}
for(const invalid of [null,"rules",[null],[{...advancedRule,id:true}],[advancedRule,advancedRule],
 [{...advancedRule,enabled:1}],[{...advancedRule,labels:[]}],[{...advancedRule,labels:["missing"]}],
 [{...advancedRule,coverage:{operator:">",value:10.5}}],[{...advancedRule,length:{operator:"=",value:500}}]]){
 assert.throws(()=>EACCore.rules(invalid,ruleLabels),/advanced/);
}
const savedMeta={sample_id:"sample_A",fingerprint:"fp",eac_path:"browser:a",eac_sha256:"hash"};
const savedState={retain_sample:true,length:500,coverage:10,maxCoverage:null,labels:["eukaryota"],
  excludedIds:new Set(["a"]),advancedRules:[advancedRule],mode:"manual",reviewed:false,reviewed_at:null,note:"",expected_kept_contigs:null,expected_kept_bp:null};
const savedContext={C:EACCore,PROJECT:{dataset_id:"dataset",samples:[savedMeta]},settings:new Map([["sample_A",savedState]]),
  metadata:new Map([["sample_A",savedMeta]]),availableLabels:["eukaryota","other","NA"],
  normalizeLabels:labels=>labels,t:key=>key,valid:()=>true};
vm.createContext(savedContext);
vm.runInContext(reviewSource.slice(reviewSource.indexOf("function documentData()"),reviewSource.indexOf("let serverReady=")),savedContext);
vm.runInContext(reviewSource.slice(reviewSource.indexOf("function restoreDocument("),reviewSource.indexOf("function listSamples(")),savedContext);
const savedReview=savedContext.documentData();
assert.equal(savedReview.schema_version,9);
assert.equal(savedReview.samples[0].advanced_rules[0].length.value,500);
assert.deepEqual(Array.from(savedReview.samples[0].excluded_seq_ids),["a"]);
savedContext.restoreDocument(savedReview);
assert.deepEqual(Array.from(savedContext.settings.get("sample_A").excludedIds),["a"]);
const legacyReview=JSON.parse(JSON.stringify(savedReview));legacyReview.schema_version=7;delete legacyReview.samples[0].excluded_seq_ids;delete legacyReview.samples[0].advanced_rules;
savedContext.restoreDocument(legacyReview);
assert.equal(savedContext.settings.get("sample_A").excludedIds.size,0);
assert.equal(savedContext.settings.get("sample_A").advancedRules.length,0);
const invalidReview=JSON.parse(JSON.stringify(savedReview));invalidReview.samples[0].excluded_seq_ids=["a","a"];
assert.throws(()=>savedContext.restoreDocument(invalidReview),/wrongFields/);

// Click and repeated rectangles toggle only their target points.
const selectionContext={selectedIds:new Set(),updateSelectionControls(){},draw(){}};
vm.createContext(selectionContext);
vm.runInContext(reviewSource.slice(reviewSource.indexOf("function togglePointSelection("),reviewSource.indexOf("function removeSelectedPoints(")),selectionContext);
selectionContext.togglePointSelection(["a"]);
assert.deepEqual(Array.from(selectionContext.selectedIds),["a"]);
selectionContext.togglePointSelection(["a"]);
assert.equal(selectionContext.selectedIds.size,0);
selectionContext.togglePointSelection(["a","b"]);
assert.equal(selectionContext.selectedIds.size,2);
selectionContext.togglePointSelection(["a","b"]);
assert.equal(selectionContext.selectedIds.size,0);
selectionContext.togglePointSelection(["a"]);
selectionContext.togglePointSelection(["a","b"]);
assert.deepEqual(Array.from(selectionContext.selectedIds),["a","b"]);
selectionContext.togglePointSelection(["c"]);
selectionContext.togglePointSelection(["a","b"]);
assert.deepEqual(Array.from(selectionContext.selectedIds),["c"]);
selectionContext.togglePointSelection([]);
assert.deepEqual(Array.from(selectionContext.selectedIds),["c"]);

// Remove/restore acts on the current selection and preserves outside exclusions.
const actionState={excludedIds:new Set(["outside","a"]),reviewed:true};
const actionControls=Object.fromEntries(["selectionCount","selectionHint","removeSelected","restoreSelected"].map(id=>[id,{addEventListener(){}}]));
const actionContext={selectedIds:new Set(["a","b"]),loaded:true,lang:"en",format:String,
  translations:{en:{selectionCounts:(n,m)=>`${n} selected · ${m} removed`}},
  t:key=>key,$:id=>actionControls[id],state:()=>actionState,
  invalidate(){actionState.reviewed=false;},refresh(){actionContext.updateSelectionControls();}};
vm.createContext(actionContext);
vm.runInContext(reviewSource.slice(reviewSource.indexOf("function updateSelectionControls("),reviewSource.indexOf("function cancelPointSelection(")),actionContext);
vm.runInContext(reviewSource.slice(reviewSource.indexOf("function removeSelectedPoints("),reviewSource.indexOf("function drawSelectedPoints(")),actionContext);
actionContext.updateSelectionControls();
assert.equal(actionControls.selectionCount.textContent,"2 selected · 2 removed");
assert.equal(actionControls.removeSelected.disabled,false);
assert.equal(actionControls.restoreSelected.disabled,false);
actionContext.removeSelectedPoints();
assert.equal(actionContext.selectedIds.size,0);
assert.deepEqual(Array.from(actionState.excludedIds),["outside","a","b"]);
assert.equal(actionControls.selectionCount.textContent,"0 selected · 3 removed");
assert.equal(actionControls.removeSelected.disabled,true);
assert.equal(actionState.reviewed,false);
actionContext.selectedIds=new Set(["a"]);
actionContext.restoreSelectedPoints();
assert.equal(actionContext.selectedIds.size,0);
assert.deepEqual(Array.from(actionState.excludedIds),["outside","b"]);
assert.equal(actionControls.selectionCount.textContent,"0 selected · 2 removed");
assert.equal(actionControls.restoreSelected.disabled,true);
actionState.reviewed=true;
actionContext.restoreSelectedPoints();
assert.equal(actionState.reviewed,true);
actionContext.selectedIds=new Set();
actionContext.updateSelectionControls();
assert.equal(actionControls.selectionCount.textContent,"0 selected · 2 removed");
assert.equal(actionControls.removeSelected.disabled,true);
assert.equal(actionControls.restoreSelected.disabled,true);

let rasterizedPoints = 0;
const normalPalette = ["#1d1d1f", "#709487", "#96969e"];
const fadedPalette = ["#bbbbbc", "#d4dfdb", "#e0e0e2"];
const renderRows = [
  {id: "c", length: 499, coverage: 10, label: "NA"},
  {id: "b", length: 500, coverage: 10, label: "eukaryota"},
  {id: "a", length: 1000, coverage: 20, label: "other"},
  {id: "d", length: 2000, coverage: 30, label: "eukaryota"}
];
function mockCanvas() {
  const canvas = {points: []};
  const context = {
    scale() {}, beginPath() {},
    clearRect() {canvas.points = [];},
    arc(x, y) {this.point = {x, y}; rasterizedPoints++;},
    fill() {canvas.points.push({...this.point, color: this.fillStyle, alpha: this.globalAlpha});}
  };
  canvas.getContext = () => context;
  return canvas;
}
const renderContext = {
  C:EACCore,pointLayers: null, hoverPoints: [], rows: renderRows,
  availableLabels: ["eukaryota", "other", "NA"], labelColor: index => normalPalette[index],
  document: {createElement: mockCanvas}
};
vm.createContext(renderContext);
vm.runInContext(reviewSource.slice(reviewSource.indexOf("  function preparePointLayers("), reviewSource.indexOf("  function drawDistributions(")), renderContext);
const geometry = {w: 300, h: 200, left: 20, right: 280, top: 20, bottom: 180,
  x: {a: 2, b: 4, f: length => 20 + (Math.log10(length) - 2) / 2 * 260},
  y: {b: 5, f: coverage => 180 - Math.asinh(coverage) / 5 * 160}};
for (const retain of [true, false]) for (const labels of [["eukaryota"], ["other", "NA"], ["eukaryota", "other", "NA"]]) for (const excluded of [[], ["a", "b"]]) for (const advancedRules of [[],[advancedRule]]) {
  const settings = {length: 500, coverage: 10, maxCoverage: 20, retain_sample: retain, labels, excludedIds: new Set(excluded),advancedRules};
  renderContext.result = EACCore.evaluate(renderRows, settings);
  renderContext.preparePointLayers(geometry, 1, settings);
  const cached = renderContext.pointLayers;
  const rendered = rasterizedPoints;
  const shown = renderRows.map(() => []);
  let rectangles = [];
  const context = {
    save() {}, restore() {}, beginPath() {rectangles = [];},
    rect(x, y, w, h) {rectangles.push({x, y, w, h});}, clip() {},
    drawImage(canvas) {
      for (const point of canvas.points) {
        const index = renderContext.hoverPoints.findIndex(item => item.x === point.x && item.y === point.y);
        const inside = rectangles.filter(r => point.x >= r.x && point.x <= r.x + r.w && point.y >= r.y && point.y <= r.y + r.h).length % 2;
        if (inside) shown[index].push([point.color, point.alpha]);
      }
    }
  };
  renderContext.drawPoints(context, geometry, settings, geometry.x.f(500), geometry.y.f(10), geometry.y.f(20));
  shown.forEach((appearances, index) => {
    const label = renderContext.availableLabels.indexOf(renderRows[index].label);
    const palette = renderContext.result.annotated[index].keep ? normalPalette : fadedPalette;
    assert.deepEqual(appearances, [[palette[label], .55]]);
  });
  const changedThresholds = {...settings, length: 800, coverage: 15, maxCoverage: 25};
  renderContext.result = EACCore.evaluate(renderRows, changedThresholds);
  renderContext.preparePointLayers(geometry, 1, changedThresholds);
  assert.equal(renderContext.pointLayers, cached);
  assert.equal(rasterizedPoints, rendered);
}
// Adding or disabling rules rebuilds selected appearances while preserving draw order.
const beforeRules={length:500,coverage:10,maxCoverage:null,retain_sample:true,labels:ruleLabels,excludedIds:new Set(),advancedRules:[]};
renderContext.result=EACCore.evaluate(renderRows,beforeRules);
renderContext.preparePointLayers(geometry,1,beforeRules);
const ruleCache=renderContext.pointLayers;
const withRules={...beforeRules,advancedRules:[advancedRule]};
renderContext.result=EACCore.evaluate(renderRows,withRules);
renderContext.preparePointLayers(geometry,1,withRules);
assert.notEqual(renderContext.pointLayers,ruleCache);
assert.deepEqual(Array.from(renderContext.pointLayers.order),Array.from(ruleCache.order));
const afterRulePaint=rasterizedPoints;
renderContext.preparePointLayers(geometry,1,{...withRules,coverage:11});
assert.equal(rasterizedPoints,afterRulePaint);
// Fully overlapping, initially label-sorted input must be interleaved on both layers.
const overlapRows = ["eukaryota", "other", "NA"].flatMap(label =>
  Array.from({length: 30}, () => ({length: 1000, coverage: 20, label})));
renderContext.rows = overlapRows;
renderContext.pointLayers = null;
const allLabels = {length: 500, coverage: 10, maxCoverage: null, retain_sample: true, labels: ["eukaryota", "other", "NA"]};
renderContext.result = EACCore.evaluate(overlapRows, allLabels);
renderContext.preparePointLayers(geometry, 1, allLabels);
const firstOrder = Array.from(renderContext.pointLayers.order);
assert.equal(new Set(firstOrder).size, overlapRows.length);
const colors = renderContext.pointLayers.layers[1].canvas.points.map(point => point.color);
assert.deepEqual(renderContext.pointLayers.layers[0].canvas.points.map(point => normalPalette[fadedPalette.indexOf(point.color)]), colors);
assert(colors.filter((color, index) => index > 0 && color !== colors[index - 1]).length > 30);
for (const color of normalPalette) assert.equal(colors.filter(value => value === color).length, 30);
const fadedCanvas = renderContext.pointLayers.layers[0].canvas;
const beforeSelection = rasterizedPoints;
const eukOnly = {...allLabels, labels: ["eukaryota"]};
renderContext.result = EACCore.evaluate(overlapRows, eukOnly);
renderContext.preparePointLayers(geometry, 1, eukOnly);
assert.equal(renderContext.pointLayers.layers[0].canvas, fadedCanvas);
assert.equal(rasterizedPoints - beforeSelection, overlapRows.length);
assert.deepEqual(renderContext.pointLayers.layers[1].canvas.points.map(point => point.color),
  colors.map(color => color === normalPalette[0] ? color : fadedPalette[normalPalette.indexOf(color)]));
for (const point of renderContext.pointLayers.layers[1].canvas.points) assert.equal(point.alpha, .55);
// Dense removed clusters must not accumulate back to the dark retained color.
function compositeOnWhite(points) {
  return points.reduce((pixel, point) => pixel.map((channel, index) =>
    channel * (1 - point.alpha) + parseInt(point.color.slice(1 + index * 2, 3 + index * 2), 16) * point.alpha), [255, 255, 255]);
}
const fadedPixel = compositeOnWhite(fadedCanvas.points);
assert(fadedPixel.every(channel => channel >= 187));
const eukPixel = compositeOnWhite(Array.from({length: 100}, () => ({color: normalPalette[0], alpha: .55})));
const removedEukPixel = compositeOnWhite(Array.from({length: 100}, () => ({color: fadedPalette[0], alpha: .55})));
assert(removedEukPixel.every((channel, index) => channel - eukPixel[index] > 150));
renderContext.preparePointLayers({...geometry, w: 320}, 2, eukOnly);
assert.deepEqual(Array.from(renderContext.pointLayers.order), firstOrder);
renderContext.pointLayers = null;
renderContext.preparePointLayers(geometry, 1, eukOnly);
assert.deepEqual(Array.from(renderContext.pointLayers.order), firstOrder);
renderContext.rows = [];
renderContext.pointLayers = null;
renderContext.result = null;
renderContext.preparePointLayers(geometry, 1, allLabels);
assert.equal(renderContext.pointLayers.layers[0].canvas.points.length, 0);
console.log("Seeded mixed-label order, cache reuse and filtering appearance passed.");

(async()=>{
  const fasta=Buffer.from(">a full header\r\n"+"A".repeat(10001)+"\r\n>b\n"+"C".repeat(500)+"\n>c\n"+"G".repeat(499)+"\n>d tail\r\n"+"T".repeat(935)+"\r\n");
  const coverage="identifier\tmean_reads\nc\t0\nb\t10\nd\t84.9343\na\t100.5\n";
  const ecm="id\tbp\tkingdom\nd\t935\tother\na\t10001\teuk\nc\t499\tNA\nb\t500\tother\n";
  fs.mkdirSync(out,{recursive:true});
  fs.writeFileSync(path.join(out,"assembly.fasta"),fasta);
  fs.writeFileSync(path.join(out,"coverage.tsv"),coverage);
  fs.writeFileSync(path.join(out,"classification.tsv"),ecm);
  const original=new File([fasta],"assembly.fasta");
  const sample={id:"sample_A",fasta:original,coverage:new File([coverage],"coverage.tsv"),ecm:new File([ecm],"classification.tsv")};
  const columns={coverageId:"identifier",coverageValue:"mean_reads",ecmId:"id",ecmLabel:"kingdom",ecmLength:"bp"};
  const prepared=await EACData.prepare(sample,columns);
  const raw=JSON.parse(await prepared.blob.text());
  assert.equal(raw.input_sha256.fasta,crypto.createHash("sha256").update(fasta).digest("hex"));
  const cases=[
    {name:"advanced",min:10,max:null,labels:ruleLabels,retain:true,rules:[advancedRule],count:2,bp:10501},
    {name:"advanced_inclusive",min:10,max:null,labels:ruleLabels,retain:true,rules:[{...advancedRule,length:{operator:">=",value:500}}],count:1,bp:10001},
    {name:"advanced_and",min:10,max:null,labels:ruleLabels,retain:true,rules:[{...advancedRule,coverage:{operator:"<",value:85}}],count:2,bp:10501},
    {name:"advanced_disabled",min:10,max:null,labels:ruleLabels,retain:true,rules:[{...advancedRule,enabled:false}],count:3,bp:11436},
    {name:"advanced_multi",min:10,max:null,labels:ruleLabels,retain:true,rules:[advancedRule,{...advancedRule,id:2,labels:["eukaryota"]}],count:1,bp:500},
    {name:"advanced_overlap",min:10,max:null,labels:ruleLabels,retain:true,rules:[advancedRule,{...advancedRule,id:2}],count:2,bp:10501},
    {name:"advanced_manual",min:10,max:null,labels:ruleLabels,retain:true,rules:[advancedRule],excluded:["a"],count:1,bp:500},
    {name:"advanced_sample",min:null,max:null,labels:ruleLabels,retain:false,rules:[advancedRule],count:0,bp:0},
    {name:"default",min:10,max:null,labels:["eukaryota","other","NA"],retain:true,count:3,bp:11436},
    {name:"fractional",min:85,max:null,labels:["eukaryota","other","NA"],retain:true,count:1,bp:10001},
    {name:"upper",min:10,max:100,labels:["eukaryota","other","NA"],retain:true,count:2,bp:1435},
    {name:"labels",min:10,max:null,labels:["other"],retain:true,count:2,bp:1435},
    {name:"excluded",min:null,max:null,labels:["eukaryota","other","NA"],retain:false,count:0,bp:0},
    {name:"manual",min:10,max:null,labels:["eukaryota","other","NA"],retain:true,excluded:["b"],count:2,bp:10936},
    {name:"manual_batch",min:10,max:null,labels:["eukaryota","other","NA"],retain:true,excluded:["a","c","d"],count:1,bp:500},
    {name:"manual_all",min:10,max:null,labels:["eukaryota","other","NA"],retain:true,excluded:["a","b","c","d"],count:0,bp:0}
  ];
  let manualDecision,advancedDecision;
  for(const item of cases){
    const decision={sample_id:"sample_A",fingerprint:prepared.meta.fingerprint,eac_path:prepared.meta.eac_path,
      eac_sha256:prepared.meta.eac_sha256,reviewed:true,retain_sample:item.retain,min_length_bp:500,min_coverage:item.min,
      max_coverage:item.max,labels:item.labels,excluded_seq_ids:item.excluded||[],advanced_rules:item.rules||[],expected_kept_contigs:item.count,expected_kept_bp:item.bp};
    if(item.name==="manual")manualDecision=decision;
    if(item.name==="advanced_multi")advancedDecision=decision;
    const result=await EACData.apply(prepared.blob,original,decision,prepared.meta);
    fs.writeFileSync(path.join(out,item.name+".fasta"),Buffer.from(await result.fasta.arrayBuffer()));
    fs.writeFileSync(path.join(out,item.name+".tsv"),await result.table.text());
    const zip=await EACZip.build([{name:"sample_A/sample_A.fasta",blob:result.fasta},{name:"sample_A/sample_A.sequences.tsv",blob:result.table}]);
    fs.writeFileSync(path.join(out,item.name+".zip"),Buffer.from(await zip.arrayBuffer()));
    await assert.rejects(()=>EACData.apply(prepared.blob,original,{...decision,reviewed:false},prepared.meta),/Confirm/);
    await assert.rejects(()=>EACData.apply(prepared.blob,original,{...decision,expected_kept_contigs:999},prepared.meta),/counts/);
  }
  await assert.rejects(()=>EACData.apply(prepared.blob,original,{...manualDecision,excluded_seq_ids:["missing"]},prepared.meta),/Unknown manually excluded/);
  await assert.rejects(()=>EACData.apply(prepared.blob,original,{...manualDecision,excluded_seq_ids:["b","b"]},prepared.meta),/Invalid manual/);
  const workerMessages = [];
  let workerHandler;
  const workerContext = {EACCore, EACData, EACZip, sha256, Blob, File,
    importScripts() {}, self: {
      addEventListener(type, handler) {workerHandler=handler;},
      postMessage(message) {workerMessages.push(message);}
    }
  };
  vm.createContext(workerContext);
  vm.runInContext(fs.readFileSync(path.join(root,"web/src/worker.js"),"utf8"),workerContext);
  await workerHandler({data:{id:1,type:"prepare",samples:[sample],columns,projectName:"Manual removal test"}});
  assert(!workerMessages.some(message=>message.type==="error"),JSON.stringify(workerMessages));
  const webManifest=workerMessages.find(message=>message.type==="result").value;
  await workerHandler({data:{id:2,type:"export",decisions:{schema_version:8,dataset_id:webManifest.dataset_id,samples:[manualDecision]}}});
  assert(!workerMessages.some(message=>message.type==="error"),JSON.stringify(workerMessages));
  const workerZip=workerMessages.find(message=>message.id===2&&message.type==="result").value;
  fs.writeFileSync(path.join(out,"worker-manual.zip"),Buffer.from(await workerZip.arrayBuffer()));
  await workerHandler({data:{id:3,type:"export",decisions:{schema_version:9,dataset_id:webManifest.dataset_id,samples:[advancedDecision]}}});
  assert(!workerMessages.some(message=>message.id===3&&message.type==="error"),JSON.stringify(workerMessages));
  const advancedZip=workerMessages.find(message=>message.id===3&&message.type==="result").value;
  fs.writeFileSync(path.join(out,"worker-advanced.zip"),Buffer.from(await advancedZip.arrayBuffer()));
  const missingRules={...advancedDecision};delete missingRules.advanced_rules;
  await workerHandler({data:{id:4,type:"export",decisions:{schema_version:9,dataset_id:webManifest.dataset_id,samples:[missingRules]}}});
  assert(workerMessages.some(message=>message.id===4&&message.type==="error"));

  await assert.rejects(()=>EACData.prepare({...sample,coverage:new File([coverage+"b\t10\n"],"duplicate.tsv")},columns),/duplicate/);
  await assert.rejects(()=>EACData.prepare({...sample,coverage:new File([coverage.replace("b\t10\n","")],"missing.tsv")},columns),/same IDs/);
  await assert.rejects(()=>EACData.prepare({...sample,ecm:new File([ecm.replace("b\t500","b\t501")],"length.tsv")},columns),/length differs/);
  await assert.rejects(()=>EACData.prepare({...sample,coverage:new File([coverage.replace("b\t10","b\tNaN")],"nan.tsv")},columns),/Invalid/);
  const noLength=await EACData.prepare(sample,{...columns,ecmLength:""});
  assert.deepEqual(JSON.parse(await noLength.blob.text()).rows,raw.rows);
  const chunks=new File([fasta],"chunked.fasta");
  chunks.stream=()=>new ReadableStream({start(controller){for(let i=0;i<fasta.length;i+=17)controller.enqueue(new Uint8Array(fasta.subarray(i,i+17)));controller.close();}});
  assert.equal((await EACData.fasta(chunks)).hash,raw.input_sha256.fasta);
  fs.writeFileSync(path.join(out,"cases.json"),JSON.stringify(cases));
  console.log("Browser parsing, boundary filters, exclusions, validation and ZIP generation passed.");
})().catch(error=>{console.error(error);process.exitCode=1;});

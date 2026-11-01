"""Build the browser edition from the shared review template."""
from pathlib import Path
import re
import shutil

root = Path(__file__).resolve().parents[1]
web = root / "web"
dist = web / "dist"
if dist.exists():
    shutil.rmtree(dist)
shutil.copytree(web / "src", dist)
template = (root / "src/eukasmchecker/review.html").read_text()
start = template.index("const C={")
end = template.index("(async () => {", start)
core = template[start:end].replace("const C=", "globalThis.EACCore=", 1)
(dist / "core.js").write_text('"use strict";\n' + core)

needle = ' const meta=metadata.get(name),tag=document.createElement("script");tag.src=meta.data_file;'
replacement = ''' if(PROJECT.browser_local){parent.postMessage({type:"eac:sample",sample:name},browserOrigin);return;}
 const meta=metadata.get(name),tag=document.createElement("script");tag.src=meta.data_file;'''
assert template.count(needle) == 1
template = template.replace(needle, replacement)
needle = 'function readReviewPayload(meta){'
replacement = '''function readReviewPayload(meta){
 if(PROJECT.browser_local)return new Promise((resolve,reject)=>{
  bulkPending={sample_id:meta.sample_id,resolve,reject};
  parent.postMessage({type:"eac:sample",sample:meta.sample_id},browserOrigin);
 }).finally(()=>{bulkPending=null;});'''
assert template.count(needle) == 1
template = template.replace(needle, replacement)

# Browser inputs are selected before review; one action exports the complete archive.
start = template.index('$("save").addEventListener("click",async()=>{')
end = template.index('window.addEventListener("pagehide",saveLocal);', start)
template = template[:start] + '''
$("save").addEventListener("click",async()=>{
 await saveLocal();
 parent.postMessage({type:"eac:export",document:documentData()},browserOrigin);
});
$("webBack").addEventListener("click",()=>{
 saveLocal();
 parent.postMessage({type:"eac:back"},browserOrigin);
});
window.addEventListener("message",event=>{
 if(event.source!==parent||event.origin!==browserOrigin)return;
 if(event.data.type==="eac:payload")window.receiveGenomeSample(event.data.payload);
 if(event.data.type==="eac:language"){lang=event.data.language==="zh"?"zh":"en";translate();}
});
''' + template[end:]
template = re.sub(r'<button id="restore"[^>]*>.*?</button>', '', template)
template = re.sub(r'<input id="restoreFile"[^>]*>', '', template)
template = template.replace('<header><div>', '<header><div class="web-header-left"><button id="webBack" data-i18n="webBack">← Inputs</button>', 1)
template = template.replace('</head>', '''<style>
header h1,.language-toggle,#projectTitle{display:none}
header{top:0;padding:8px 20px}
.web-header-left{display:flex;align-items:center;gap:12px}
#webBack{height:32px;padding:6px 10px;font-size:12px}
.header-meta #saveStatus:not(:empty){border:0;padding:0}
</style></head>''')
template = template.replace('const t=k=>translations[lang][k]', '''
Object.assign(translations.en,{
 webBack:"← Inputs",save:"Save results",saveBackup:"Save results to download your review.",
 allDone:"All samples reviewed. Save results to download them.",
 autoSaved:"Save results before refreshing or closing this page.",
 restoreHelp:" Return to Inputs and choose Import review (.json).",
 ruleOutput:"Save results downloads sequence tables, a summary and decisions. Uploaded FASTA files are filtered into the same archive. Example profiles contain statistics only."
});
Object.assign(translations.zh,{
 webBack:"← 输入文件",save:"保存结果",saveBackup:"保存结果可下载审阅记录。",
 allDone:"全部样本已审阅，可保存结果。",autoSaved:"刷新或关闭页面前请保存结果。",restoreHelp:" 请返回输入页导入审阅记录 (.json)。",
 ruleOutput:"保存结果会下载序列筛选表、汇总和 decisions。上传的 FASTA 会筛选后一并打包。示例仅包含序列统计。"
});
const t=k=>translations[lang][k]''')
# Web preferences and progress live only for the current review.
template = template.replace("localStorage.getItem", "reviewPreferences.get").replace("localStorage.setItem", "reviewPreferences.set")
template = template.replace("(async () => {", "(async () => {\nconst reviewPreferences=new Map();\nconst browserOrigin=new URL(document.baseURI).origin;")
start = template.index("let serverReady=")
end = template.index("function restoreDocument(d){", start)
template = template[:start] + '''function saveLocal(){setStatus("autoSaved");}
function scheduleSave(){saveLocal();}
''' + template[end:]
template = template.replace('await initializeProgress();', '''
if(PROJECT.initial_decisions){
 try{restoreDocument(PROJECT.initial_decisions);await saveLocal();setStatus("restored");}
 catch(e){parent.postMessage({type:"eac:restore-error",message:e.message},browserOrigin);return;}
}
parent.postMessage({type:"eac:ready"},browserOrigin);''')
template = template.replace('selectSample(PROJECT.samples.find(m=>!settings.get(m.sample_id).reviewed)?.sample_id||PROJECT.samples[0].sample_id);',
                            'selectSample(PROJECT.initial_sample||PROJECT.samples.find(m=>!settings.get(m.sample_id).reviewed)?.sample_id||PROJECT.samples[0].sample_id);')
(dist / "review-template.html").write_text(template)
print(dist)
